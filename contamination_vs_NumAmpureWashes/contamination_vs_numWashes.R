rm(list=ls())
require(data.table)
require(scales)
require(cowplot)
require(lattice)
require(dplyr)
require(tidyr)
require(Matrix)

setwd('~/singlecell-pipeline/mHSC/subtype_xContamination_analysis/contamination_vs_NumAmpureWashes/')

# munge 2-wash data ----
load('data/fetal_venule_exp1_HS.ReadsPerGene.matrix.RData')
mat.HiSeq.2r = copy(mat)
rm(mat)
load('data/fetal_venule_exp1_NS.ReadsPerGene.matrix.RData')
mat.NextSeq.2r = copy(mat)
rm(mat)

# get basic stats (not sure if I will need)
cell.info.NS = data.table(cell.name = rownames(mat.NextSeq.2r),
                          numreads = rowSums(mat.NextSeq.2r),
                          numgenes = rowSums(mat.NextSeq.2r>0))
cell.info.NS[, seq.well := substr(cell.name, regexpr("_", cell.name)+1, nchar(cell.name))]
cell.info.NS[, row := substr(seq.well, 1, 1)]
cell.info.NS[, col := substr(seq.well, 2, 3)]

# add sequencing plate info
fv.p1.info=fread('data/fetal_venule_exp1_plateLayout1.csv', header = T)
fv.p1.info=data.table(gather(fv.p1.info, column, cell.name.dil, -row))
fv.p1.info[, cell.name :=trimws(substr(cell.name.dil, 1, regexpr("\\s", cell.name.dil)))]
fv.p1.info[, column := as.numeric(as.character(column))]
fv.p1.info[, cell.name := paste0(cell.name, "_", row, sprintf("%02d", column))]

fv.p2.info=fread('data/fetal_venule_exp1_plateLayout2.csv', header = T)
fv.p2.info=data.table(gather(fv.p2.info, column, cell.name.dil, -row))
fv.p2.info[, cell.name :=trimws(substr(cell.name.dil, 1, regexpr("\\s", cell.name.dil)))]
fv.p2.info[, column := as.numeric(as.character(column))]
fv.p2.info[, cell.name := paste0(cell.name, "_", row, sprintf("%02d", column))]

cell.info.NS <- cell.info.NS[cell.name %in% c(fv.p1.info[, cell.name], fv.p2.info[, cell.name])]
cell.info.NS[cell.name %in% fv.p1.info[, cell.name], seq.plate := 'p1']
cell.info.NS[cell.name %in% fv.p2.info[, cell.name], seq.plate := 'p2']

# Quantify swapping in 2-wash data ------

# function to get nth non-zero max
max.n <- function(x, n){
  return(sort(x, decreasing = T)[n])
}
# function to get diff between 1st and n+1th max
diff.n <- function(x, n){
  return(max.n(x, 1) - max.n(x, n+1))
}
mat.p1.NS <- mat.NextSeq.2r[cell.info.NS[seq.plate=='p1', cell.name], ]
mat.p1.NS <- mat.p1.NS[, colSums(mat.p1.NS) > 0]
mat.p1.NS <- log10(mat.p1.NS + 1)

mat.p1.HS <- mat.HiSeq.2r[cell.info.NS[seq.plate=='p1', cell.name], ]
mat.p1.HS <- mat.p1.HS[, colSums(mat.p1.HS) > 0]
mat.p1.HS <- log10(mat.p1.HS + 1)

# Use NextSeq data to find low-expressed genes
stats.p1.ns <- data.table(gene = colnames(mat.p1.NS),
                          numcells = colSums(mat.p1.NS > 0),
                          diff.10 = apply(mat.p1.NS, 2, diff.n, n=10))

stats.p1.ns[order(diff.10), ]
g.plot = stats.p1.ns[order(diff.10, decreasing = T)][!gene %like% '^N_', gene][1]

# Plot expression in plate, NextSeq data
dt.gene.NS <- dcast.data.table(data.table(inner_join(data.table(cell.name = rownames(mat.p1.NS), count=mat.p1.NS[, g.plot]),
                      cell.info.NS[, .(cell.name, row, col)])), row ~ col, value.var = 'count')
plate.expr.NS <- as.matrix(dt.gene.NS[, 2:ncol(dt.gene.NS)])
rownames(plate.expr.NS) <- dt.gene.NS[, row]
color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

pdf(paste0("plots/plateMatrix_",g.plot,"_counts_2wash_NextSeq.pdf"),width = 8,height = 6)
levelplot(t(plate.expr.NS),col.regions=color.palette(100),scales=list(x=list(rot=45)),
          main=paste0('Expression of ',g.plot,' (log10 cpm) 2-wash NextSeq plate'))
dev.off()

# Plot expression in plate, HiSeq data
dt.gene.HS <- dcast.data.table(data.table(inner_join(data.table(cell.name = rownames(mat.p1.HS), count=mat.p1.HS[, g.plot]),
                      cell.info.NS[, .(cell.name, row, col)])), row ~ col, value.var = 'count')
plate.expr.HS <- as.matrix(dt.gene.HS[, 2:ncol(dt.gene.HS)])
rownames(plate.expr.HS) <- dt.gene.HS[, row]
color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

pdf(paste0("plots/plateMatrix_",g.plot,"_counts_2wash_HiSeq.pdf"),width = 8,height = 6)
levelplot(t(plate.expr.HS),col.regions=color.palette(100),scales=list(x=list(rot=45)),
          main=paste0('Expression of ',g.plot,' (log10 cpm) in 2-wash HiSeq plate'))
dev.off()

# Plot false positive expression level w/in row J and column 19
row.use = 'J'
col.use = 19
cells.falsepos <- union(cell.info.NS[row==row.use, cell.name],cell.info.NS[col==col.use, cell.name])
cells.falsepos <- cells.falsepos[cells.falsepos %in% rownames(mat.p1.NS)[mat.p1.NS[, g.plot] == 0]]
max.expr <- max(max(mat.p1.HS[, g.plot]))

hist(mat.p1.HS[cells.falsepos, g.plot] - max.expr)

falsepos.2washes <- mat.p1.HS[cells.falsepos, g.plot] - max.expr


# Quantify swapping rate in 1-wash library -------
load('data/mHSC_plate1_HS.counts_matrix.RData')
load('data/mHSC_plate1_NS.counts_matrix.RData')

rm(cell.info.NS)
cell.info.NS = data.table(cell.name = rownames(mat.mHSC_p1_NS),
                          numreads = rowSums(mat.mHSC_p1_NS),
                          numgenes = rowSums(mat.mHSC_p1_NS>0))
cell.info.NS[, seq.well := substr(cell.name, nchar(cell.name) - 2, nchar(cell.name))]
cell.info.NS[, row := substr(seq.well, 1, 1)]
cell.info.NS[, col := substr(seq.well, 2, 3)]

max.n <- function(x, n){
  return(sort(x, decreasing = T)[n])
}
# function to get diff between 1st and n+1th max
diff.n <- function(x, n){
  return(max.n(x, 1) - max.n(x, n+1))
}
rm(mat.p1.NS,mat.p1.HS)
mat.p1.NS <- mat.mHSC_p1_NS[, colSums(mat.mHSC_p1_NS) > 0]
mat.p1.NS <- log10(mat.p1.NS + 1)
mat.p1.HS <- mat.mHSC_p1_HS[, colSums(mat.mHSC_p1_HS) > 0]
mat.p1.HS <- log10(mat.p1.HS + 1)

# Use NextSeq data to find low-expressed genes
stats.p1.ns <- data.table(gene = colnames(mat.p1.NS),
                          numcells = colSums(mat.p1.NS > 0),
                          max.1 = apply(mat.p1.NS, 2, max.n, n=1),
                          diff.10 = apply(mat.p1.NS, 2, diff.n, n=10))

g.plot = stats.p1.ns[order(diff.10, decreasing = T)][!gene %like% '^N_', gene][1]

# Plot expression in plate, NextSeq data
dt.gene.NS <- dcast.data.table(data.table(inner_join(data.table(cell.name = rownames(mat.p1.NS), count=mat.p1.NS[, g.plot]),
                      cell.info.NS[, .(cell.name, row, col)])), row ~ col, value.var = 'count')
plate.expr.NS <- as.matrix(dt.gene.NS[, 2:ncol(dt.gene.NS)])
rownames(plate.expr.NS) <- dt.gene.NS[, row]
color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

pdf(paste0("plots/plateMatrix_",g.plot,"_counts_1wash_NextSeq.pdf"),width = 8,height = 6)
levelplot(t(plate.expr.NS),col.regions=color.palette(100),scales=list(x=list(rot=45)),
          main=paste0('Expression of ',g.plot,' (log10 cpm) 2-wash NextSeq plate'))
dev.off()

# Plot expression in plate, HiSeq data
dt.gene.HS <- dcast.data.table(data.table(inner_join(data.table(cell.name = rownames(mat.p1.HS), count=mat.p1.HS[, g.plot]),
                      cell.info.NS[, .(cell.name, row, col)])), row ~ col, value.var = 'count')
plate.expr.HS <- as.matrix(dt.gene.HS[, 2:ncol(dt.gene.HS)])
rownames(plate.expr.HS) <- dt.gene.HS[, row]
color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")

pdf(paste0("plots/plateMatrix_",g.plot,"_counts_1wash_HiSeq.pdf"),width = 8,height = 6)
levelplot(t(plate.expr.HS),col.regions=color.palette(100),scales=list(x=list(rot=45)),
          main=paste0('Expression of ',g.plot,' (log10 cpm) in 1-wash HiSeq plate'))
dev.off()

# Plot false positive expression level w/in row B and column 05 (column w/high expressor)
row.use = 'B'
col.use = '05'
cells.falsepos <- union(cell.info.NS[row==row.use, cell.name],cell.info.NS[col==col.use, cell.name])
cells.falsepos <- cells.falsepos[cells.falsepos %in% rownames(mat.p1.NS)[mat.p1.NS[, g.plot] == 0]]
max.expr <- max(max(mat.p1.HS[, g.plot]))

falsepos.1washes <- mat.p1.HS[cells.falsepos, g.plot] - max.expr

# 1 vs 2 washes swapping rate ------
swap_rates <- rbindlist(list(data.table(cell.name = names(falsepos.1washes),
                         log_rate = falsepos.1washes),
                         data.table(cell.name = names(falsepos.2washes),
                         log_rate = falsepos.2washes)),
                        use.names = T, fill = F)
swap_rates[cell.name %in% names(falsepos.1washes), num.washes := "1"]
swap_rates[cell.name %in% names(falsepos.2washes), num.washes := "2"]
swap_rates[, rate := 10^log_rate]

# Swapping rate between samples
swap_meds <- swap_rates[, median(rate), by = num.washes]
p=ggplot(swap_rates, aes(num.washes, rate, fill = num.washes)) +
  geom_boxplot() + 
  ylab('Sample-sample swap rate') +
  theme(legend.position = 'none', panel.grid.major.y = element_line(color='grey80')) +
  xlab('AMPure washes') +
  scale_y_continuous(trans = log10_trans(),
                      breaks = trans_breaks("log10", function(x) 10^x, n=3),
                      labels = trans_format("log10", math_format(10^.x))) +
  geom_text(data=swap_meds, mapping=aes(num.washes,0.03, label = sprintf("%.5f%%", V1)),
            size = 6)
p
save_plot('plots/swapRate_byWashes_sampleToSample.pdf',p, base_height = 3.5)


swap_rates[, total_misassigned := ( rate*(15+23) / (rate*(15+23) + 1) )]
swap_meds <- swap_rates[, median(total_misassigned), by = num.washes]
p=ggplot(swap_rates, aes(num.washes, total_misassigned, fill = num.washes)) +
  geom_boxplot() + 
  ylab('Total misassigned rate') +
  theme(legend.position = 'none', panel.grid.major.y = element_line(color='grey80')) +
  xlab('AMPure washes') +
  scale_y_continuous(trans = log10_trans(),
                      breaks = trans_breaks("log10", function(x) 10^x, n = 5),
                      labels = trans_format("log10", math_format(10^.x))) +
  geom_text(data=swap_meds, mapping=aes(num.washes,0.8, label = sprintf("%.4f%%", V1)),
            size = 6)
p
save_plot('plots/swapRate_byWashes_totalMisassigned.pdf',p, base_height = 3.5)
