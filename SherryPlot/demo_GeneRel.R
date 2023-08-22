source('D:/analysis_eng/SherryPlot/function_circos.R')
load('data/WES_gene_bed.RData')
###
gene_rel_file <- 'data/geneRel.xlsx'
count_thre <- 1
####
pre_define <- c('blue', 'red', 'yellow', 'green','yellow', 'green')
names(pre_define) <- c('WNT', 'SHH', 'Group3', 'Group4','GroupC', 'GroupD')
gene_rel <- read.xlsx(gene_rel_file)
pdf('plot/geneRel.pdf',width=8,height = 8)
par(mar=c(2,2,2,2))
draw_circos_genoRel(gene_rel,pre_define,count_thre,WES_gene)
dev.off()

