
source('D:/analysis_eng/SherryPlot/function_circos.R')
## option
load('data/WES_gene_bed.RData')
dis2gene_file <- 'data/dis2gene.xlsx'
##
dis2gene <- read.xlsx(dis2gene_file)
pdf('plot/phenoGeno.pdf',width=11,height=10)
par(mar=c(4,4,4,10))
draw_circos_phenoGeno(dis2gene,gene_bed=WES)
dev.off()







