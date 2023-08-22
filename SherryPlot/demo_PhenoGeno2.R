library(NetBID2)
source('D:/analysis_eng/SherryPlot/function_circos.R')
source('D:/analysis_eng/SherryPlot/function_basic.R')
## option
load('D:/analysis_eng/SherryPlot/data/WES_gene_bed.RData')

##
dis2gene_file <- 'data/dis2gene-HYPERBILIRUBINEMIA20210917.xlsx';
out_pdf <- 'plot/HYPERBILIRUBINEMIA20210917.pdf'
#
dis2gene <- clean_table(read.xlsx(dis2gene_file))
pdf(out_pdf,width=14,height=12)
par(mar=c(4,4,4,15))
draw_circos_phenoGeno(dis2gene,gene_bed=WES_gene,class1='',class2='',
                      gene_cex=0.8,num_cex=0.8,radius_chr=1.2,radius_system=1.25,
                      radius_gene=1,legend_col='Main.category',count_col='patient.number')
dev.off()

