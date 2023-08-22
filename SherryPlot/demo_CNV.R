library(NetBID2)
use_cnv_data <- read.xlsx('data/use_cnv_data.xlsx')
source('functions_CNV.R')
#######################################
pdf('plot/CNV.pdf',width=12,height=8)
draw_CNV(data=use_cnv_data,Type_col='Type',Sample_col='Sample',Chr_col='chromosome',Start_col='start_position',
         End_col='end_position',pos_align='verticle',
         cytoband_width=0.1,each_width=0.015,cex_chr=1,draw_bandname=FALSE,cex_bandname=0.3,
         dup_del_col=brewer.pal(9,'Set1')[2:1])
draw_CNV(data=use_cnv_data,Type_col='Type',Sample_col='Sample',Chr_col='chromosome',Start_col='start_position',
         End_col='end_position',pos_align='horizon',
         cytoband_width=0.1,each_width=0.015,cex_chr=1,draw_bandname=FALSE,cex_bandname=0.3,
         dup_del_col=brewer.pal(9,'Set1')[2:1])
dev.off()








