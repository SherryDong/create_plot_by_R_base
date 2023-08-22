library(NetBID2)
source('D:/analysis_eng/SherryPlot/functions_pheno.R')
##
table_file <- 'data/table_diagPheno.xlsx'
##
d1 <- read.xlsx(table_file)
d1$System_Undiagnosed <- d1$System-d1$System_Diagnosed
p1 <- paste0(round(10000*d1$Diagnosed/d1$Total)/100,'%')
pdf('plot/diagpheno_f1.pdf',width=10,height=10)
draw_phenoDist(d1,Name_col='Label',Total_col='System',
               Sub_col=c('System_Diagnosed','System_Undiagnosed'),
               cex.label=1.1,cex.legend = 1.6,color_name = 'Set3')
dev.off()  
#########
total <- d1$System[1]+d1$UnSystem[1]
d1$Diagnosed <- d1$System_Diagnosed+d1$UnSystem_Diagnosed
d1$System_UnDiagnosed <- d1$System-d1$System_Diagnosed
d1$system_per <- d1$System/total
d1$system_diagnose_per <- d1$System_Diagnosed/d1$Diagnosed
d1$system_undiagnose_per <- d1$System_UnDiagnosed/(total-d1$Diagnosed)
d1$OR <- d1$system_diagnose_per/d1$system_undiagnose_per
##
pdf('plot/diagpheno_f2.pdf',width=10,height=10)
par(mar=c(8,8,8,8))
draw_phenoDist_Two(d1,Name_col='Label',Total_num=4214,Total_col='System',
                   Sub_num=294,Sub_col='System_Diagnosed',
                   cex.label=1,cex.legend = 1.1,color_name = 'Set1',
                   Total_legend_text='Total Patients',Sub_legend_text='Diagnosed',
                   legend_text='Patients with involvement\nof this organ')
dev.off()

#########









