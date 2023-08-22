library(NetBID2)
source('D:/analysis_eng/SherryPlot/function_circos.R')
##
pheno_rel_file  <- 'data/PheoRel.xlsx'
pv_thre <- 0.01

##
d1  <- read.xlsx(pheno_rel_file)
tab <- colSums(d1,na.rm=T)
tab <- tab[which(tab>0)]
d3  <- d1[,names(tab)]
rownames(d3) <- sprintf('S%s',1:nrow(d3))
ov <- matrix(0,nrow=length(tab),ncol=length(tab))
colnames(ov) <- rownames(ov) <- names(tab)
for(i in 1:ncol(d3)){
  for(j in 1:ncol(d3)){
    if(i==j) next 
    ov_count <- apply(d3[,c(i,j)],1,function(x)sum(x,na.rm=T))
    w1 <- which(ov_count==2)
    t1 <- table(list(d3[,i],d3[,j]))
    pv <- fisher.test(t1,alternative = 'greater')$p.value
    if(pv<pv_thre) ov[i,j] <- length(w1)
  }
}
## if 中文
library(Cairo)
CairoPDF('plot/phenoRel.pdf',family='SimSun',width=8,height=8)
draw_circos_complicate(tab,ov_mat=ov,ori_mat=d3,use_cairo=TRUE,
                       fig_lim=1.4,ov_text=F)
dev.off()


