draw_brainspan_pattern <- function(use_gene,use_mat,phe_info,choose='activity',
                                   use_phe=c('class1_mod'),main='',
                                   use_boxplot=FALSE){
  use_gene <- intersect(use_gene,rownames(use_mat))
  pp <- get_obs_label(phe_info,use_col = use_phe)
  pp <- factor(pp,levels=unique(pp))
  draw_mat <- use_mat[use_gene,,drop=F]
  tmp1 <- aggregate(t(draw_mat),list(pp),c)
  tmp2 <- aggregate(t(draw_mat),list(pp),median)
  draw_mat_1 <- tmp2[,-1,drop=F]; pp1 <- tmp2$Group.1
  cc <- get.class.color(pp)
  c1 <- get.class.color(use_gene)
  if(use_boxplot==TRUE){
    par(mar=c(8,3,4,2))
    pp2 <- levels(pp1)
    boxplot(t(draw_mat_1),col=cc[as.character(pp1)],names=NA,outline=F,main=main)
    axis(side=1,at=1:length(pp2),labels=NA)
    par1 <- par()$usr
    text(1:length(pp2),par1[3]-(par1[4]-par1[3])/20,pp2,xpd=T,adj=1,srt=60,cex=0.9)
  }else{
    par(mar=c(8,3,4,10))
    plot(y=draw_mat_1[,1],x=1:length(pp1),
         pch=16,ylab=choose,xlab='',
         cex.lab=1.2,main=main,
         col=c1[1],xaxt='n',cex=1,
         ylim=c(min(draw_mat_1),max(draw_mat_1)))
    for(g1 in use_gene){
      points(y=draw_mat_1[,g1],x=1:length(pp1),col=c1[g1],pch=16)
      lines(y=draw_mat_1[,g1],x=1:length(pp1),col=c1[g1])
    }
    axis(side=1,at=1:length(pp1),labels=NA)
    par1 <- par()$usr
    text(1:length(pp1),par1[3]-(par1[4]-par1[3])/20,pp1,xpd=T,adj=1,srt=60,cex=0.9)
    #
    legend(par1[2],par1[4],fill=c1,legend=use_gene,
           border=NA,bty='n',ncol=1,xpd=T,xjust=0,cex=0.7)
    
  }
}