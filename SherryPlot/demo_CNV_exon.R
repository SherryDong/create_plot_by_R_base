library(quantsmooth)
library(RColorBrewer)

### plot for exon-level
draw_CNV_exon <- function(normalize_read_path,normalize_read_dat=NULL,
                          sample.names,
                          or_thre = c(0.8,1.2),min_mean_count=10,
                          cytoband_width=0.15,plot_add=0.2,plot_max_ratio=1,
                          point_cex=0.2,
                          cex_chr=1,draw_bandname=FALSE,cex_bandname=0.3,
                          dup_del_col=brewer.pal(9,'Set1')[1:2]){
  if(is.null(normalize_read_dat)==TRUE){
    normalize_read_dat <- read.delim(normalize_read_path,stringsAsFactors = F)
  }
  rownames(normalize_read_dat) <- normalize_read_dat[,4]
  counts_mat <- normalize_read_dat[,5:ncol(normalize_read_dat)]
  gene_info <- normalize_read_dat[,1:4]
  sample_file <- read.delim(normalize_read_path,
                            stringsAsFactors = F,header = F,nrow=1)[5:ncol(normalize_read_dat)]
  sample_idx<-which(sample_file %in% sample.names) # 210826-wx
  if(length(sample_idx) != 1){
    print(sprintf('%s not in %s, please check and retry !',sample_idx,normalize_read_path));return()
  }
  ## normalize
  mean_counts_mat <- apply(counts_mat,1,mean)
  w1 <- which(mean_counts_mat>=min_mean_count) ## remove
  counts_mat <- counts_mat[w1,];gene_info <- gene_info[w1,]
  norm_counts_mat <- t(apply(counts_mat,1,function(x){x/mean(x,na.rm=T)}))
  target_counts <- norm_counts_mat[,sample_idx]
  bg_counts <- norm_counts_mat[,-sample_idx]
  bg_counts_mean <- apply(bg_counts,1,mean)
  or <- target_counts/bg_counts_mean
  # filter by or
  w1 <- which(or<=or_thre[1] | or>=or_thre[2])
  if(length(w1) == 0){
    print('check or_thre and retry!');return()
  }
  or <- or[w1]; gene_info <- gene_info[w1,,drop=F]
  ## draw cytoband
  pos_align <- 'horizon'
  if(pos_align=='horizon'){
    CHR<-c(1:22,'X','Y')  # Chromosomes
    chrlen<-lengthChromosome(CHR,"bases")
    chr_part1 <- c(1:12);
    chr_part2 <- c(22:13,'X','Y');
    chr_part1 <- as.character(chr_part1)
    chr_part2 <- as.character(chr_part2)
    chr_part  <- cbind(chr_part1,chr_part2)
    max_w <- max(apply(chr_part,1,function(x)chrlen[x[1]]+chrlen[x[2]]))*1.05
    max_h <- max(length(chr_part1),length(chr_part2))
    #
    par(mar=c(4,4,4,4))
    plot(1,col='white',xlim=c(0,max_w),ylim=c(0,max_h),xaxt='n',yaxt='n',xaxs='i',yaxs='i',bty='n',xlab='',ylab='')
    for(i in 1:length(chr_part1)){
      paintCytobands(chr_part1[i],pos=c(0,max_h-i+1),width=cytoband_width,xpd=TRUE,
                     legend=draw_bandname,cex.leg = cex_bandname)
      text(0,max_h-i+1-cytoband_width/2,chr_part1[i],cex=cex_chr,xpd=T,pos=2)
    }
    for(i in 1:length(chr_part2)){
      paintCytobands(chr_part2[i],pos=c(max_w-chrlen[chr_part2[i]],max_h-i+1),width=cytoband_width,xpd=TRUE,
                     legend=draw_bandname,cex.leg = cex_bandname)
      text(max_w,max_h-i+1-cytoband_width/2,chr_part2[i],cex=cex_chr,xpd=T,pos=4)
    }
  }
  ## color
  or_log <- log2(or+1)-1;
  bb <- cut(or_log,breaks = c(log2(1+c(0,0.25,0.6,1/0.6,1/0.25))-1,Inf),include.lowest = T,right=F) ## 
  cc <- colorRampPalette(c(dup_del_col[2],'light grey',dup_del_col[1]))(length(levels(bb)))
  cc <- adjustcolor(cc,alpha.f = 0.8);names(cc) <- levels(bb);
  bb_max <- c(-2,-1,0,1,2)/(3*plot_max_ratio); names(bb_max) <- levels(bb)
  ## draw exon
  gene_info$mid <- gene_info$Start/2+gene_info$End/2
  gene_info$noChr <- gsub('chr','',gene_info$Chr)
  for(i in 1:length(chr_part1)){
    w1 <- which(gene_info$noChr==chr_part1[i])
    gene_use <- gene_info[w1,,drop=F]; bb_use <- bb[w1]
    #segments(x0=gene_use$mid,x1=gene_use$mid,y0=max_h-i+1-cytoband_width-plot_add,
    #         y1=max_h-i+1-cytoband_width-plot_add+plot_add*bb_max[bb_use],
    #         col=cc[bb_use],lwd=point_cex,xpd=T)
    points(x=gene_use$mid,y=max_h-i+1-cytoband_width-plot_add+plot_add*bb_max[bb_use],
           col=cc[bb_use],cex=point_cex,pch=16,xpd=T)
  }
  for(i in 1:length(chr_part2)){
    w1 <- which(gene_info$noChr==chr_part2[i])
    if(length(w1)==0) next
    gene_use <- gene_info[w1,,drop=F]; bb_use <- bb[w1];
    #segments(x0=max_w-chrlen[chr_part2[i]]+gene_use$mid,
    #         x1=max_w-chrlen[chr_part2[i]]+gene_use$mid,y0=max_h-i+1-cytoband_width-plot_add,
    #         y1=max_h-i+1-cytoband_width-plot_add+plot_add*bb_max[bb_use],
    #         col=cc[bb_use],lwd=point_cex,xpd=T)
    points(x=max_w-chrlen[chr_part2[i]]+gene_use$mid,
           y=max_h-i+1-cytoband_width-plot_add+plot_add*bb_max[bb_use],
           col=cc[bb_use],cex=point_cex,pch=16,xpd=T)
  }
  return(NULL)
}
# case
normalize_read_path <- 'data/AH7JC2DSX2_BH7JC5DSX2_autochr_normalized.reads.txt'
case_id <- '21Y74206'
normalize_read_dat <- read.delim(normalize_read_path,stringsAsFactors = F)
pdf('CNV_exon_21Y74206.pdf',width=12,height=8)
draw_CNV_exon(normalize_read_path=normalize_read_path,
              normalize_read_dat=normalize_read_dat,sample.names=case_id,
              or_thre=c(0.8,1/0.8),
              plot_max_ratio=2)
dev.off()
pdf('CNV_exon_case21Y74278.pdf',width=12,height=8)
draw_CNV_exon(normalize_read_path=normalize_read_path,
              normalize_read_dat=normalize_read_dat,sample.names='21Y74278',
              or_thre=c(0.8,1/0.8),
              plot_max_ratio=2)
dev.off()



