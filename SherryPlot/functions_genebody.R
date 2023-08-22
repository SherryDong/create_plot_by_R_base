# Q1: CFTR gene structure
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggbio)
library(Homo.sapiens)
library(biovizBase)
library(VariantAnnotation)
class(Homo.sapiens)
library(EnsDb.Hsapiens.v75)
ensdb <- EnsDb.Hsapiens.v75
#library(EnsDb.Hsapiens.v86)
#ensdb <- EnsDb.Hsapiens.v86 ##hg38
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
data(genesymbol, package = "biovizBase")
library(RColorBrewer)
library(readxl);library(XML);library(NetBID2)
###
get_pfam <- function(x){
  ori_pf1 <- readHTMLTable(sprintf('http://pfam.xfam.org/protein/%s',x),stringsAsFactors=F)
  pf1 <- ori_pf1$imageKey[,1:4];colnames(pf1) <- c('Source','Domain','Start','End')
  pf1$Domain[which(pf1$Domain=='n/a')] <- pf1$Source[which(pf1$Domain=='n/a')]
  pf1
}
find_pos2exon <- function(exon_tab,x){
  w1 <- which(exon_tab$start-x<=0)
  w2 <- which(exon_tab$end-x>=0)
  w3 <- intersect(w1,w2)
  if(length(w3)>0) return(exon_tab[w3[1],'exon_id'])
  w1 <- min(which(exon_tab$start-x>=0))-1
  w2 <- max(which(exon_tab$end-x<=0))+1
  return(sprintf('%s-%s',exon_tab[w2,'exon_id'],exon_tab[w2,'exon_id']))
}
find_pos2intron <- function(intron_tab,x){
  w1 <- which(intron_tab$start-x<=0)
  w2 <- which(intron_tab$end-x>=0)
  w3 <- intersect(w1,w2)
  return(intron_tab[w3,])
}
##
type2color <- c('Missense mutation',
                "Splicing mutation",
                "Stopgain mutation",
                "Synonymous mutation",
                "Small insertion/deletion",
                "Nonframeshift substitution","Frameshift mutation",
                "Startlost SNV",'Other')
cc1 <- colorRampPalette(brewer.pal(12,'Paired'))(length(type2color));
cc1 <- adjustcolor(cc1,alpha=0.8);names(cc1) <- type2color;
type2color <- cc1
transfer_mutation_type <- function(x){
  use_type <- c('Missense mutation',
                "Splicing mutation",
                "Stopgain mutation",
                "Synonymous mutation",
                "Small insertion/deletion",
                "Nonframeshift substitution","Frameshift mutation",
                "Startlost SNV",'Other')
  x[which(x %in% 'nonsynonymous SNV')] <- 'Missense mutation'
  x[which(x %in% c('splice_region_variant','splicing'))] <- 'Splicing mutation'
  x[which(x %in% c('stopgain SNV','stopgain'))] <- 'Stopgain mutation'
  x[which(x %in% c('frameshift substitution'))] <- 'Frameshift mutation'
  x[which(x %in% c('nonframeshift substitution'))] <- 'Nonframeshift substitution'
  x[which(x %in% c('synonymous SNV'))] <- 'Synonymous mutation' #
  x[which(x %in% c('startloss'))] <- 'Startlost SNV' #
  x[grep('indel',x)] <- 'Small insertion/deletion' #
  x[grep('insertion',x)] <- 'Small insertion/deletion' #
  x[grep('deletion',x)] <- 'Small insertion/deletion' #
  x[which(!x %in% use_type)] <- 'Other'
  x
}

##
draw_gene <- function(use_gene,use_ensg,use_enst,all_variant=NULL,count_log=FALSE,
                      log_use=2,mark_top_num=3,
                      mark_top_cex=0.3,use_domain=NULL,
                      mark_common=TRUE,
                      only_use_trascript=FALSE,mark_lwd=0.5,
                      step_h=2,max_count_val=NULL,use_srt=45,
                      draw_chr=TRUE,draw_ENSG=TRUE,
                      draw_ENST=TRUE,draw_mRNA=TRUE,
                      draw_CDS=TRUE,draw_domain=TRUE){
  transcript_tab <- as.data.frame(transcripts(ensdb,filter = ~ gene_id == use_ensg &
                                                tx_biotype == "protein_coding"))
  exon_tab <- as.data.frame(exons(ensdb,filter = ~ gene_id == use_ensg &
                                    tx_biotype == "protein_coding"))
  intron_tab <- as.data.frame(intronsByTranscript(ensdb,filter = ~ gene_id == use_ensg &
                                                    tx_biotype == "protein_coding"))
  cds_tab <- as.data.frame(cdsBy(ensdb,filter=~ gene_id == use_ensg &
                                   tx_biotype == "protein_coding"))
  domain2color <- c(); mutationType2color <- c()
  ## draw exon structure for a gene
  if(only_use_trascript==TRUE){
    transcript_tab <- transcript_tab[which(transcript_tab$tx_name==use_enst),]
  }
  strand <- as.character(transcript_tab[1,'strand'])
  range_min <- min(exon_tab$start);
  range_max <- max(exon_tab$end);
  dd <- (range_max-range_min)/20;
  n_transcript <- nrow(transcript_tab)
  plot(1,xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i',bty='n',col='white',
       xlim=c(range_min-dd,range_max+dd*3),
       ylim=c(-6-length(all_variant)*step_h,n_transcript/2+1))
  pp <- par()$usr ## merge exon
  # draw exon
  ss <- seq(range_min,range_max,length.out=11);
  rr <- floor(min(log10(ss)))-5
  ss1 <- round(ss/(10^rr))*10^rr
  if(draw_chr==TRUE) axis(side=3,at=ss,labels=get_label_manual(ss1),xpd=TRUE,cex.axis=0.45)
  #segments(x0=range_min,x1=range_max,y0=n_transcript+0.5,y1=n_transcript+0.5,lty=2,lwd=0.5,col='grey')
  if(draw_chr==TRUE) text(range_min-dd/4,(n_transcript+2)/2,sprintf('Chr%s',transcript_tab[1,1]),xpd=TRUE,cex=0.6,adj=1)
  #text(range_min-dd/4,(n_transcript+1)/2,use_ensg,xpd=TRUE,cex=0.8,adj=1)
  if(draw_ENSG==TRUE) rect(xleft=exon_tab$start,xright=exon_tab$end,ybottom=n_transcript/2,ytop=(n_transcript+1)/2,col='black',border = 'black')
  if(draw_ENSG==TRUE){
  for(i in 1:nrow(intron_tab)){
    segments(x0=intron_tab[i,'start'],x1=intron_tab[i,'end'],y0=(n_transcript+0.5)/2,y1=(n_transcript+0.5)/2)
    #segments(x0=intron_tab[i,'start'],x1=intron_tab[i,'end']/2+intron_tab[i,'start']/2,y0=n_transcript+0.5,y1=n_transcript+1)
    #segments(x1=intron_tab[i,'end'],x0=intron_tab[i,'end']/2+intron_tab[i,'start']/2,y0=n_transcript+1,y1=n_transcript+0.5)
  }}
  #
  all_enst <- unique(transcript_tab$tx_id)
  all_enst <- c(use_enst,setdiff(all_enst,use_enst))
  if(draw_ENST==TRUE){
  for(i in 1:length(all_enst)){
    each_enst <- all_enst[i]
    if(each_enst==use_enst){
      text(range_min-dd/4,(i-0.75)/2,each_enst,adj=1,xpd=TRUE,cex=0.6,col='blue')
    }else{
      text(range_min-dd/4,(i-0.75)/2,each_enst,adj=1,xpd=TRUE,cex=0.5)
    }
    transcript_tab_1 <- as.data.frame(exons(ensdb,filter = ~ tx_id == each_enst &
                                              tx_biotype == "protein_coding"))
    intron_tab_1 <- intron_tab[which(intron_tab$group_name==each_enst),]
    for(j in 1:nrow(transcript_tab_1)){
      rect(xleft=transcript_tab_1[j,'start'],xright=transcript_tab_1[j,'end'],ybottom=(i-0.5)/2,ytop=(i-1)/2,col='black',border = 'black')
    }
    for(j in 1:nrow(intron_tab_1)){
      segments(x0=intron_tab_1[j,'start'],x1=intron_tab_1[j,'end']/2+intron_tab_1[j,'start']/2,y0=(i-1)/2,y1=(i-0.5)/2,col='dark grey')
      segments(x1=intron_tab_1[j,'end'],x0=intron_tab_1[j,'end']/2+intron_tab_1[j,'start']/2,y0=(i-0.5)/2,y1=(i-1)/2,col='dark grey')
    }
  }
  }
  # merge exon, need to specify transcript
  exon_tab <- as.data.frame(exons(ensdb,filter = ~ tx_id == use_enst &
                                    tx_biotype == "protein_coding"))
  intron_tab <- as.data.frame(intronsByTranscript(ensdb,filter = ~ tx_id == use_enst &
                                                    tx_biotype == "protein_coding"))
  cds_tab <- as.data.frame(cdsBy(ensdb,filter=~ tx_id == use_enst &
                                   tx_biotype == "protein_coding"))
  ## merge exon
  total_exon_len <- sum(exon_tab$width) # range_min->range_max
  rr <- (range_max-range_min)/total_exon_len
  len_t <- function(x){
    w1 <- which(exon_tab$start-x<=0)
    w2 <- which(exon_tab$end-x>=0)
    w3 <- intersect(w1,w2)
    if(length(w1)==0) return(len_t(exon_tab$start[1]))
    if(length(w3)==0) return(len_t(exon_tab$end[max(w1)]))
    #print(x);print(w3)
    if(w3>1) cumsum(exon_tab$width)[(w3-1)]+(x-exon_tab$start[w3]) else (x-exon_tab$start[w3])
  }
  ss1 <- rr*unlist(lapply(exon_tab$start,len_t))+range_min
  ee1 <- rr*unlist(lapply(exon_tab$end,len_t))+range_min
  mRNA_start_pos=ss1;mRNA_end_pos=ee1
  if(draw_mRNA==TRUE) rect(xleft=ss1,xright=ee1,ybottom=-1.5,ytop=-1,col='light grey',border = 'blue')
  if(draw_ENST==TRUE) segments(x0=exon_tab$start/2+exon_tab$end/2,x1=ss1/2+ee1/2,y0=0,y1=-1,lwd=0.3,col=adjustcolor('blue',0.3),lty=1)
  if(draw_mRNA==TRUE) text(range_min-dd/4,-1.25,'mRNA',xpd=TRUE,adj=1,cex=0.7)
  if(strand=='+' & draw_mRNA==TRUE)text(ss1/2+ee1/2,-2,1:length(ss1),cex=0.5)
  if(strand=='-' & draw_mRNA==TRUE)text(ss1/2+ee1/2,-2,rev(1:length(ss1)),cex=0.5)
  ## mark cds
  ss1 <- rr*unlist(lapply(cds_tab$start,len_t))+range_min
  ee1 <- rr*unlist(lapply(cds_tab$end,len_t))+range_min
  if(draw_CDS==TRUE) rect(xleft=ss1,xright=ee1,ybottom=-3,ytop=-2.5,col='light grey',border = 'red')
  if(draw_CDS==TRUE) text(range_min-dd/4,-2.75,'CDS',xpd=TRUE,adj=1,cex=0.7)
  ##
  total_cds_len <- sum(cds_tab$width) # range_min->range_max
  ss <- seq(len_t(min(cds_tab$start)),len_t(max(cds_tab$end)),length.out=11);
  ss1 <- round(ss); ss1 <- ss1-min(ss1)+1
  if(draw_CDS==TRUE) segments(x0=ss*rr+range_min,x1=ss*rr+range_min,y0=-3.75,y1=-3.65,xpd=TRUE,cex=0.6)
  if(draw_CDS==TRUE) segments(x0=min(ss*rr+range_min),x1=max(ss*rr+range_min),y0=-3.75,y1=-3.75)
  if(strand == '+'){
    if(draw_CDS==TRUE) text(ss*rr+range_min,-4,get_label_manual(floor(ss1/3)-1),xpd=TRUE,cex=0.5)
    if(draw_CDS==TRUE) text(ss*rr+range_min,-3.5,get_label_manual(ss1),xpd=TRUE,cex=0.5)
  }
  if(strand == '-'){
    if(draw_CDS==TRUE) text(ss*rr+range_min,-4,rev(get_label_manual(floor(ss1/3)-1)),xpd=TRUE,cex=0.5)
    if(draw_CDS==TRUE) text(ss*rr+range_min,-3.5,rev(get_label_manual(ss1)),xpd=TRUE,cex=0.5)
  }
  ## mark domain
  if(is.null(use_domain)==FALSE){
    u1 <- unique(use_domain$Domain)
    cc1 <- colorRampPalette(brewer.pal(8,'Set3'))(length(u1));
    cc1 <- adjustcolor(cc1,alpha=0.8);names(cc1) <- u1;
    aa_start <- ss[1]+3; aa_end <- max(ss); aa_pos <- seq(aa_start,aa_end,length.out=max(ss1)/3-1)
    if(strand == '-')aa_pos <- rev(aa_pos)
    sss <- rep(1:5,length.out=nrow(use_domain))*0.1; sss <- 0
    if(draw_domain==TRUE) rect(xleft=aa_pos[as.numeric(use_domain$Start)]*rr+range_min,xright=aa_pos[as.numeric(use_domain$End)]*rr+range_min,
         ybottom=-5.5+sss,ytop = -5.5+sss+0.25,col=cc1[use_domain$Domain],border=cc1[use_domain$Domain])
    if(draw_domain==TRUE) legend(range_max,-4,fill=cc1,names(cc1),xpd=TRUE,border=NA,bty='n',cex=0.6,yjust = 1,
           title='Domain')
    #text(range_min-dd/4,-5.5,'Domain',xpd=TRUE,adj=1,cex=0.7)
    domain2color <- cc1
  }

  ## mark variants
  if(is.null(all_variant)==TRUE){
    return(list(domain2color=domain2color,mutationType2color=mutationType2color,
                                             mRNA_start_pos=mRNA_start_pos,mRNA_end_pos=mRNA_end_pos))
  }
  all_variant <- all_variant[rev(names(all_variant))]
  yy_s <- -6;  yy_e <- yy_s-step_h*length(all_variant);
  yy_u <- seq(yy_e,yy_s,by=step_h); names(yy_u) <- names(all_variant)
  # same variant
  cohort_name <- names(all_variant)[1]
  use_variant <- all_variant[[cohort_name]]
  if(count_log==TRUE) use_variant$count <- log(use_variant$count)/log(log_use)+1
  p1 <- unlist(lapply(use_variant$Pos,len_t))
  all_variant_mark <- lapply(all_variant,function(x)apply(x[,c('Chr','Pos','Ref','Alt')],1,function(x1)paste(x1,collapse = '-')))
  names(all_variant_mark) <- names(all_variant)
  variant_mark_count <- table(unlist(lapply(all_variant_mark,unique)));
  common_mark <- names(variant_mark_count)[which(variant_mark_count==length(all_variant))]
  use_variant <- all_variant[[cohort_name]]; w1 <- which(all_variant_mark[[cohort_name]] %in% common_mark);print(w1)
  if(mark_common==TRUE){
    segments(x0=(p1*rr+range_min)[w1],x1=(p1*rr+range_min)[w1],
             y0=rep(yy_s-2,length=length(w1)),y1=rep(yy_e,length=length(w1)),col=adjustcolor('grey',0.3),lwd=0.2)
  }
  # draw for each cohort
  u1 <- unique(unlist(lapply(all_variant,function(x)x$Mutation_type)))
  cc1 <- type2color[which(names(type2color) %in% u1)]
  mutationType2color <- cc1;
  ##
  mm_len <- max(unlist(lapply(all_variant,function(x)max(x$Pos))));
  mm_len <- max(mm_len,max(cds_tab$end))
  #-6-length(all_variant)*step_h
  legend(range_max,-6-length(all_variant)*step_h,
         fill=cc1,names(cc1),xpd=TRUE,border=NA,bty='n',cex=0.6,yjust = 1,
         title='Mutation Type')
  ##
  #p1_min <- min(unlist(lapply(all_variant,function(x)min(x$Pos))))
  for(cohort_name in names(all_variant)){
    use_variant <- all_variant[[cohort_name]]
    use_variant$count[which(is.na(use_variant$count)==T)] <- 0
    if(count_log==TRUE) use_variant$count <- log(use_variant$count)/log(log_use)+1
    p1 <- unlist(lapply(use_variant$Pos,len_t))
    max_count <- max(use_variant$count,na.rm=T)
    if(is.null(max_count_val)==TRUE | is.na(max_count_val)==T) max_count_val <- max_count
    segments(x0=p1*rr+range_min,x1=p1*rr+range_min,y0=yy_u[cohort_name],
             y1=yy_u[cohort_name]+0.01+0.9*step_h/2*use_variant$count/max_count_val,
             col=cc1[use_variant$Mutation_type],lwd=mark_lwd) ##
    text(range_min-dd/4,yy_u[cohort_name]+0.2,cohort_name,xpd=TRUE,adj=1,cex=0.7)
    segments(x0=len_t(mm_len)*rr+range_min,x1=len_t(mm_len)*rr+range_min,y0=yy_u[cohort_name],
             y1=yy_u[cohort_name]+0.9*step_h/2,lwd=0.5,col='grey')
    #segments(x0=range_min,x1=len_t(mm_len)*rr+range_min,y0=yy_u[cohort_name]+0.45,
    #         y1=yy_u[cohort_name]+0.45,lwd=0.3,lty=2,col='grey')
    segments(x0=range_min,x1=len_t(mm_len)*rr+range_min,y0=yy_u[cohort_name],
             y1=yy_u[cohort_name],lwd=0.5,col='grey')
    ss <- seq(yy_u[cohort_name],yy_u[cohort_name]+0.9*step_h/2,length.out = 3)
    segments(x0=len_t(mm_len)*rr+range_min,x1=len_t(mm_len)*rr+range_min+dd/50,y0=ss,y1=ss,lwd=0.5,col='grey')
    if(count_log==FALSE) text(len_t(mm_len)*rr+range_min+dd/40,ss,signif(seq(0,max_count_val,length.out=3),3),adj=0,xpd=TRUE,cex=0.35)
    if(count_log==TRUE)  text(len_t(mm_len)*rr+range_min+dd/40,ss,round(log_use^(seq(0,max_count_val,length.out=3)-1)),adj=0,xpd=TRUE,cex=0.35)
    if(mark_top_num>0){
      w1 <- rank(use_variant$count);w2 <- which(w1>nrow(use_variant)-mark_top_num)
      text(x=(p1*rr+range_min)[w2],
           y=yy_u[cohort_name]+0.05+0.95*step_h/2*use_variant$count[w2]/max_count_val,
           use_variant$Report[w2],adj=0,cex=mark_top_cex,srt=use_srt)
    }
    if('Pathogenity' %in% colnames(use_variant)){
      w1 <- which(use_variant$Pathogenity==1)
      text((p1*rr+range_min)[w1],yy_u[cohort_name]-0.05,'*',
           cex=0.3,xpd=TRUE,col=adjustcolor('red',0.8))
    }
  }
  return(list(domain2color=domain2color,mutationType2color=mutationType2color,
              mRNA_start_pos=mRNA_start_pos,mRNA_end_pos=mRNA_end_pos))
}
##
draw_gene_genome <- function(chr,pos_stt,pos_end,
                      only_use_trascript=FALSE,min_y=0,
                      tx_biotype_protein_coding=FALSE,
                      one_gene_one_tx=FALSE,transfer_tab=NULL){
  filter <- list(seqnames = c(chr))
  if(tx_biotype_protein_coding==TRUE){
    transcript_tab <- as.data.frame(transcripts(ensdb,filter= ~SeqNameFilter(chr) &
                                                  tx_biotype == "protein_coding"))
    exon_tab <- as.data.frame(exons(ensdb,filter= ~SeqNameFilter(chr) &
                                      tx_biotype == "protein_coding"))
    intron_tab <- as.data.frame(intronsByTranscript(ensdb,filter= ~SeqNameFilter(chr) &
                                                      tx_biotype == "protein_coding"))
  }else{
    transcript_tab <- as.data.frame(transcripts(ensdb,filter= ~SeqNameFilter(chr)))
    exon_tab <- as.data.frame(exons(ensdb,filter= ~SeqNameFilter(chr)))
    intron_tab <- as.data.frame(intronsByTranscript(ensdb,filter= ~SeqNameFilter(chr)))
  }
  transcript_tab <- transcript_tab[which(transcript_tab$start >= pos_stt
                                         & transcript_tab$end <= pos_end),]
  exon_tab <- exon_tab[which(exon_tab$start >= pos_stt
                                         & exon_tab$end <= pos_end),]
  intron_tab <- intron_tab[which(intron_tab$start >= pos_stt
                                         & intron_tab$end <= pos_end),]
  all_g <- unique(transcript_tab$gene_id)
  if(one_gene_one_tx==TRUE){
    tmp1 <- lapply(all_g,function(x){
      x1 <- transcript_tab[which(transcript_tab$gene_id==x),]
      x1$len <- x1$tx_cds_seq_end-x1$tx_cds_seq_start
      x1 <- x1[which.max(x1$len),]
    })
    transcript_tab <- do.call(rbind,tmp1)
    transcript_tab$tx_name <- transcript_tab$gene_id
    if(is.null(transfer_tab)==FALSE){
      transcript_tab$tx_name <- transfer_tab[transcript_tab$tx_name,]$external_gene_name
    }
  }
  ## draw exon structure for genes
  range_min <- pos_stt;
  range_max <- pos_end;
  dd <- 0;
  n_transcript <- nrow(transcript_tab)
  plot(1,xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i',bty='n',col='white',
       xlim=c(range_min-dd,range_max+dd*3),
       ylim=c(min_y,n_transcript/2+1))
  pp <- par()$usr ## merge exon
  # draw exon
  ss <- seq(range_min,range_max,length.out=11);
  rr <- floor(min(log10(ss)))-5
  ss1 <- round(ss/(10^rr))*10^rr
  axis(side=3,at=ss,labels=get_label_manual(ss1),xpd=TRUE,cex.axis=0.45)
  text(range_min-dd/4,(n_transcript+2)/2,sprintf('Chr%s',transcript_tab[1,1]),xpd=TRUE,cex=0.6,adj=1)
  rect(xleft=exon_tab$start,xright=exon_tab$end,ybottom=n_transcript/2,ytop=(n_transcript+1)/2,col='black',border = 'black')
  for(i in 1:nrow(intron_tab)){
    segments(x0=intron_tab[i,'start'],x1=intron_tab[i,'end'],y0=(n_transcript+0.5)/2,y1=(n_transcript+0.5)/2)
  }
  #
  all_enst <- transcript_tab$tx_id
  for(i in 1:length(all_enst)){
    each_enst <- all_enst[i]
    text(range_min-dd/4,(i-0.75)/2,transcript_tab$tx_name[i],adj=1,xpd=TRUE,cex=0.5)
    if(tx_biotype_protein_coding==TRUE){
      transcript_tab_1 <- as.data.frame(exons(ensdb,filter = ~ tx_id == each_enst &
                                                tx_biotype == "protein_coding"))
    }else{
      transcript_tab_1 <- as.data.frame(exons(ensdb,filter = ~ tx_id == each_enst))
    }
    intron_tab_1 <- intron_tab[which(intron_tab$group_name==each_enst),]
    if(nrow(intron_tab_1)>0){
      for(j in 1:nrow(intron_tab_1)){
        segments(x0=intron_tab_1[j,'start'],x1=intron_tab_1[j,'end']/2+intron_tab_1[j,'start']/2,y0=(i-1)/2,y1=(i-0.5)/2,col='dark grey')
        segments(x1=intron_tab_1[j,'end'],x0=intron_tab_1[j,'end']/2+intron_tab_1[j,'start']/2,y0=(i-0.5)/2,y1=(i-1)/2,col='dark grey')
      }
    }
    for(j in 1:nrow(transcript_tab_1)){
      rect(xleft=transcript_tab_1[j,'start'],
           xright=transcript_tab_1[j,'end'],
           ybottom=(i-0.5)/2,ytop=(i-1)/2,col='black',border = 'black')
    }
  }
  ##
}
##############################################
draw.heatmap.local <- function(mat,inner_line=FALSE,out_line=TRUE,inner_col='black',
                               n=20,col_srt=60,display_text_mat=NULL,row_cex=1,column_cex=1,text_cex=1,text_col='up',
                               bb_max=NULL,main='main'){
  if(ncol(mat)>1) mat1 <- mat[nrow(mat):1,] else mat1 <- as.matrix(mat[nrow(mat):1,])
  if(is.null(display_text_mat)==FALSE){
    if(ncol(display_text_mat)>1) display_text_mat <- display_text_mat[nrow(display_text_mat):1,] else display_text_mat <- as.matrix(display_text_mat[nrow(display_text_mat):1,])
  }
  colnames(mat1) <- colnames(mat)
  #cc1 <- grDevices::colorRampPalette(c(brewer.pal(9,'Blues')[8],'white',brewer.pal(9,'Reds')[7]))(n*2)
  cc1 <- c('white',colorRampPalette(rev(brewer.pal(11,'Spectral')[1:7]))(n-1))
  if(is.null(bb_max)==TRUE) bb_max <- base::max(abs(mat1),na.rm=TRUE)
  bb1 <- c(-1,base::seq(0,bb_max,length.out=n))
  #mat1[which(is.na(mat1)==TRUE)] <- 0;
  if(out_line==TRUE) graphics::image(t(mat1),col=cc1,breaks=bb1,xaxt='n',yaxt='n',main=main)
  if(out_line==FALSE) graphics::image(t(mat1),col=cc1,breaks=bb1,xaxt='n',yaxt='n',bty='n',main=main)
  pp <- par()$usr
  xx <- base::seq(pp[1],pp[2],length.out=1+ncol(mat1))
  yy <- base::seq(pp[3],pp[4],length.out=1+nrow(mat1))
  xxx <- xx[1:(base::length(xx)-1)]/2+xx[2:base::length(xx)]/2
  yyy <- yy[1:(base::length(yy)-1)]/2+yy[2:base::length(yy)]/2
  if(inner_line==TRUE){
    graphics::abline(v=xx,col=inner_col)
    graphics::abline(h=yy,col=inner_col)
  }
  graphics::text(pp[1]-par.char2pos()[1],yyy,rownames(mat1),adj=1,xpd=TRUE,cex=row_cex) ## distance for one character
  if(text_col=='up'){
    if(col_srt==0){
      graphics::text(xxx,pp[4],colnames(mat1),pos=3,xpd=TRUE,srt=col_srt,cex=column_cex) ## do not use pos, for srt?
    }else{
      graphics::text(xxx,pp[4]+base::max(strheightMod(colnames(mat1))/2)+base::max(strheightMod(colnames(mat1))/10),colnames(mat1),adj=0.5-col_srt/90,xpd=TRUE,srt=col_srt,cex=column_cex) ## do not use pos, for srt?
    }
  }
  if(text_col=='down'){
    if(col_srt!=0) graphics::text(xxx,pp[3]-0.1*(pp[4]-pp[3])/nrow(mat1),colnames(mat1),adj=1,xpd=TRUE,srt=col_srt,cex=column_cex)
    if(col_srt==0) graphics::text(xxx,pp[3]-0.1*(pp[4]-pp[3])/nrow(mat1),colnames(mat1),adj=0.5,xpd=TRUE,srt=col_srt,cex=column_cex)
  }
  if(is.null(display_text_mat)==FALSE){
    for(i in 1:nrow(display_text_mat)){
      for(j in 1:ncol(display_text_mat)){
        graphics::text(xxx[j],yyy[i],display_text_mat[i,j],cex=text_cex)
      }
    }
  }
  return(TRUE)
}
draw.colorbar <- function(col=NULL,min_val=NULL,max_val=NULL,n=5,digit_num=2,direction='vertical',xleft=0,xright=1,ytop=1,ybottom=0,cex=1){
  val_bar<-base::seq(max_val,min_val,length.out = n)
  if(is.null(col)==TRUE) col <- z2col(val_bar)
  if(direction=='vertical'){
    y_pos <- base::seq(ytop,ybottom,length.out=n+1)
    yy <- y_pos[2:base::length(y_pos)]/2+y_pos[1:(base::length(y_pos)-1)]/2
    graphics::rect(xleft=xleft,xright=xright,ytop=y_pos[2:base::length(y_pos)],ybottom=y_pos[1:(base::length(y_pos)-1)],col=col,border='light grey',xpd=TRUE)
    graphics::text(xright,yy,signif(val_bar,digits = digit_num),xpd=TRUE,cex=cex,pos=4)
    graphics::text(xleft/2+xright/2,ytop,'Z value',cex=cex,xpd=TRUE,pos=3)
  }
}

par.lineHeight2inch <- function(){
  lheight <- par()$lheight
  y1 <- par.char2inch()[2]*lheight ## line height in inches
  y1
}
par.char2pos <- function(){par()$cxy}
strheightMod <- function(s, units = "inch", cex = 1,ori=TRUE,mod=FALSE){
  if(ori==TRUE) return(strheight(s=s,units=units,cex=cex))
  if(units=='user') return(par.char2pos()[2]*cex)
  if(units=='inch' | units=='inches') return(par.char2inch()[2]*cex)
}
strwidthMod <- function(s, units = "inch", cex = 1,ori=TRUE,mod=FALSE){
  if(ori==TRUE) return(strwidth(s=s,units=units,cex=cex))
  if(mod==TRUE){
    plot.new()
    rt <- strwidth(s,units=units)/strwidth('W',units=units); rt <- ceiling(rt)
    if(units=='user') r1 <- par.char2pos()[1]*cex*rt
    if(units=='inch') r1 <- par.char2inch()[1]*cex*rt
    dev.off(); return(r1)
  }else{
    if(units=='user') return(par.char2pos()[1]*cex*nchar(s))
    if(units=='inch'| units=='inches') return(par.char2inch()[1]*cex*nchar(s))
  }
}
##
draw_haploptype <- function(use_gene,use_ensg,use_enst,use_domain,all_hap,all_maf,all_hap_mark,mark_site){
  transcript_tab <- as.data.frame(transcripts(ensdb,filter = ~ gene_id == use_ensg &
                                                tx_biotype == "protein_coding"))
  exon_tab <- as.data.frame(exons(ensdb,filter = ~ gene_id == use_ensg &
                                    tx_biotype == "protein_coding"))
  intron_tab <- as.data.frame(intronsByTranscript(ensdb,filter = ~ gene_id == use_ensg &
                                                    tx_biotype == "protein_coding"))
  cds_tab <- as.data.frame(cdsBy(ensdb,filter=~ gene_id == use_ensg &
                                   tx_biotype == "protein_coding"))
  ## draw exon structure for a gene
  range_min <- min(exon_tab$start);
  range_max <- max(exon_tab$end);
  dd <- (range_max-range_min)/10;
  n_transcript <- nrow(transcript_tab)
  plot(1,xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i',bty='n',col='white',xlim=c(range_min-dd,range_max+dd*0.5),
       ylim=c(-2-length(all_hap)*2,n_transcript/2+1))
  pp <- par()$usr ## merge exon
  # draw exon
  ss <- seq(range_min,range_max,length.out=11);
  rr <- floor(min(log10(ss)))-5
  ss1 <- round(ss/(10^rr))*10^rr
  axis(side=3,at=ss,labels=get_label_manual(ss1),xpd=TRUE,cex=0.8)
  text(range_min-dd/4,(n_transcript+2)/2,sprintf('Chr%s',transcript_tab[1,1]),xpd=TRUE,cex=0.8,adj=1)
  text(range_min-dd/4,(n_transcript+1)/2,use_ensg,xpd=TRUE,cex=0.8,adj=1)
  rect(xleft=exon_tab$start,xright=exon_tab$end,ybottom=n_transcript/2,ytop=(n_transcript+1)/2,col='black',border = 'black')
  for(i in 1:nrow(intron_tab)){
    segments(x0=intron_tab[i,'start'],x1=intron_tab[i,'end'],y0=(n_transcript+0.5)/2,y1=(n_transcript+0.5)/2)
  }
  #
  all_enst <- unique(transcript_tab$tx_id)
  all_enst <- c(use_enst,setdiff(all_enst,use_enst))
  for(i in 1:length(all_enst)){
    each_enst <- all_enst[i]
    if(each_enst==use_enst){
      text(range_min-dd/4,(i-0.75)/2,each_enst,adj=1,xpd=TRUE,cex=0.7,col='blue')
    }else{
      text(range_min-dd/4,(i-0.75)/2,each_enst,adj=1,xpd=TRUE,cex=0.6)
    }
    transcript_tab_1 <- as.data.frame(exons(ensdb,filter = ~ tx_id == each_enst &
                                              tx_biotype == "protein_coding"))
    intron_tab_1 <- intron_tab[which(intron_tab$group_name==each_enst),]
    for(j in 1:nrow(transcript_tab_1)){
      rect(xleft=transcript_tab_1[j,'start'],xright=transcript_tab_1[j,'end'],
           ybottom=(i-0.5)/2,ytop=(i-1)/2,col='black',border = 'black')
    }
    for(j in 1:nrow(intron_tab_1)){
      segments(x0=intron_tab_1[j,'start'],x1=intron_tab_1[j,'end']/2+intron_tab_1[j,'start']/2,y0=(i-1)/2,y1=(i-0.5)/2,col='dark grey')
      segments(x1=intron_tab_1[j,'end'],x0=intron_tab_1[j,'end']/2+intron_tab_1[j,'start']/2,y0=(i-0.5)/2,y1=(i-1)/2,col='dark grey')
    }
  }
  ## pfam
  # merge exon, need to specify transcript
  exon_tab <- as.data.frame(exons(ensdb,filter = ~ tx_id == use_enst &
                                    tx_biotype == "protein_coding"))
  intron_tab <- as.data.frame(intronsByTranscript(ensdb,filter = ~ tx_id == use_enst &
                                                    tx_biotype == "protein_coding"))
  cds_tab <- as.data.frame(cdsBy(ensdb,filter=~ tx_id == use_enst &
                                   tx_biotype == "protein_coding"))
  ## merge exon
  total_exon_len <- sum(exon_tab$width) # range_min->range_max
  rr <- (range_max-range_min)/total_exon_len
  len_t_p2d <- function(x){
    #print(x)
    x <- as.numeric(x)
    x1 <- cumsum(cds_tab$width)/3
    w1 <- max(which(x1<x))+1
    cds_tab$start[w1]+(x-x1[w1-1])*3
  }
  ## mark domain
  u1 <- unique(use_domain$Domain)
  cc1 <- colorRampPalette(brewer.pal(8,'Set3'))(length(u1));cc1 <- adjustcolor(cc1,alpha=0.8);names(cc1) <- u1;
  aa_start <- ss[1]+3; aa_end <- max(ss); aa_pos <- seq(aa_start,aa_end,length.out=max(ss1)/3-1)
  sss <- rep(1:5,length.out=nrow(use_domain))*0.1; sss <- 0;
  rect(xleft=unlist(lapply(pf1$Start,len_t_p2d)),xright=lapply(pf1$End,len_t_p2d),
       ybottom=-1.5+sss,ytop = -1.5+sss+0.3,col=cc1[use_domain$Domain],border=cc1[use_domain$Domain])
  legend(range_max,-1.5+sss,fill=cc1,names(cc1),xpd=TRUE,border=NA,bty='n',cex=1.2,
         yjust = 0,xjust=1,
         title='Domain')
  #text(range_min-dd/4,-1.5,'Domain',xpd=TRUE,adj=1,cex=1.4)
  ## haplotypes
  cc1 <- brewer.pal(8,'Set1')[5:4]
  cc2 <- c('green','blue','yellow','red');names(cc2) <- c('A','T','C','G')
  yy_s <- -2;  yy_e <- yy_s-2*length(all_hap); yy_u <- seq(yy_e,yy_s,by=2); names(yy_u) <- names(all_hap)
  for(i in 1:length(all_hap)){
    use_cohort <- names(all_hap)[i]
    h1 <- all_hap[[i]]
    rect(xleft=h1$BP1,xright=h1$BP2,ybottom=yy_u[i]+0.1,ytop = yy_u[i]+0.2,border=NA,lwd=0.3,col=adjustcolor('black',0.6),xpd=TRUE)
    text(range_min-dd/4,yy_u[i]+0.2,names(all_hap)[i],xpd=TRUE,adj=1,cex=1.2)
    mf1 <- all_maf[[use_cohort]]
    mm1 <- all_hap_mark[[use_cohort]]
    cc_r <- unlist(lapply(mf1$POS,function(x){
      w1 <- which(exon_tab$start>x)
      w2 <- which(exon_tab$end<x)
      w3 <- setdiff(1:nrow(exon_tab),c(w1,w2))
      if(length(w3)>0) cc1[1] else cc1[2]
    }))
    segments(x0=mf1$POS,x1=mf1$POS,y0=yy_u[i]+0.2,y1=yy_u[i]+mf1$ALT_AF+0.2,lwd=ifelse(cc_r == cc1[1],0.9,0.6),
             col=ifelse(cc_r == cc1[1],adjustcolor(cc_r,0.8),adjustcolor(cc_r,0.3)))
    w1 <- which(h1$KB>1)
    text(x=(h1$BP1/2+h1$BP2/2)[w1],y=yy_u[i]-0.05,sprintf('Block%s',w1),cex=0.8,xpd=T)
    #points(mf1$POS,yy_u[i]+mf1$ALT_AF+0.2,pch=16,cex=0.05,col=1)
    #segments(x0=mm1$POS,x1=mm1$POS,y0=yy_u[i]-0.05,y1=yy_u[i]+0.15,lwd=0.2,col=adjustcolor(cc2[mm1$A2],0.4),xpd=T)
    #points(mm1$POS,rep(yy_u[i]+0.05,length.out=length(mm1$POS)),pch=16,col=adjustcolor(cc2[mm1$A2],0.4),xpd=T)
    #for(j in 1:nrow(mm1)) draw.ellipse(x=mm1$POS[j],y=yy_u[i]-0.05,a=250,b=0.08,col=adjustcolor(cc2[mm1$A2[j]],0.4),xpd=T,border=NA)
    for(j in 1:length(mark_site)){
      w1 <- which(mf1$POS %in% mark_site[[j]])
      if(names(mark_site)[j]==names(all_hap)[i]) points(mf1$POS[w1],yy_u[i]+mf1$ALT_AF[w1]+0.4,pch='*',cex=0.8,col=2)
    }
  }
  legend(range_max,-7,fill=cc1[1:2],c('exon SNP','intron SNP'),xpd=TRUE,border=NA,bty='n',cex=1.3,yjust = 1,xjust=1)
  #legend(range_max,-9,fill=cc2,names(cc2),xpd=TRUE,border=NA,bty='n',cex=0.8,yjust = 1)
  legend(range_max,-12,pch='*','Diff SNP \n(CCGT Vs. EUR)',col=2,xpd=TRUE,border=NA,bty='n',cex=1.3,yjust = 1,xjust=1)
  ## ACTG green blue yellow red
}
##
par.pos2inch <- function(){
  user.range <- par("usr")[c(2,4)] - par("usr")[c(1,3)]
  region.pin <- par("pin")
  return(region.pin/user.range)
}

add_m <- function(x,log='',cex=2,xadjust=TRUE){
  pp <- par()$usr;ppm <- par()$mai
  ddx <- strwidth(x,cex=cex)*2
  ddy <- strheight(x,cex=cex)
  if(xadjust==TRUE){
    if(log=='') text(pp[1]-ppm[2]/par.pos2inch()[1]+ddx,pp[4]+ppm[3]/par.pos2inch()[2]-ddy,x,xpd=TRUE,font=2,cex=cex,pos=2)
    if(log=='x')  text(10^(pp[1]-ppm[2]/par.pos2inch()[1]+ddx),pp[4]+ppm[3]/par.pos2inch()[2]-ddy,x,xpd=TRUE,font=2,cex=cex,pos=2)
    if(log=='y')  text(pp[1]-ppm[2]/par.pos2inch()[1]+ddx,10^(pp[4]+ppm[3]/par.pos2inch()[2]-ddy),x,xpd=TRUE,font=2,cex=cex,pos=2)
  }else{
    if(log=='') text(pp[1],pp[4]+ppm[3]/par.pos2inch()[2]-ddy,x,xpd=TRUE,font=2,cex=cex,pos=2)
    if(log=='x')  text(10^(pp[1]),pp[4]+ppm[3]/par.pos2inch()[2]-ddy,x,xpd=TRUE,font=2,cex=cex,pos=2)
    if(log=='y')  text(pp[1],10^(pp[4]+ppm[3]/par.pos2inch()[2]-ddy),x,xpd=TRUE,font=2,cex=cex,pos=2)
  }
}
##
draw_AF_tab <- function(r1,min_val=7,text_cex=0.5){
  bb <- 10^(-min_val:0)
  cc <- colorRampPalette(brewer.pal(8,'Reds'))(length(bb)-1)
  par(mar=c(8,10,3,3))
  image(t(r1),xaxt='n',yaxt='n',bty='n',col=cc,breaks=bb)
  pp <- par()$usr;
  xx1 <- seq(pp[1],pp[2],length.out=ncol(r1)+1);xx<-xx1[1:(length(xx1)-1)]/2+xx1[2:length(xx1)]/2
  yy1 <- seq(pp[3],pp[4],length.out=nrow(r1)+1);yy<-yy1[1:(length(yy1)-1)]/2+yy1[2:length(yy1)]/2
  text(x=pp[1],y=yy,rownames(r1),adj=1,xpd=T)
  text(y=pp[3],x=xx,colnames(r1),adj=1,xpd=T,srt=60)
  abline(h=yy1);abline(v=xx1)
  for(i in 1:nrow(r1)){
    for(j in 1:ncol(r1)){
      if(r1[i,j]>0)text(y=yy[i],x=xx[j],signif(r1[i,j],2),cex=text_cex)
    }
  }
}
draw_incidence <- function(inc_each,choose='AF'){
  par(mar=c(8,10,3,3))
  cc <- brewer.pal(8,'Paired')[1];
  r1<-inc_each[,1];r2 <- inc_each[,2];r3 <- inc_each[,3]
  if(choose=='AF'){
    a <- barplot(r1,las=2,beside = T,axisnames = F,col=cc,border = cc,
                 ylim=c(0,1.25*max(inc_each)),cex.axis=0.8,yaxt='n')
    axis(side=2,at=c(0,1e-7,1e-6,1e-5,2e-5,5e-5),
         labels=format(c(0,1e-7,1e-6,1e-5,2e-5,5e-5),scientific=T),
         las=2,cex.axis=0.8)
    mtext(side=2,line=4,'Estimated affected frequency by\nBayesian framework')
  }
  if(choose=='count'){
    a <- barplot(r1,las=2,beside = T,axisnames = F,col=cc,border = cc,
                 ylim=c(0,1.25*max(r1)),cex.axis=0.8)
    mtext(side=2,line=4,'Estimated Patient number to\nobserve one affected case by\nBayesian framework')
  }
  pp <- par()$usr
  text(a,0,rownames(inc_each),xpd=T,adj=1,cex=0.7,srt=60)
  segments(x0=a,x1=a,y0=r2,y1=r3,xpd=T,col='dark grey')
  segments(x0=a,x1=a,y0=r2,y1=r3,xpd=T,col='dark grey')
}
###
Incidence.get_risk <- function(n,x){
  c(x/n,(x/n)^2,(x/n)^2/4)
}
Incidence.get_inc <- function(m,x){
  q2P_025=qbeta(0.025,x,m-x)^2
  q2P_975=qbeta(0.975,x,m-x)^2
  alphaP=x;betaP=m-x
  q2P=(alphaP/(alphaP+betaP))^2+alphaP*betaP/((alphaP+betaP)^2*(alphaP+betaP+1))
  c(q2P,q2P_025,q2P_975)
}
Incidence.get_p <- function(m,x1,x2){ # x1>x2, m: AN
  r1 <- get_inc(m,x1)
  r2 <- get_inc(m,x2)
  p1 <- pbeta(sqrt(r2[1]),x1,m-x1)
  p2 <- 1-pbeta(sqrt(r1[1]),x2,m-x2)
  print(c(p1,p2))
  if(p1>0.5){p1<- 1-p1}
  if(p2>0.5){p2<- 1-p2}
  mean(c(p1,p2))
  #combinePvalVector(c(p1,p2))[2]
}

