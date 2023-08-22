library(XML)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
ensdb <- EnsDb.Hsapiens.v75
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
###
source('D:/analysis_eng/SherryPlot/function_basic.R')
get_domain <- function(use_uniprot, only_pfam=FALSE){
  ori_pf1 <- readHTMLTable(sprintf('http://pfam.xfam.org/protein/%s',use_uniprot),
                           stringsAsFactors=F)
  pf1 <- ori_pf1$imageKey[,1:4];colnames(pf1) <- c('Source','Domain','Start','End')
  pf1$Domain[which(pf1$Domain=='n/a')] <- pf1$Source[which(pf1$Domain=='n/a')]
  pf1$Start <- as.numeric(pf1$Start)
  pf1$End <- as.numeric(pf1$End)
  if(only_pfam==TRUE){
    pf1 <- pf1[which(pf1$Source == 'Pfam'),]
  }
  pf1
}
get_exon <- function(use_enst){
  exon_tab <- as.data.frame(exons(ensdb,filter = ~ tx_id == use_enst &
                                    tx_biotype == "protein_coding"))
  exon_tab
}
get_cds <- function(use_enst){
  cds_tab <- as.data.frame(cdsBy(ensdb,filter=~ tx_id == use_enst &
                                   tx_biotype == "protein_coding"))
  cds_tab
}
get_transcript <- function(use_enst){
  transcript_tab <- as.data.frame(transcripts(ensdb,filter = ~ tx_id == use_enst &
                                                tx_biotype == "protein_coding"))
  transcript_tab
}
Pos2RNA_exon_start <- function(cds_start,cds_end,first_pos,inner_size=0,strand='+'){
  cds_start <- sort(cds_start)
  cds_end <- sort(cds_end)
  if(strand=='+'){
    cds_width <- cds_end-cds_start+1
    cds_cum   <- cumsum(cds_width)
    cds_cum   <- c(1,cds_cum[1:(length(cds_cum)-1)])
    r1 <- cds_cum+first_pos
    r2 <- r1+inner_size*(0:(length(r1)-1))
    r2
  }else{
    cds_width <- rev(cds_end-cds_start+1)
    cds_cum   <- cumsum(cds_width)
    cds_cum   <- c(1,cds_cum[1:(length(cds_cum)-1)])
    r1 <- cds_cum+first_pos
    r2 <- r1+inner_size*(0:(length(r1)-1))
    r2
  }
}
Pos2RNA_exon_end <- function(cds_start,cds_end,first_pos,inner_size=0,strand='+'){
  cds_start <- sort(cds_start)
  cds_end <- sort(cds_end)
  if(strand=='+'){
    cds_width <- cds_end-cds_start+1
    cds_cum   <- cumsum(cds_width)
    r1 <- cds_cum+first_pos
    r2 <- r1+inner_size*(0:(length(r1)-1))
    r2
  }else{
    cds_width <- rev(cds_end-cds_start+1)
    cds_cum   <- cumsum(cds_width)
    r1 <- cds_cum+first_pos
    r2 <- r1+inner_size*(0:(length(r1)-1))
    r2
  }
}
## pos: position from cds start --> position from rna start
Pos2RNA <- function(pos,cds_start,cds_end,first_pos,inner_size=0,strand='+'){
  cds_start <- sort(cds_start)
  cds_end <- sort(cds_end)
  res <- c()
  for(i in pos){
    cds_width <- cds_end-cds_start+1
    if(strand == '-') cds_width <- rev(cds_width)
    cds_cum   <- cumsum(cds_width)
    w1 <- min(which(cds_cum>=i))
    w2 <- i-cds_cum[w1]
    r1 <- cds_cum[w1]+w2+first_pos
    r2 <- r1+inner_size*(w1-1)
    res <- c(res,r2)
  }
  res
}
## pos: position from mRNA start --> output genome position
Pos2Genome <- function(pos,cds_start,cds_end,strand='+'){
  cds_start <- sort(cds_start)
  cds_end <- sort(cds_end)
  res <- c()
  cds_width <- cds_end-cds_start+1
  if(strand == '-'){cds_width <- rev(cds_width);
  cds_start <- sort(cds_start,decreasing = T);
  cds_end <- sort(cds_end,decreasing = T)}
  cds_cum   <- cumsum(cds_width)
  if(strand == '+'){
    for(i in pos){
      w1 <- min(which(cds_cum>=i))
      if(w1!=1){
        w2 <- i-cds_cum[w1-1]
      }else{
        w2 <- i
      }
      res <- c(res,cds_start[w1]+w2)
    }
  }
  if(strand == '-'){
    for(i in pos){
      w1 <- min(which(cds_cum>=i))
      if(w1!=1){
        w2 <- i-cds_cum[w1-1]
      }else{
        w2 <- i
      }
      res <- c(res,cds_end[w1]-w2)
    }
  }
  res
}
query_pos <- function(x,pos_start,pos_end,ID){
  w1 <- which(x>=pos_start & x<=pos_end)
  if(length(w1)>0) return(ID[w1]) else return('other')
}
## pos: genome position --> position from RNA start
Pos.Genome2RNA <- function(pos,exon_start,exon_end,strand='+'){
  exon_start <- sort(exon_start)
  exon_end <- sort(exon_end)
  exon_width <- exon_end-exon_start+1
  if(strand == '-') exon_width <- rev(exon_width)
  exon_cum   <- cumsum(exon_width)
  res <- list()
  for(ii in 1:length(pos)){
    i <- pos[ii]
    print(i)
    if(strand=='+'){
      w1 <- query_pos(ii,exon_start,exon_end,1:length(exon_end))
      if(w1=='other'){
        w21 <- abs(i-exon_start)
        w22 <- abs(i-exon_end)
        if(min(w21)<min(w22)){
          if(which.min(w21)==1) res[[ii]] <- 1
          if(which.min(w21)>1) res[[ii]] <- exon_cum[which.min(w21)-1]
        }else{
          res[[ii]] <- exon_cum[which.min(w22)]
        }
      }else{
        if(w1>1) res[[ii]] <- exon_cum[w1-1]+i-exon_start[w1]+1
        if(w1==1) res[[ii]] <- i-exon_start[w1]+1
      }
    }else{
      w1 <- query_pos(i,exon_start,exon_end,length(exon_end):1)
      if(w1=='other'){
        w21 <- abs(i-rev(exon_start))
        w22 <- abs(i-rev(exon_end))
        if(min(w21)<min(w22)){
          res[[ii]] <- exon_cum[which.min(w21)]
        }else{
          if(which.min(w22)==1) res[[ii]] <- 1
          if(which.min(w22)>1) res[[ii]] <- exon_cum[which.min(w22)-1]
        }
      }else{
        if(w1==1) res[[ii]] <- rev(exon_end)[w1]-i+1
        if(w1>1) res[[ii]] <- exon_cum[w1-1]+rev(exon_end)[w1]-i+1
      }
    }
  }
  unlist(res)
}
get_use_cdstab <- function(transcript_tab,cds_tab){
  strand <- transcript_tab[1,'strand']
  tss <- transcript_tab[1,]$tx_cds_seq_start
  ted <- transcript_tab[1,]$tx_cds_seq_end
  use_cds_tab <- cds_tab[which(cds_tab$end >= tss & cds_tab$start <= ted),]
  n <- nrow(use_cds_tab)
  if(use_cds_tab$start[1]<tss){use_cds_tab$start[1] <- tss}
  if(use_cds_tab$end[n]>ted){use_cds_tab$end[n] <- ted} #
  use_cds_tab
}
annotate_domain <- function(transcript_tab,exon_tab,cds_tab,pf,inner_size=0){
  strand <- as.character(cds_tab[1,'strand'])
  use_cds_tab <- get_use_cdstab(transcript_tab,cds_tab)
  if(strand == '+'){
    first_pos <- min(use_cds_tab$start)-min(exon_tab$start)+1
  }
  if(strand == '-'){
    first_pos <- max(exon_tab$end)-max(use_cds_tab$end)+1
  }
  print(first_pos)
  pf$RNA_Start_pos <- Pos2RNA(pf$Start*3+1,use_cds_tab$start,use_cds_tab$end,
                              first_pos,inner_size=inner_size,strand=strand)
  pf$RNA_End_pos <- Pos2RNA(pf$End*3+1,use_cds_tab$start,use_cds_tab$end,
                            first_pos,inner_size=inner_size,strand=strand)
  pf$Start_pos <- Pos2Genome(pf$RNA_Start_pos,use_cds_tab$start,use_cds_tab$end,strand=strand)
  pf$End_pos <- Pos2Genome(pf$RNA_End_pos,use_cds_tab$start,use_cds_tab$end,strand=strand)
  return(pf)
}
mod_variant_info <- function(variant_info,
                             Sample_col='',
                             Var_col='',
                             CDS_pos_col='',
                             Outcome_col='',
                             MutationType_col='',
                             Zygosity_col=''){
  r1 <- data.frame(Sample=variant_info[,Sample_col],
                   Variant=variant_info[,Var_col],
                   CDS_pos=as.numeric(variant_info[,CDS_pos_col]),
                   Outcome=variant_info[,Outcome_col],
                   MutationType=variant_info[,MutationType_col],
                   Zygosity=variant_info[,Zygosity_col],
                   stringsAsFactors = F)
  r1
}
draw_domain_BySample <- function(transcript_tab,exon_tab,cds_tab,pf,
                        variant_info,mark_outcome=NULL,
                        unit=0.35,unit_s=0.4,inner_size=0){
  strand <- as.character(cds_tab[1,'strand'])
  use_cds_tab <- get_use_cdstab(transcript_tab,cds_tab)
  if(strand == '+'){
    first_pos <- min(use_cds_tab$start)-min(exon_tab$start)+1
  }
  if(strand == '-'){
    first_pos <- max(exon_tab$end)-max(use_cds_tab$end)+1
  }
  CRNA_Start_pos <- Pos2RNA_exon_start(use_cds_tab$start,use_cds_tab$end,first_pos,inner_size=inner_size,strand=strand)
  CRNA_End_pos <- Pos2RNA_exon_end(use_cds_tab$start,use_cds_tab$end,first_pos,inner_size=inner_size,strand=strand)
  ERNA_Start_pos <- Pos2RNA_exon_start(exon_tab$start,exon_tab$end,first_pos=0,inner_size=inner_size,strand=strand)
  ERNA_End_pos <- Pos2RNA_exon_end(exon_tab$start,exon_tab$end,first_pos=0,inner_size=inner_size,strand=strand)
  pf <- annotate_domain(transcript_tab=transcript_tab,exon_tab=exon_tab,
                        cds_tab=cds_tab,pf=pf,inner_size=inner_size)
  ## draw
  all_s <- unique(variant_info$Sample)
  n <- nrow(exon_tab)
  ud <- unique(pf$Domain)
  stt = ERNA_Start_pos[1]
  end = ERNA_End_pos[n]
  mm <- min(pf$RNA_Start_pos)
#  unit <- 0.35
#  unit_s <- 0.4
  par(mar=c(1,2,1,4))
  cc <- brewer.pal(9,'Pastel1')
  plot_new(xlim=c(stt,end),ylim = c(-1.1-length(ud)*unit,
                                    length(all_s)*unit_s+unit_s))
  rect(xleft=ERNA_Start_pos,
       xright=ERNA_End_pos,
       ybottom=-0.5,ytop=0.1,col=cc[1],border = 'black')
  if(inner_size>0){
    ww <- 2:length(ERNA_End_pos)
    segments(x0=ERNA_End_pos[ww-1],x1=ERNA_Start_pos[ww],
             y0=-0.2,y1=-0.2)
  }
  rect(xleft=CRNA_Start_pos[1],
       xright=max(CRNA_End_pos),
       ybottom=-0.8,ytop=-0.6,col=cc[2],border = cc[2])
  text(ERNA_Start_pos/2+ERNA_End_pos/2,y=-0.2,
       sprintf('E%s',1:length(ERNA_Start_pos)),cex=0.5)
  text(CRNA_Start_pos[1]/2+max(CRNA_End_pos)/2,y=-0.7,
       'CDS',cex=0.5)
  cc <- get.class.color(ud)
  for(i in 1:length(ud)){
    rect(xleft=pf$RNA_Start_pos[which(pf$Domain == ud[i])],
         xright=pf$RNA_End_pos[which(pf$Domain == ud[i])],
         ybottom=-0.9-unit*i,ytop=-0.9-unit*i+unit,col=cc[i],border = cc[i])
    text(x=mm,y=-0.9-unit*i+unit/2,pos=2,xpd=TRUE,ud[i],cex=0.5)
  }
  ###
  all_mut <- unique(variant_info$MutationType)
  cc1 <- get.class.color(all_mut,use_color = brewer.pal(8,'Set2'))
  yy <- seq(from=unit_s,by=unit_s,length.out=length(all_s))
  mm <- min(variant_info$CDS_pos)
  for(i in 1:length(all_s)){
    x1 <- variant_info[which(variant_info$Sample == all_s[i]),]
    if(nrow(x1)!=1){
      segments(x0=min(x1$CDS_pos+first_pos),x1=max(x1$CDS_pos+first_pos),
               y0=yy[i],y1=yy[i],lwd=0.5,lty=2,col='grey')
    }
    points(x=x1$CDS_pos+first_pos,
           y=rep(yy[i],nrow(x1)),
           pch=ifelse(x1$Zygosity == 'Heterozygous',17,16),
           xpd=TRUE,col=cc1[x1$MutationType],
           cex=ifelse(x1$Zygosity == 'Heterozygous',0.8,1.1),)
    if(is.null(mark_outcome)==FALSE){
      if(mark_outcome %in% x1$Outcome){
        text(mm+first_pos,yy[i],sprintf('%s*',all_s[i]),pos=2,cex=0.5,
             xpd=T,col='red')
      }else{
        text(mm+first_pos,yy[i],all_s[i],pos=2,cex=0.5,xpd=T)
      }
    }else{
      text(mm+first_pos,yy[i],all_s[i],pos=2,cex=0.5,xpd=T)
    }
  }
  ## legend
  mm <- max(CRNA_End_pos)
  legend(x=mm,y=length(all_s)/2*unit_s,fill=cc1,names(cc1),
         xjust=0,yjust=1,cex=0.5,xpd=T,border = NA,bty='n')
  legend(x=mm,y=length(all_s)/2*unit_s,pch=c(17,16),
         c('Heterozygous','Homozygous'),
         xjust=0,yjust=0,cex=0.5,xpd=T,border = NA,bty='n')
  mm <- min(variant_info$CDS_pos)
  if(is.null(mark_outcome)==FALSE){
    text(x=mm,y=length(all_s)/2*unit_s,
         sprintf('*%s',mark_outcome),col=2,adj=1,
         xpd=T,cex=0.5)
  }
}
## mark_type, sample/variant
draw_domain_BySite <- function(transcript_tab,exon_tab,cds_tab,pf,
                                 variant_info,mark_outcome=NULL,mark_type='variant',
                                 unit=0.35,unit_s=0.4,inner_size=0,use_color=TRUE){
  variant_info <- variant_info[order(variant_info$CDS_pos),]
  strand <- as.character(cds_tab[1,'strand'])
  use_cds_tab <- get_use_cdstab(transcript_tab,cds_tab)
  if(strand == '+'){
    first_pos <- min(use_cds_tab$start)-min(exon_tab$start)+1
  }
  if(strand == '-'){
    first_pos <- max(exon_tab$end)-max(use_cds_tab$end)+1
  }
  CRNA_Start_pos <- Pos2RNA_exon_start(use_cds_tab$start,use_cds_tab$end,
                                                  first_pos,inner_size=inner_size,strand=strand)
  CRNA_End_pos <- Pos2RNA_exon_end(use_cds_tab$start,use_cds_tab$end,
                                              first_pos,inner_size=inner_size,strand=strand)
  ERNA_Start_pos <- Pos2RNA_exon_start(exon_tab$start,exon_tab$end,first_pos=0,
                                               inner_size=inner_size,strand=strand)
  ERNA_End_pos <- Pos2RNA_exon_end(exon_tab$start,exon_tab$end,first_pos=0,
                                           inner_size=inner_size,strand=strand)
  pf <- annotate_domain(transcript_tab=transcript_tab,exon_tab=exon_tab,
                        cds_tab=cds_tab,pf=pf,inner_size=inner_size)
  ## draw
  all_s <- unique(variant_info$Sample)
  all_v <- unique(variant_info$Variant)
  table_v <- table(variant_info$Variant)
  n <- nrow(exon_tab)
  ud <- unique(pf$Domain)
  stt = ERNA_Start_pos[1]
  end = ERNA_End_pos[n]
  mm <- min(pf$RNA_Start_pos)
  #  unit <- 0.35
  #  unit_s <- 0.4
  par(mar=c(1,2,1,4))
  cc <- brewer.pal(9,'Pastel1')
  if(use_color==FALSE) cc <- colorRampPalette(c('grey','black'))(20)
  plot_new(xlim=c(stt,end),ylim = c(-1.1-length(ud)*unit,
                                    max(table_v)*unit_s+2))
  rect(xleft=ERNA_Start_pos,
       xright=ERNA_End_pos,
       ybottom=-0.5,ytop=0.1,col=cc[1],border = 'black')
  if(inner_size>0){
    ww <- 2:length(ERNA_End_pos)
    segments(x0=ERNA_End_pos[ww-1],x1=ERNA_Start_pos[ww],
             y0=-0.2,y1=-0.2)
  }
  rect(xleft=CRNA_Start_pos[1],
       xright=max(CRNA_End_pos),
       ybottom=-0.8,ytop=-0.6,col=cc[2],border = cc[2])
  text(ERNA_Start_pos/2+ERNA_End_pos/2,y=-0.2,
       sprintf('E%s',1:length(ERNA_Start_pos)),cex=0.5)
  text(CRNA_Start_pos[1]/2+max(CRNA_End_pos)/2,y=-0.7,
       'CDS',cex=0.5)
  cc <- get.class.color(ud)
  if(use_color==FALSE) cc <- colorRampPalette(c('dark grey','black'))(20)
  for(i in 1:length(ud)){
    rect(xleft=pf$RNA_Start_pos[which(pf$Domain == ud[i])],
         xright=pf$RNA_End_pos[which(pf$Domain == ud[i])],
         ybottom=-0.9-unit*i,ytop=-0.9-unit*i+unit,col=cc[i],border = cc[i])
    text(x=mm,y=-0.9-unit*i+unit/2,pos=2,xpd=TRUE,ud[i],cex=0.5)
  }
  ###
  all_mut <- unique(variant_info$MutationType)
  cc1 <- get.class.color(all_mut,use_color = brewer.pal(8,'Dark2'))
  cc1 <- adjustcolor(cc1,0.7);names(cc1) <- all_mut
  if(use_color==FALSE){
    cc1 <- rep(1,length.out=length(cc1)); names(cc1) <- all_mut
  }
  yy <- seq(from=unit_s,by=unit_s,length.out=max(table_v))
  mm <- min(variant_info$CDS_pos)
  all_v2pos <- sapply(all_v,function(x)unique(variant_info[which(variant_info$Variant == x),]$CDS_pos))
  xx <- seq(min(variant_info$CDS_pos),max(variant_info$CDS_pos),
            length.out = length(all_v))+first_pos
  prev_pos <- 0
  ddd <- (end-stt)/300 ### minimum distance between adjacent points
  for(i in 1:length(all_v)){
    x1 <- variant_info[which(variant_info$Variant == all_v[i]),]
    if(i > 1){
      if(x1$CDS_pos[1]-prev_pos<ddd/2){
        x1$CDS_pos <- x1$CDS_pos+ddd
      }
    }
    #segments(x0=x1$CDS_pos+first_pos,x1=x1$CDS_pos+first_pos,
    #         y0=yy[1],y1=yy[nrow(x1)],lwd=2,col='light grey')
    segments(x0=x1$CDS_pos+first_pos,x1=x1$CDS_pos+first_pos,
             y0=yy[nrow(x1)],y1=max(yy)+unit_s,lwd=0.5,
             col='dark grey',lty=2)
    segments(x0=x1$CDS_pos+first_pos,x1=xx[i],
             y0=max(yy)+unit_s,y1=max(yy)+unit_s*2,lwd=0.5,
             col='dark grey',lty=1)
    text(x=xx[i],y=max(yy)+unit_s*2,
         all_v[i],adj=0,srt=90,cex=0.5)
    prev_pos <- x1$CDS_pos[1]
  }
  for(i in 1:length(all_v)){
    x1 <- variant_info[which(variant_info$Variant == all_v[i]),]
    if(i > 1){
      if(x1$CDS_pos[1]-prev_pos<ddd/2){
        x1$CDS_pos <- x1$CDS_pos+ddd
      }
    }
    ## mark_outcome
    if(is.null(mark_outcome)==FALSE){
      w1 <- which(x1$Outcome == mark_outcome)
      if(length(w1)>0){
        if(mark_type=='sample'){
          boxtext(x1$CDS_pos[w1],yy[w1],label='    ',border=2)
        }
        if(mark_type=='variant'){
          boxtext(xx[i],max(yy)+unit_s*2,label=all_v[i],
                  border=2,adj=0,srt=90,cex=0.5,add_per = 10)
        }
      }
    }
    if(use_color==TRUE){
      pch1 <- 17; pch2 <- 16;
      points(x=x1$CDS_pos+first_pos,
             y=yy[1:nrow(x1)],
             pch=ifelse(x1$Zygosity == 'Heterozygous',pch2,pch1),
             xpd=TRUE,col=cc1[x1$MutationType],
             cex=ifelse(x1$Zygosity == 'Heterozygous',0.8,1.1))
    }else{
      all_p1 <- c(0,1,2,5,6)[1:length(cc1)]; names(all_p1) <- names(cc1)
      all_p2 <- c(15,16,17,23,25)[1:length(cc1)]; names(all_p2) <- names(cc1)
      # ppch: 0 15; 1 16; 2 17; 5 23; 6 25;
      points(x=x1$CDS_pos+first_pos,
             y=yy[1:nrow(x1)],
             pch=ifelse(x1$Zygosity == 'Heterozygous',all_p1[x1$MutationType],
                        all_p2[x1$MutationType]),
             xpd=TRUE,col=1,
             cex=ifelse(x1$Zygosity == 'Heterozygous',0.8,1.1))
    }
    prev_pos <- x1$CDS_pos[1]
  }
  ## legend
  if(use_color==TRUE){
    all_zyg <- c('Heterozygous','Homozygous')
    mm <- max(ERNA_End_pos)
    width_label=names(cc1)[which.max(nchar(names(cc1)))]
    pp <- par()$usr
    yy1 <- seq(min(yy),pp[4],length.out=2+length(all_zyg)+length(cc1))
    for(i in 1:length(cc1)){
      boxtext(x=mm+strwidth(width_label,cex=0.5),
              #y=length(all_s)/2*unit_s-unit_s*i,
              y=yy1[i+length(all_zyg)],
              fill=cc1[i],label = names(cc1)[i],
              width_label=width_label,
              adj=0.5,cex=0.5,border=NA,add_per = 20)
    }
    legend(x=mm,y=yy1[1]/2+yy1[2]/2,pch=c(2,1),
           all_zyg,
           xjust=0,yjust=0,cex=0.5,xpd=T,border = NA,bty='n')
    if(is.null(mark_outcome)==FALSE){
      boxtext(x=mm+strwidth(width_label,cex=0.5),
              y=length(all_s)*unit_s*0.7,
              mark_outcome,cex=0.5,adj=0.5,
              width_label=width_label,
              col=2,border=2,add_per = 20)
    }
  }else{
    all_zyg <- c('Heterozygous','Homozygous')
    mm <- max(RNA_End_pos)
    width_label=names(cc1)[which.max(nchar(names(cc1)))]
    for(i in 1:length(cc1)){
      text(x=mm+strwidth(width_label,cex=0.5),
              y=length(all_s)/2*unit_s-unit_s*i,
              label = names(cc1)[i],
              adj=0,cex=0.5,xpd=T)
      points(x=mm+strwidth(width_label,cex=0.5)/2,
             y=length(all_s)/2*unit_s-unit_s*i,adj=1,pch=all_p1[i],xpd=T)
    }
    legend(x=mm,y=length(all_s)/2*unit_s,pch=c(all_p1[1],all_p2[1]),
           all_zyg,
           xjust=0,yjust=0,cex=0.5,xpd=T,border = NA,bty='n')
    if(is.null(mark_outcome)==FALSE){
      boxtext(x=mm+strwidth(width_label,cex=0.5),
              y=length(all_s)*unit_s*0.7,
              mark_outcome,cex=0.5,adj=0.5,
              width_label=width_label,
              col=2,border=2,add_per = 20)
    }
  }
}
pos_convert <- function(x,x_start,x_end,range_min,range_max){
  rr <- range_max-range_min
  (x-x_start)/(x_end-x_start+1)*rr+range_min
}

