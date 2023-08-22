##
library(quantsmooth)
##
draw_CNV <- function(data,Type_col,Sample_col,Chr_col,Start_col,End_col,pos_align='verticle',
                     cytoband_width=0.15,region_width=0.8,
                     each_width=0.02,cex_chr=1,draw_bandname=FALSE,cex_bandname=0.3,dup_del_col=brewer.pal(9,'Set1')[1:2]){
  CHR<-c(1:22,'X','Y')  # Chromosomes
  chrlen<-lengthChromosome(CHR,"bases")
  ##
  if(pos_align=='verticle'){
    chr_part1 <- c(1:12);
    chr_part2 <- c(13:22,'X','Y');
    chr_part1 <- as.character(chr_part1)
    chr_part2 <- as.character(chr_part2)
    chr_part  <- cbind(chr_part1,chr_part2)
    max_h <- max(apply(chr_part,1,function(x)chrlen[x[1]]+chrlen[x[2]]))*1.05
    dd <- max_h-max_h/1.05
    max_w <- max(length(chr_part1),length(chr_part2))
    #
    par(mar=c(4,2,4,4))
    plot(1,col='white',xlim=c(0,max_w),ylim=c(0,max_h),xaxt='n',yaxt='n',xaxs='i',yaxs='i',bty='n',xlab='',ylab='')
    for(i in 1:length(chr_part1)){
      paintCytobands(chr_part1[i],pos=c(i,max_h),width=cytoband_width,xpd=TRUE,
                     legend=draw_bandname,cex.leg = cex_bandname,orientation = 'v')
      text(i-cytoband_width/2,max_h,chr_part1[i],cex=cex_chr,xpd=T,pos=3)
    }
    for(i in 1:length(chr_part2)){
      paintCytobands(chr_part2[i],pos=c(i,chrlen[chr_part2[i]]),width=cytoband_width,xpd=TRUE,
                     legend=draw_bandname,cex.leg = cex_bandname,orientation = 'v')
      text(i-cytoband_width/2,0,chr_part2[i],cex=cex_chr,xpd=T,pos=1)
    }
    ##
    ## step for sample
#    region_width <- 1-cytoband_width-0.1
    cc <- dup_del_col
    cc <- rgb(t(col2rgb(cc)/255),alpha=0.8)
    data[,Chr_col] <- gsub('chr','',data[,Chr_col],ignore.case = T)
    ss <- (apply(cbind(data[,Chr_col],data[,Sample_col]),1,function(x)paste(x,collapse='_')))
    step <- rep(0.1,length=length(ss));names(step) <- paste0('chr',ss);
    for(i in CHR){
      w1 <- grep(paste0('^chr',i,'_'),names(step))
      w1 <- w1[order(table(ss)[names(step)[w1]], decreasing = TRUE)]
      w11 <- w1[grep('DUP',data[w1,Type_col])]
      w12 <- w1[grep('DEL',data[w1,Type_col])]
      step[w11] <- seq(cytoband_width/2+0.05,region_width/2-0.05,length.out=length(w11))
      step[w12] <- seq(cytoband_width/2+0.05,region_width/2-0.05,length.out=length(w12))
    }
    #
    chr2x_1 <- 1:length(chr_part1)-cytoband_width/2; names(chr2x_1) <- chr_part1
    chr2x_2 <- 1:length(chr_part2)-cytoband_width/2; names(chr2x_2) <- chr_part2
    chr2y_1 <- rep(max_h,length.out=length(chr_part1));names(chr2y_1) <- chr_part1
    chr2y_2 <- chrlen[chr_part2]
    chr2x <- c(chr2x_1,chr2x_2);chr2y <- c(chr2y_1,chr2y_2)
    #
    for(i in 1:nrow(data)){
      ch <- data[i,Chr_col]
      tt <- data[i,Type_col]
      nn <- step[i]
      ss <- data[i,Start_col]
      ee <- data[i,End_col]
      if(ss=='whole chromosome') ss <- 1
      if(ee=='whole chromosome') ee <- chrlen[ch]
      ss <- as.numeric(ss);ee <- as.numeric(ee);
      if(grepl('DUP',tt)){
        polygon(y=chr2y[ch]-c(ss,ee,ee,ss),
                x=chr2x[ch]-c(nn+each_width,nn+each_width,nn,nn),col=cc[1],border = NA,xpd=T)
      }else{
        polygon(y=chr2y[ch]-c(ss,ee,ee,ss),
                x=chr2x[ch]+c(nn+each_width,nn+each_width,nn,nn),col=cc[2],border = NA,xpd=T)
      }
    }
    pp <- par()$usr
    legend(pp[2],pp[3]-max_h*0.05,fill=cc[1:2],c('Duplication','Deletion'),border = NA,bty='n',xpd=T,horiz = T,cex=1.2,xjust=1,yjust=1)
  }
  if(pos_align=='horizon'){
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
    ##
    ## step for sample
#    region_width <- 1-cytoband_width-0.1
    cc <- dup_del_col
    cc <- rgb(t(col2rgb(cc)/255),alpha=0.8)
    data[,Chr_col] <- gsub('chr','',data[,Chr_col],ignore.case = T)
    ss <- (apply(cbind(data[,Chr_col],data[,Sample_col]),1,function(x)paste(x,collapse='_')))
    step <- rep(0.1,length=length(ss));names(step) <- paste0('chr',ss);
    for(i in CHR){
      w1 <- grep(paste0('^chr',i,'_'),names(step))
      w1 <- w1[order(table(ss)[names(step)[w1]], decreasing = TRUE)]
      step[w1] <- seq(cytoband_width/2+0.05,region_width-0.05,length.out=length(w1))
    }
    #
    chr2x_1 <- length(chr_part1):1-cytoband_width/2; names(chr2x_1) <- chr_part1
    chr2x_2 <- length(chr_part2):1-cytoband_width/2; names(chr2x_2) <- chr_part2
    chr2y_1 <- rep(0,length.out=length(chr_part1));names(chr2y_1) <- chr_part1
    chr2y_2 <- max_w-chrlen[chr_part2]
    chr2y <- c(chr2x_1,chr2x_2);chr2x <- c(chr2y_1,chr2y_2)
    #
    for(i in 1:nrow(data)){
      ch <- data[i,Chr_col]
      tt <- data[i,Type_col]
      nn <- step[i]
      ss <- data[i,Start_col]
      ee <- data[i,End_col]
      if(ss=='whole chromosome') ss <- 1
      if(ee=='whole chromosome') ee <- chrlen[ch]
      ss <- as.numeric(ss);ee <- as.numeric(ee);
      if(grepl('DUP',tt)){
        polygon(x=chr2x[ch]+c(ss,ee,ee,ss),
                y=chr2y[ch]-c(nn+each_width,nn+each_width,nn,nn),col=cc[1],border = NA,xpd=T)
      }else{
        polygon(x=chr2x[ch]+c(ss,ee,ee,ss),
                y=chr2y[ch]-c(nn+each_width,nn+each_width,nn,nn),col=cc[2],border = NA,xpd=T)
      }
    }
    pp <- par()$usr
    legend(pp[2],pp[3],fill=cc[1:2],c('Duplication','Deletion'),border = NA,bty='n',xpd=T,horiz = T,cex=1.2,xjust=1,yjust=1)
  }
##
}
draw_CNV_circos <- function(bone,data,Type_col,Sample_col,Chr_col,Mark_col,
                            Start_col,End_col,query=NULL,dup_del_col=brewer.pal(9,'Set1')[1:2]){
  #source('D:/analysis_eng/SherryPlot/function_circos.R')
  #source('D:/analysis_eng/SherryPlot/functions_genebody.R')
  dd <- 100
  max_radius <- 0.8
  min_radius <- 0.3
  bone$len <- (bone[,End_col]-bone[,Start_col])/dd
  bone$use_len <- bone$len
  thre_len <- max(bone$len)/nrow(bone) ## set minimun len to the 1/nrow(bone) maximum len
  bone$use_len[which(bone$len<thre_len)] <- thre_len
  total_len <- sum(bone$use_len)
  interval_len <- (total_len/10)/nrow(bone)
  total_len1 <- total_len+interval_len*nrow(bone)
  ## draw backbone
  plot(1,xlab='',ylab='',xaxt='n',yaxt='n',main='',
       xlim=c(-1,1),ylim=c(-1,1),bty='n',col='white')
  interval_p <- interval_len/total_len1
  bone$len_p <- bone$use_len/total_len1
  bone$end_p <- cumsum(bone$use_len)/total_len1+interval_p*(1:nrow(bone))
  bone$stt_p <- bone$end_p-bone$len_p
  start_pos <- bone$stt_p
  end_pos   <- bone$end_p
  bone$each_p <- bone$len_p/bone$len
  gc <- adjustcolor('light grey',0.5)
  gd <- adjustcolor('dark grey',1)
  for(i in 1:nrow(bone)){
    each_circle_circle(start_pos[i],end_pos[i],0.9,0.85,col=gd,
                       border = gd)
    each_circle_circle(start_pos[i],end_pos[i],max_radius,min_radius,
                       col=gc,
                       border = gd)
    x1 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.9)
    x2 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.95)
    each_text(x2$x[500],x2$y[500],sprintf('%s',bone[i,1]),style=1,cex=max_radius)
    each_text(x1$x[10],x1$y[10],sprintf('%s',get_label_manual(bone[i,Start_col])),
              style=2,cex=0.3)
    each_text(x1$x[1000-10],x1$y[1000-10],sprintf('%s',get_label_manual(bone[i,End_col])),
              style=2,cex=0.3)
    if(is.null(query)==FALSE){
      b1 <- bone[i,]
      q1 <- query[which(query[,Mark_col] == bone[i,Mark_col]),]
      if(nrow(q1)>0){
        for(xxx in 1:nrow(q1)){
          st <- b1$stt_p+b1$each_p*(q1[xxx,Start_col]-b1[,Start_col])/dd
          ed <- b1$stt_p+b1$each_p*(q1[xxx,End_col]-b1[,Start_col])/dd
          if(st<b1$stt_p) st = b1$stt_p
          if(ed>b1$end_p) ed = b1$end_p
          each_circle_circle(st,ed,0.88,0.87,col='black',
                             border = 'black')
          each_circle_circle(st,ed,max_radius,min_radius,col=gc,
                             border =gc)
          x1 <- t2xy(seq(st,ed,length.out = 1000),radius=0.865)
          each_text(x1$x[10],x1$y[10],sprintf('%s',get_label_manual(q1[xxx,Start_col])),
                    style=3,cex=0.25)
          each_text(x1$x[1000-10],x1$y[1000-10],sprintf('%s',get_label_manual(q1[xxx,End_col])),
                    style=3,cex=0.25)
        }
      }
    }
  }
  ## draw cnv
  cc <- dup_del_col
  cc <- rgb(t(col2rgb(cc)/255),alpha=0.8)
  names(cc) <- c('DUP','DEL')
  data_sort <- lapply(bone[,Mark_col],function(x){
    x1 <- data[which(data[,Mark_col]==x),]
    x1 <- x1[order(x1[,Start_col]+x1[,End_col]/10-x1[,Start_col]/10),]
    x2 <- lapply(unique(x1[,Sample_col]),function(xx){
      x1[which(x1[,Sample_col]==xx),]
    })
    names(x2) <- unique(x1[,Sample_col])
    x2
  })
  names(data_sort) <- bone[,Mark_col]
  chr_count <- unlist(lapply(data_sort,length))
  each_p_radius <- -(max_radius-min_radius)/max(chr_count)
  for(x in names(data_sort)){
    x1 <- data_sort[[x]]
    b1 <- bone[which(bone[,Mark_col]==x),]
    r1 <- seq(max_radius,by=each_p_radius,length.out=length(x1))
    r2 <- seq(max_radius+each_p_radius,by=each_p_radius,
              length.out=length(x1))-each_p_radius/3
    names(r1) <- names(r2) <- names(x1)
    for(xx in names(x1)){
      x2 <- x1[[xx]]
      for(xxx in 1:nrow(x2)){
        st <- b1$stt_p+b1$each_p*(x2[xxx,Start_col]-b1[,Start_col])/dd
        ed <- b1$stt_p+b1$each_p*(x2[xxx,End_col]-b1[,Start_col])/dd
        if(ed-st<0.0005){ed=st+0.0005}
        each_circle_circle(start_pos = st,
                           end_pos = ed,
                           radius1 = r1[xx],radius2 = r2[xx],
                           col=cc[x2[xxx,Type_col]],border=NA)
      }
    }
  }
  ##

}
