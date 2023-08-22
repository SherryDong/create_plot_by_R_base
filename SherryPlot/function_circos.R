library(NetBID2)
library(RCircos)
library(quantsmooth)
library(RCircos)
library(annotables) #
# Human build 38 (grch38)
# Human build 37 (grch37)
# Mouse (grcm38)
# Rat (rnor6)
# Chicken (galgal5)
# Worm (wbcel235)
# Fly (bdgp6)
data(UCSC.HG19.Human.CytoBandIdeogram);
source('D:/analysis_eng/SherryPlot/function_basic.R')
##########
t2xy <- function(t,init.angle=90,radius=1) {
  t2p <- 2*pi * t + init.angle * pi/180
  list(x = -radius * cos(t2p), y = radius * sin(t2p))
}
twopi <-  2 * pi;
##########
draw_circos_each <- function(x1=1,x2=1.2,r=1,h=0.1,col='red',tag='',tag_pos=1.5,tag_adj=1,tag_cex=1,init.angle=90,use_tag=TRUE,border=NA){
  P1 <- t2xy(seq.int(x1,x2,length.out = 100),r,init.angle=init.angle)
  P2 <- t2xy(seq.int(x1,x2,length.out = 100),r+h,init.angle=init.angle)
  polygon(c(P1$x,rev(P2$x)), c(P1$y,rev(P2$y)), density = NA, col=col,border=border)
  P3 <- t2xy((x1+x2)/2,r+h*tag_pos,init.angle=init.angle)
  if(use_tag==TRUE){
    text(P3$x,P3$y,tag,xpd=TRUE,adj=tag_adj,cex=tag_cex)
  }else{
    return(P3)
  }
}
adj_pos <- function(x,s=range(x)[1],e=range(x)[2]){
  seq(from=s,to=e,length.out=length(x))
}
adj_pos_auto <- function(x,min_w=0.0025){
  x_ori <- x
  if(length(x)==1) return(x)
  dd <- x[2:length(x)]-x[1:(length(x)-1)]
  w1 <- which(dd<min_w)
  if(length(w1)==0) return(x)
  for(i in 2:length(x)){
    d1 <- x[i]-x[i-1]
    if(d1<min_w){
      x[i] <- x[i-1]+min_w
    }
  }
  return(x)
}
adj_pos_auto_2 <- function(x,min_w=0.0025){
  dd <- x[2:length(x)]-x[1:(length(x)-1)]
  x1 <- x
  w1 <- c()
  d1 <- dd[1]
  for(i in 1:length(dd)){
    print(x[w1])
    if(d1<min_w*(length(w1)+1)){
      print(c(names(dd)[i],d1,i,length(w1)+1))
      w1 <- unique(c(w1,c(i-1,i)))
      d1 <- d1+dd[i]
    }else{
      if(length(w1)>1){
        print(x1[w1])
        x1[w1] <- adj_pos_auto_each(x1[w1],min_w=min_w)
      }
      w1 <- c();d1 <- dd[i];
    }
  }
  return(x1)
}
adj_pos_auto_each <- function(x,min_w=0.0025){
  m_x <- mean(x)
  s <- m_x-floor(length(x)/2)*min_w
  e <- m_x+floor(length(x)/2)*min_w
  seq(from=s,to=e,length.out=length(x))
}
###
bleacher <- function(x) {
  bleach <- 0
  (x * (1 - bleach))
}
RCircos.Link.Line<-function (line.start = NULL, line.end = NULL)
{
  P0 <- line.start
  P2 <- line.end
  numOfBCPoint <- 1000
  t <- seq(0, 1, 1/numOfBCPoint)
  linkX <- (1 - t)^2 * P0[1] + t^2 * P2[1]
  linkY <- (1 - t)^2 * P0[2] + t^2 * P2[2]
  return(list(pos.x = linkX, pos.y = linkY))
}
each_circle_circle <- function(start_pos,end_pos,radius1,radius2,col,border){
  x1 <- t2xy(seq(start_pos,end_pos,length.out = 1000),radius=radius1)
  x2 <- t2xy(seq(start_pos,end_pos,length.out = 1000),radius=radius2)
  polygon(x=c(x1$x,rev(x2$x)),y=c(x1$y,rev(x2$y)),border=border,
          col=col,lwd=0.2)
}
each_text <- function(x,y,text,style=1,cex=1){
  if(style==1){
    ss <- ifelse(x>0,-90,90)
    text(x=x,y=y,text,xpd=TRUE,
         srt=ss+180*(atan(y/x))/pi,adj=0.5,cex=cex)
  }
  if(style==2){
    ss <- ifelse(x>0,0,1)
    text(x=x,y=y,text,xpd=TRUE,srt=180*(atan(y/x))/pi,
         adj=ss,cex=cex)
  }
  if(style==3){
    ss <- ifelse(x>0,1,0)
    text(x=x,y=y,text,xpd=TRUE,srt=180*(atan(y/x))/pi,
         adj=ss,cex=cex)
  }
}
draw_circos <- function(tab,ov_mat,main='',ov_text=TRUE){
  ####
  plot(1,col='white',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',
       xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),main=main)
  all_pos <- cumsum(tab)/sum(tab)
  start_pos <- c(0,all_pos[-length(all_pos)])
  end_pos <- all_pos
  cc2chr <- colorRampPalette(brewer.pal(9,'Spectral'))(length(tab));
  names(cc2chr) <- names(tab)
  for(i in 1:length(tab)){
    x1 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=1)
    x2 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.8)
    polygon(x=c(x1$x,rev(x2$x)),y=c(x1$y,rev(x2$y)),border='white',
            col=cc2chr[names(tab)[i]])
    x1 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.9)
    x2 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=1.1)
    text(x=x2$x[500],y=x2$y[500],names(tab)[i],xpd=TRUE,
         srt=180*(atan(mean(x2$y)/mean(x2$x)))/pi,adj=ifelse(mean(x1$x)>0,0,1),cex=1.2)
    text(x=x1$x[500],y=x1$y[500],tab[i],xpd=TRUE,
         srt=180*(atan(mean(x1$y)/mean(x1$x)))/pi,adj=0.5,cex=1)
  }
  ## interaction : f1_info_noS
  x0 <- t2xy(seq(0,1,length.out = 1001),radius=0.8)
  ov_start_pos <- start_pos
  ov_end_pos <- end_pos
  for(i in 1:(nrow(ov_mat)-1)){
    for(j in (i+1):nrow(ov_mat)){
      v = ov_mat[i,j]
      if(v>0){
        v1 <- v/sum(tab)
        s1 <- 1+round(1000*(ov_start_pos[i]))
        s2 <- 1+round(1000*(ov_start_pos[j]))
        s3 <- 1+round(1000*(ov_start_pos[i]+v1))
        s4 <- 1+round(1000*(ov_start_pos[j]+v1))
        x1 <- RCircos.Link.Line(c(x0$x[s1],x0$y[s1]),c(x0$x[s4],x0$y[s4]))
        x2 <- RCircos.Link.Line(c(x0$x[s3],x0$y[s3]),c(x0$x[s2],x0$y[s2]))
        polygon(x=c(x1$pos.x,rev(x0$x[s2:s4]),rev(x2$pos.x),rev(x0$x[s1:s3])),
                y=c(x1$pos.y,rev(x0$y[s2:s4]),rev(x2$pos.y),rev(x0$y[s1:s3])),border='white',
                col=colorRampPalette(cc2chr[names(tab)[i:j]])(3)[2],xpd=T)
        ov_start_pos[i] <- ov_start_pos[i]+v1
        ov_start_pos[j] <- ov_start_pos[j]+v1
        text(x=x0$x[mean(s2:s4)],y=x0$y[mean(s2:s4)],v)
        text(x=x0$x[mean(s1:s3)],y=x0$y[mean(s1:s3)],v)
        mid_x <- x1$pos.x[length(x1$pos.x)/2]/2+x2$pos.x[length(x2$pos.x)/2]/2
        mid_y <- x1$pos.y[length(x1$pos.y)/2]/2+x2$pos.y[length(x2$pos.y)/2]/2
        if(ov_text==TRUE) text(mid_x,mid_y,sprintf('%s+%s',names(tab)[i],names(tab)[j]))
      }
    }
  }
  ##
}
##
draw_circos_complicate <- function(tab,ov_mat,ori_mat,main='',ov_text=TRUE,
                                   fig_lim=1.1,use_cairo=FALSE){
  ####
  plot(1,col='white',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',
       xlim=c(-fig_lim,fig_lim),ylim=c(-fig_lim,fig_lim),main=main)
  all_pos <- cumsum(tab)/sum(tab)
  start_pos <- c(0,all_pos[-length(all_pos)])
  end_pos <- all_pos
  cc2chr <- colorRampPalette(brewer.pal(9,'Spectral'))(length(tab));
  names(cc2chr) <- names(tab)
  for(i in 1:length(tab)){
    x1 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=1)
    x2 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.8)
    polygon(x=c(x1$x,rev(x2$x)),y=c(x1$y,rev(x2$y)),border='white',
            col=cc2chr[names(tab)[i]])
    x1 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.9)
    x2 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=1.1)
    if(use_cairo==TRUE){
      text(x=x2$x[500],y=x2$y[500],names(tab)[i],xpd=TRUE,
           srt=180*(atan(mean(x2$y)/mean(x2$x)))/pi,adj=ifelse(mean(x1$x)>0,0,1),
           cex=1.2,family="SimHei")
    }else{
      text(x=x2$x[500],y=x2$y[500],names(tab)[i],xpd=TRUE,
           srt=180*(atan(mean(x2$y)/mean(x2$x)))/pi,adj=ifelse(mean(x1$x)>0,0,1),cex=1.2)
    }
    text(x=x1$x[500],y=x1$y[500],tab[i],xpd=TRUE,
         srt=180*(atan(mean(x1$y)/mean(x1$x)))/pi,adj=0.5,cex=1)
  }
  ## interaction : f1_info_noS
  cate2ind <- apply(ori_mat,2,function(x)rownames(ori_mat)[which(x==1)])
  use2ind <- lapply(rep('',length.out=length(cate2ind)),c);
  names(use2ind) <- names(cate2ind)
  x0 <- t2xy(seq(0,1,length.out = 1001),radius=0.8)
  ov_start_pos <- start_pos
  ov_end_pos <- end_pos
  for(i in 1:(nrow(ov_mat)-1)){
    for(j in (i+1):nrow(ov_mat)){
      v = ov_mat[i,j]
      if(v>0){
        v1 <- v/sum(tab)
        dd <- intersect(cate2ind[[i]],cate2ind[[j]])
        diffv_i <- setdiff(setdiff(use2ind[[i]],dd),'')
        diffv_j <- setdiff(setdiff(use2ind[[j]],dd),'')
        ov_start_pos[i] <- length(diffv_i)/sum(tab)+ov_start_pos[i]
        ov_start_pos[j] <- length(diffv_j)/sum(tab)+ov_start_pos[j]
        s1 <- 1+round(1000*(ov_start_pos[i]))
        s2 <- 1+round(1000*(ov_start_pos[j]))
        s3 <- 1+round(1000*(ov_start_pos[i]+v1))
        s4 <- 1+round(1000*(ov_start_pos[j]+v1))
        x1 <- RCircos.Link.Line(c(x0$x[s1],x0$y[s1]),c(x0$x[s4],x0$y[s4]))
        x2 <- RCircos.Link.Line(c(x0$x[s3],x0$y[s3]),c(x0$x[s2],x0$y[s2]))
        polygon(x=c(x1$pos.x,rev(x0$x[s2:s4]),rev(x2$pos.x),rev(x0$x[s1:s3])),
                y=c(x1$pos.y,rev(x0$y[s2:s4]),rev(x2$pos.y),rev(x0$y[s1:s3])),
                border='white',
                col=adjustcolor(colorRampPalette(cc2chr[names(tab)[i:j]])(3)[2],0.5),xpd=T)
        text(x=x0$x[mean(s2:s4)],y=x0$y[mean(s2:s4)],v)
        text(x=x0$x[mean(s1:s3)],y=x0$y[mean(s1:s3)],v)
        mid_x <- x1$pos.x[length(x1$pos.x)/2]/2+x2$pos.x[length(x2$pos.x)/2]/2
        mid_y <- x1$pos.y[length(x1$pos.y)/2]/2+x2$pos.y[length(x2$pos.y)/2]/2
        if(ov_text==TRUE) text(mid_x,mid_y,sprintf('%s+%s',names(tab)[i],names(tab)[j]))
        dd <- intersect(cate2ind[[i]],cate2ind[[j]])
        use2ind[[i]] <- dd
        use2ind[[j]] <- dd
      }
    }
  }
  ##
}
curve_circos <- function (x0, x1, y0, y1,nsteps = 100){
  tt <- seq(0, 1, 1/nsteps)
  xx <- (1 - tt)^2 * x0 + tt^2 * x1
  yy <- (1 - tt)^2 * y0 + tt^2 * y1
  list(x=xx,y=yy)
}
####
draw_circos_phenoGeno <- function(dis2gene,gene_bed=WES_gene,
                                  radius_chr=1.4,radius_system=1.4,radius_gene=1,
                                  class1='pathogenic',class2='VUS',
                                  gene_cex=0.8,num_cex=0.8,system_cex=1.4,legend_col=NULL,
                                  count_col=NULL,sort_disease=T){
  all_disease <- table(dis2gene[,1])
  if(sort_disease==T) all_disease <- sort(all_disease)
  all_gene <- unique(dis2gene[,2])
  ##
  CHR<-c(1:22,'X','Y')  # Chromosomes
  MapInfo<-lengthChromosome(CHR,"bases")/1000
  total_len <- sum(MapInfo)
  ##
  gene_bed <- gene_bed[which(gene_bed$V4 %in% all_gene),]
  gene2pos <- lapply(all_gene,function(x){
    x1 <- gene_bed[which(gene_bed$V4 %in% x),]
    c(x1[1,1],(max(c(x1$V2,x1$V3))+min(c(x1$V2,x1$V3)))/2000,x)
  })
  gene2pos <- do.call(rbind,gene2pos)
  gene2pos <- gene2pos[order(gene2pos[,1]),]
  gene2pos[,1] <- gsub('chr','',gene2pos[,1])
  rownames(gene2pos) <- gene2pos[,3]
  ### deal with gene too close
  #min_len <- 1000
  #gene2pos_tmp <- as.data.frame(gene2pos,stringsAsFactors=F)
  #gene2pos_tmp$V2 <- as.numeric(gene2pos_tmp$V2)
  ###############################################
  col1 <- colorRampPalette((brewer.pal(11,'Spectral')))(24)
  col2 <- rev(colorRampPalette((brewer.pal(10,'Paired')))(length(all_disease)))
  names(col1) <- CHR
  names(col2) <- names(all_disease)
  ##
  x <- MapInfo
  x <- c(0, cumsum(x)/sum(x));
  x <- x*0.7
  xw <- x[2:length(x)]-x[1:(length(x)-1)]
  x1 <- x[1:(length(x)-1)]
  x2 <- x1+xw-0.0025
  names(x1) <- names(x)[2:length(x)]
  names(x2) <- names(x1)
  CHR_START <- x1
  CHR_END <- x2
  gene2pos_r <- apply(gene2pos,1,function(x){
    CHR_START[x[1]]+as.numeric(x[2])/MapInfo[x[1]]*(CHR_END[x[1]]-CHR_START[x[1]])
  })
  UCSC.HG19.Human.CytoBandIdeogram$Chromosome <- as.character(UCSC.HG19.Human.CytoBandIdeogram$Chromosome)
  UCSC.HG19.Human.CytoBandIdeogram$Chromosome <- gsub("chr","",UCSC.HG19.Human.CytoBandIdeogram$Chromosome)
  # acen    gneg gpos100  gpos25  gpos50  gpos75    gvar   stalk
  cyto2pos_r <- t(apply(UCSC.HG19.Human.CytoBandIdeogram,1,function(x){
    s <- CHR_START[x[1]]+as.numeric(x[2])/(1000*MapInfo[x[1]])*(CHR_END[x[1]]-CHR_START[x[1]])
    e <- CHR_START[x[1]]+as.numeric(x[3])/(1000*MapInfo[x[1]])*(CHR_END[x[1]]-CHR_START[x[1]])
    c(s,e,x[5],x[1])
  }))
  names(gene2pos_r) <- gene2pos[,3]
  gene2pos_r <- sort(gene2pos_r)
  w1 <- names(which(gene2pos[names(gene2pos_r),1]=='X'));gene2pos_r[w1] <- adj_pos(gene2pos_r[w1],s=CHR_START['X'],e=CHR_END['X'])
  #for(i in 1:10){
    gene2pos_r<- adj_pos_auto(gene2pos_r,min_w=0.003)
  #}
  ########
  xlim <- ylim <- c(-1.5, 1.5)
  plot_new(xlim,ylim)
  ######## draw chromosome
  for(i in 1:length(x1)){
    cc <- col1[i]
    cc1 <- rgb(t(col2rgb(cc))/255,alpha=0.1)
    draw_circos_each(x1=-x1[i],x2=-x2[i],r=radius_chr,h=0.1,
                     col=cc1,tag=names(x1)[i],tag_pos=1.8,tag_adj=0,init.angle=90,border=1)
  }
  type.b <- match(cyto2pos_r[,3], c("acen", "gneg", "gpos", "gvar", "stalk", "gpos25", "gpos50", "gpos75", "gpos100","gpos33", "gpos66"))
  bandcol <- bleacher(c(0.5, 1, 0.2, 0.6, 0.75, 0.7,0.5, 0.3, 0.1, 0.6, 0.4))[type.b]
  for(i in 1:nrow(cyto2pos_r)){
    cc <- col1[cyto2pos_r[i,4]]
    cc <- colorRampPalette(c('white',cc,'black'))(300)[60:160]
    cc1 <- cc[101-bandcol[i]*100]
    if(cyto2pos_r[i,3] == 'acen') cc1<-'red'
    draw_circos_each(x1=-as.numeric(cyto2pos_r[i,1]),x2=-as.numeric(cyto2pos_r[i,2]),
                     r=radius_chr,h=0.1,col=cc1,use_tag=FALSE,border=NA)
  }
  ######## draw gene
  gene2pos_p <- t2xy(gene2pos_r,radius=radius_gene,init.angle=90) # gene
  gene2pos_pp <- t2xy(gene2pos_r,radius=radius_gene-0.2,init.angle=90)
  gene2pos_ppp <- t2xy(gene2pos_r,radius=radius_gene-0.15,init.angle=90) # count
  if(is.null(count_col)==T){
    gene2pos_p <- t2xy(gene2pos_r,radius=radius_gene-0.15,init.angle=90) # gene
  }
  for(i in 1:length(gene2pos_r)){
    angle <- atan(-gene2pos_p$y[i]/gene2pos_p$x[i])*180/pi
    if(-gene2pos_p$x[i]>0) adj <- 0 else adj <- 1
    text(-gene2pos_p$x[i],gene2pos_p$y[i],names(gene2pos_r)[i],xpd=TRUE,
         cex=gene_cex,srt=angle,adj=adj)
    w1 <- which(dis2gene$Gene==names(gene2pos_r)[i])
    w2 <- dis2gene[w1,,drop=F];print(w2)
    if(class1 != '' & class2 != ''){
      c1 <- length(which(w2$classification==class1))
      c2 <- length(which(w2$classification==class2))
      if(is.null(count_col)==F){
        c1 <- sum(as.numeric(w2[which(w2$classification==class1),count_col]))
        c2 <- sum(as.numeric(w2[which(w2$classification==class2),count_col]))
      }
      c3 <- paste0('(',c1,',',c2,')')
    }else{
      c1 <- nrow(w2)
      if(is.null(count_col)==F){
        c1 <- sum(as.numeric(w2[,count_col]))
      }
      c3 <- paste0('(',c1,')')
    }
    if(is.null(count_col)==T){c3<-''}
    text(-gene2pos_ppp$x[i],gene2pos_ppp$y[i],c3,xpd=TRUE,
         cex=num_cex,srt=angle,adj=adj)
    points(-gene2pos_pp$x[i],gene2pos_pp$y[i],pch=16,col=1,cex=0.12)
  }
  ######## draw disease
  x <- all_disease
  x <- c(0, cumsum(x)/sum(x));
  x <- x*0.28+0.01
  xw <- x[2:length(x)]-x[1:(length(x)-1)]
  x1 <- x[1:(length(x)-1)]
  x2 <- x1+xw-0.0025
  names(x1) <- names(x)[2:length(x)]
  ########
  init.angle <- 90;
  for(i in 1:length(x1)){
    P3 <- draw_circos_each(x1=x1[i],x2=x2[i],r=radius_system,
                           h=0.1,col=col2[i],tag=names(x1)[i],
                           tag_pos=2,tag_adj=1,tag_cex=0.8,use_tag=FALSE)
    P1 <- t2xy((x1[i]+x2[i])/2,r=radius_system+0.1,init.angle=init.angle)
    if(i==1) {text(P3$x+0.1,P3$y+0.1,names(x1)[i],xpd=TRUE,adj=0,cex=system_cex);
      segments(x0=P3$x+0.05,y0=P3$y+0.05,x1=P1$x,y1=P1$y,xpd=TRUE);}
    if(i>=2) {text(P3$x,P3$y,names(x1)[i],xpd=TRUE,adj=0,cex=system_cex);
      segments(x0=P3$x+0.01,y0=P3$y,x1=P1$x,y1=P1$y,xpd=TRUE);}
  }
  disease2pos_pp <- t2xy((x1+x2)/2,radius=radius_system-0.1,init.angle=90)
  for(i in 1:length(gene2pos_r)){
    points(disease2pos_pp$x[i],disease2pos_pp$y[i],pch=16,col=1,cex=0.2)
  }
  ######## draw lines
  for(i in 1:nrow(dis2gene)){
    pp <- RCircos.Link.Line(line.start=c(disease2pos_pp$x[dis2gene[i,1]],disease2pos_pp$y[dis2gene[i,1]]),
                            line.end=c(-gene2pos_pp$x[dis2gene[i,2]],gene2pos_pp$y[dis2gene[i,2]]))
    lines(pp$pos.x,pp$pos.y,lwd=2,xpd=TRUE,col=col2[dis2gene[i,1]])
  }
  ##
  # add legend
  if(is.null(legend_col)==FALSE){
    #col2 <- rev(colorRampPalette((brewer.pal(10,'Paired')))(length(all_disease)))
    #names(col2) <- names(all_disease)
    tmp1 <- dis2gene
    all_c <- unique(tmp1[,legend_col])
    pp <- par()$usr;
    yy <- seq(pp[3],length.out = length(all_c),
              by=strheight('W')*(8+max(table(tmp1[,legend_col]))));
    names(yy) <- all_c
    tmp1 <- unique(dis2gene[,c('Disease.category',legend_col)]);
    for(i in all_c){
      w1 <- which(tmp1[,legend_col]==i)
      legend(x=pp[2],y=yy[i],xjust=0,yjust=0,
             legend=tmp1[w1,,drop=F]$Disease.category,
             fill=col2[tmp1[w1,,drop=F]$Disease.category],title=i,xpd=T,border=NA,bty='n')
    }
  }
}
###
draw_circos_genoRel <- function(gene_rel,pre_define,count_thre,WES_gene){
  cc2chr <- colorRampPalette(brewer.pal(9,'Spectral'))(24);
  names(cc2chr) <- sprintf("chr%s",c(1:22,'X','Y'))
  gene_rel <- unique(gene_rel)
  use_gene <- unique(c(gene_rel$gene1,gene_rel$gene2))
  use_gene <- intersect(use_gene,WES_gene$V4)
  gene_rel <- gene_rel[which(gene_rel$gene1 %in% use_gene & gene_rel$gene2 %in% use_gene),]
  tmp1 <- apply(gene_rel[,c(3,4,2)],1,function(x)paste(x,collapse=';'))
  tmp2 <- data.frame(do.call(rbind,lapply(names(table(tmp1)),function(x)strsplit(x,';')[[1]])),
                     Count=as.numeric(table(tmp1)),stringsAsFactors = FALSE)
  colnames(tmp2) <- c('gene1','gene2','group','count')
  gene_rel <- tmp2
  gene_rel <- gene_rel[which(gene_rel$count>=count_thre),]
  use_gene <- unique(c(gene_rel$gene1,gene_rel$gene2))
  use_gene_pos <- WES_gene[which(WES_gene$V4 %in% use_gene),]
  use_gene_pos <- use_gene_pos[order(use_gene_pos$V1),]
  use_gene_pos$count <- table(use_gene)[use_gene]
  use_gene_pos$pos <- cumsum(use_gene_pos$count)/sum(use_gene_pos$count)
  rownames(use_gene_pos) <- use_gene_pos$V4
  ##
  plot_new(xlim=c(-1.1,1.1),ylim=c(-1.1,1.1))
  start_pos <- c(0,use_gene_pos$pos[-nrow(use_gene_pos)])
  end_pos <- use_gene_pos$pos
  names(start_pos) <- names(end_pos) <- use_gene_pos[,4]
  f_g2p<-list()
  for(i in names(start_pos)){
    x1 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=1)
    x2 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.9)
    f_g2p[[use_gene_pos[i,4]]] <- c(mean(x2$x),mean(x2$y))
    polygon(x=c(x1$x,rev(x2$x)),y=c(x1$y,rev(x2$y)),border='white',
            col=cc2chr[use_gene_pos[i,1]])
    x1 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=1.02)
    x2 <- t2xy(seq(start_pos[i],end_pos[i],length.out = 1000),radius=0.925)
    text(x=mean(x2$x),y=mean(x2$y),use_gene_pos[i,1],xpd=TRUE,
         srt=180*(atan(mean(x2$y)/mean(x2$x)))/pi,adj=ifelse(mean(x1$x)>0,0,1),cex=0.5)
    text(x=mean(x1$x),y=mean(x1$y),i,xpd=TRUE,
         srt=180*(atan(mean(x1$y)/mean(x1$x)))/pi,adj=ifelse(mean(x1$x)>0,0,1),cex=0.8)
  }
  cc <- get.class.color(unique(gene_rel$group),pre_define=pre_define)
  for(i in 1:nrow(gene_rel)){
    l1 <- RCircos.Link.Line(f_g2p[[gene_rel[i,1]]],f_g2p[[gene_rel[i,2]]])
    if(is.na(l1$pos.x[1])==TRUE) print(str(l1))
    if(gene_rel[i,1]==gene_rel[i,2]){gene_rel[i,4]=gene_rel[i,4]*2}
    lines(l1$pos.x,l1$pos.y,col=cc[gene_rel[i,3]],lwd=1+(gene_rel[i,4]-1)*6)
  }
  legend('bottomright',fill=cc,legend=names(cc),xpd=TRUE,border=NA,bty='n',cex=0.9)
  #########
}
##
RCircos.Tile.Plot.mod <-
function (tile.data = NULL, track.num = NULL, side = c("in",
                                                       "out"), inside.pos = NULL, outside.pos = NULL, genomic.columns = 3,
          is.sorted = TRUE)
{
  if (is.null(tile.data))
    stop("Genomic data missing in RCircos.Tile.Plot().\n")
  if (is.null(genomic.columns) || genomic.columns < 3)
    stop(paste("Genomic position must include chromosome name,",
               "start and end position.\n"))
  boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos,
                                        outside.pos, FALSE)
  outerPos <- boundary[1]
  innerPos <- boundary[2]
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  tile.data <- RCircos.Get.Paired.Points.Positions(tile.data,
                                                   genomic.columns, plot.type = "tile")
  tile.layers <- RCircos.Get.Plot.Layers(tile.data, genomic.columns)
  layer.height <- RCircos.Par$track.height/RCircos.Par$max.layers
  num.layers <- max(tile.layers)
  if (num.layers > RCircos.Par$max.layers) {
    if (side == "in") {
      innerPos <- outerPos - layer.height * num.layers
    }
    else {
      outerPos <- innerPos + layer.height * num.layers
    }
    message(paste("Tiles plot may use more than one track.",
                  "Please select correct area for next track if necessary.\n"))
  }
  if (num.layers < RCircos.Par$max.layers) {
    layer.height <- RCircos.Par$track.height/num.layers
  }
  tile.colors <- RCircos.Get.Plot.Colors(tile.data, RCircos.Par$tile.color)
  RCircos.Track.Outline.mod(outerPos, innerPos, num.layers,lwd=0.5)
  the.loc <- ncol(tile.data)
  for (a.row in seq_len(nrow(tile.data))) {
    tile.start <- tile.data$LinkStart[a.row]
    tile.end <- tile.data$LinkEnd[a.row]
    layer.bot <- innerPos + layer.height * (tile.layers[a.row] -
                                              1)
    layer.top <- layer.bot + layer.height * 0.8
    polygon.x <- c(RCircos.Pos[tile.start:tile.end, 1] *
                     layer.top, RCircos.Pos[tile.end:tile.start, 1] *
                     layer.bot)
    polygon.y <- c(RCircos.Pos[tile.start:tile.end, 2] *
                     layer.top, RCircos.Pos[tile.end:tile.start, 2] *
                     layer.bot)
    polygon(polygon.x, polygon.y, col = tile.colors[a.row],lwd=0.2)
  }
}
RCircos.Track.Outline.mod <- function (inside.pos = NULL, outside.pos = NULL, num.layers = 1,
                                       chrom.list = NULL, track.colors = NULL,lwd=0.5)
{
  if (is.null(outside.pos) || is.null(inside.pos))
    stop("Missing outside.pos/inside.pos in RCircos.Track.Outline().\n")
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  subtrack.height <- (outside.pos - inside.pos)/num.layers
  chromosomes <- unique(as.character(RCircos.Cyto$Chromosome))
  if (!is.null(chrom.list)) {
    if (sum(chrom.list %in% chromosomes) != length(chrom.list)) {
      stop(paste("One or more chromosome is not",
                 "in chromosome ideogram data.\n"))
    }
    chromosomes <- chrom.list
  }
  if (is.null(track.colors)) {
    track.colors <- rep(RCircos.Par$track.background, length(chromosomes))
  }
  else {
    if (length(track.colors) != length(chromosomes))
      track.colors <- rep(track.colors, length(chromosomes))
  }
  for (aChr in seq_len(length(chromosomes))) {
    chr.rows <- which(RCircos.Cyto$Chromosome == chromosomes[aChr])
    the.chr <- RCircos.Cyto[chr.rows, ]
    plot.start <- min(RCircos.Cyto$StartPoint[chr.rows])
    plot.end <- max(RCircos.Cyto$EndPoint[chr.rows])
    polygon.x <- c(RCircos.Pos[plot.start:plot.end, 1] *
                     outside.pos, RCircos.Pos[plot.end:plot.start, 1] *
                     inside.pos)
    polygon.y <- c(RCircos.Pos[plot.start:plot.end, 2] *
                     outside.pos, RCircos.Pos[plot.end:plot.start, 2] *
                     inside.pos)
    polygon(polygon.x, polygon.y, col = track.colors[aChr],lwd=lwd)
    if (num.layers > 1) {
      for (a.line in seq_len(num.layers - 1)) {
        height <- outside.pos - a.line * subtrack.height
        lines(RCircos.Pos[plot.start:plot.end, 1] * height,
              RCircos.Pos[plot.start:plot.end, 2] * height,
              col = RCircos.Par$grid.line.color,lwd=lwd)
      }
    }
  }
}
