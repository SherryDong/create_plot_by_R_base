# grouplabel: e.g gene name
# groupid: e.g id
# 
#pdf('figure2.pdf',width=10,height=5.5)
#par(mar=c(1,1,2,2))
source('D:/analysis_eng/SherryPlot/function_basic.R')
###### only deal with class=2, if multiple need to revise !
draw_twoClass_fiveColRel <- function(dat,groupclass,groupx,groupy,groupid,grouplabel,
                            thre = 3, radius_in = 0.1,radius_out = 0.45,
                            circle_size_adj=4,
                            cexlabel=0.35,cexx=0.6,cexy=0.6,cexpv=0.5,cexgroup=0.8,
                            cexaxis=0.5,cextotal=0.6,
                            classcol=brewer.pal(9,'Pastel1')[c(2,4)],
                            xycol=adjustcolor(brewer.pal(11,'Set3'),0.5)[1]){
  d1 <- data.frame(class=dat[,groupclass],
                   x=dat[,groupx],
                   y=dat[,groupy],
                   id=dat[,groupid],
                   label=dat[,grouplabel],stringsAsFactors = F);
  u1 <- apply(d1,2,function(x)sort(unique(x)))
  u1$y <- u1$y[length(u1$y):1]
  ##
  nx <- length(u1$x)
  ny <- length(u1$y)
  c1 <- classcol;#brewer.pal(9,'Pastel1')[c(2,4,1,3,5:9)][1:length(u1$class)]
  c2 <- xycol;#adjustcolor(brewer.pal(11,'Set3'),0.5)
  u1 <- apply(d1,2,function(x)sort(unique(x)))
  u1$y <- u1$y[length(u1$y):1]
  ################################################
  plot_new(xlim=c(0,nx+2),ylim=c(0,ny+2))
  #rect(xleft=1,xright=nx+1,ybottom=1,ytop = ny+1)
  #segments(x0=1:(nx+1),x1=1:(nx+1),y0=1,y1=ny+1,col='blue')
  #segments(x0=1,x1=nx+1,y0=1:(ny+1),y1=1:(ny+1),col='blue')
  rect(xleft=1:nx,xright=2:(nx+1),ybottom = 0.4,ytop=0.8,border='grey',col=c2)
  rect(xleft=0.2,xright=0.8,ybottom = 1:ny,ytop=2:(ny+1),border='grey',col=c2)
  text(x=c(1:nx)+0.5,y=0.6,u1$x,xpd=T,cex=cexx,font=2)
  text(y=c(1:ny)+0.5,x=0.5,u1$y,xpd=T,cex=cexy,font=2)
  ####
  mm <- max(table(d1$x,d1$y))
  for(i in 1:nx){
    for(j in 1:ny){
      w1 <- which(d1$x==u1$x[i] &
                    d1$y==u1$y[j])
      t2 <- table(d1$`label`[w1],d1$class[w1])
      r1 <- t2[which(rowSums(t2)>0),,drop=F]
      if(nrow(r1)==0) next
      ## merge to other
      r2 <- rowSums(r1);r1 <- r1[order(r2,decreasing = T),,drop=F]
      w21 <- which(r2<=thre);w22 <- which(r2>thre)
      if(length(w22)>0 & length(w21)>0){
        r3 <- rbind(r1[names(w22),,drop=F],
                    'other'=colSums(r1[names(w21),,drop=F]))
      }else{r3<-r1;}
      r3 <- r3[order(rowSums(r3),decreasing = T),,drop=F]
      ## plot
      nr <- nrow(r3)+1; 
      dd1 <- sum(r3)/5
      if(nrow(r3)==1) dd1 <- 0
      dd2 <- 1/(sum(r3)+dd1)
      cr3 <- cumsum(rowSums(r3)+dd1/nrow(r3))
      nrs <- c(0,cr3[1:(length(cr3)-1)])*dd2
      nrs_start <- nrs;
      nrs_end <- nrs_start+rowSums(r3)*dd2;
      nrs_end_each <- nrs_start+r3[,1]*dd2;
      r_in <- rep(radius_in,length.out=nrow(r3))
      #r_out <- rep(radius_out,
      #             length.out=nrow(r3))
      r_out <- rep(radius_out*(sum(r3)/mm)^(1/circle_size_adj),
                   length.out=nrow(r3))
      #(radius_out-radius_in)*rowSums(r3)/mm+radius_in
      midx <- i+0.5; midy <- j+0.5
      ## polygon for gene
      for(ii in 1:nrow(r3)){
        p1 <- t2xy(seq(nrs_start[ii],nrs_end[ii],
                       length.out=100),
                   radius=r_in[ii])
        p2 <- t2xy(seq(nrs_start[ii],nrs_end[ii],
                       length.out=100),
                   radius=r_out[ii])
        polygon(x=midx+c(p1$x,rev(p2$x)),
                y=midy+c(p1$y,rev(p2$y)),
                col=c1[1],
                border = 'light grey')
      }
      ## polygon for CCGT
      for(ii in 1:nrow(r3)){
        p1 <- t2xy(seq(nrs_start[ii],nrs_end_each[ii],
                       length.out=100),
                   radius=r_in[ii])
        p2 <- t2xy(seq(nrs_start[ii],nrs_end_each[ii],
                       length.out=100),
                   radius=r_out[ii])
        polygon(x=midx+c(p1$x,rev(p2$x)),
                y=midy+c(p1$y,rev(p2$y)),
                col=c1[2],
                border = NA)
      }
      ## gene label
      for(ii in 1:nrow(r3)){
        p2 <- t2xy(seq(nrs_start[ii],nrs_end[ii],
                       length.out=100),
                   radius=r_out[ii])
        p3 <- t2xy(seq(nrs_start[ii],nrs_end[ii],
                       length.out=100),
                   radius=radius_out)
        mw <- length(p2$y)/2;
        mx <- median(p3$x)/2+median(p2$x)/2
        text(x=midx+mx,
             y=midy+p2$y[mw],
             rownames(r3)[[ii]],cex=cexlabel,
             xpd=T,adj=ifelse(mx>0,0.4,0.6))
      }
      ## total number
      text(midx,midy,sum(r3),adj=0.5,cex=cextotal,col=2)
    }
  }
  # barplot for zygosity
  t0 <- table(d1$y,d1$class)[u1$y,];
  pv <- signif(p.adjust(unlist(lapply(1:nrow(t0),function(x){
    fisher.test(rbind(t0[x,],colSums(t0)-t0[x,]))$p.value
  })),'fdr'),2)
  t1 <- t(t(t0)/colSums(t0))
  mm <- max(t1)
  ystart1 <- c(1:ny)-0.25+0.5
  yend1 <- c(1:ny)-0.025+0.5
  ystart2 <- c(1:ny)+0.025+0.5
  yend2 <- c(1:ny)+0.25+0.5
  xstart <- 1+nx+0.2
  for(i in 1:ny){
    if(t1[i,1]>0){
      rect(xleft=xstart,xright=1+nx+0.2+0.9*t1[i,1]/mm,
           ybottom = ystart1[i],ytop=yend1[i],col=c1[2])
    }
    if(t1[i,2]>0){
      rect(xleft=xstart,xright=1+nx+0.2+0.9*t1[i,2]/mm,
           ybottom = ystart2[i],ytop=yend2[i],col=c1[1])
    }
    mx <- max(1+nx+0.2+0.9*t1[i,]/mm)+0.05
    segments(x0=mx,x1=mx,y0=i-0.2+0.5,y1=i+0.2+0.5)
    if(pv[i]<=0.001)text(y=i+0.5,x=mx+0.15,sprintf('P=%s***',pv[i]),cex=cexpv,xpd=T,adj=0)
    if(pv[i]<=0.01&pv[i]>0.001)text(y=i+0.5,x=mx+0.15,sprintf('P=%s**',pv[i]),cex=cexpv,xpd=T,adj=0)
    if(pv[i]<=0.05&pv[i]>0.01)text(y=i+0.5,x=mx+0.15,sprintf('P=%s*',pv[i]),cex=cexpv,xpd=T,adj=0)
    if(pv[i]>0.05)text(y=i+0.5,x=mx+0.15,'ns',cex=cexpv,xpd=T,adj=0)
    ##
  }
  # legend
  segments(y0=1,y1=ny+1,x0=nx+1.2,x1=nx+1.2)
  segments(x0=nx+1.2,x1=nx+2.1,y0=1,y1=1)
  ss <- seq(nx+1.2,nx+2.1,length.out=5)
  bb <- round(seq(0,mm,length.out=5)*1000)/10
  segments(x0=ss,x1=ss,y0=1,y1=0.95)
  text(x=ss,y=0.9,sprintf('%s%s',bb,'%'),
       srt=45,adj=1,xpd=T,cex=cexaxis)
  
  # barplot for IUIS class
  t0 <- table(d1$x,d1$class)[u1$x,];
  pv <- signif(p.adjust(unlist(lapply(1:nrow(t0),function(x){
    fisher.test(rbind(t0[x,],colSums(t0)-t0[x,]))$p.value
  })),'fdr'),2)
  t1 <- t(t(t0)/colSums(t0))
  mm <- max(t1)
  ystart1 <- c(1:nx)-0.25+0.5
  yend1 <- c(1:nx)-0.025+0.5
  ystart2 <- c(1:nx)+0.025+0.5
  yend2 <- c(1:nx)+0.25+0.5
  xstart <- 1+ny+0.2
  for(i in 1:nx){
    if(t1[i,1]>0){
      rect(ybottom=xstart,ytop=1+ny+0.2+0.9*t1[i,1]/mm,
           xleft = ystart1[i],xright=yend1[i],col=c1[2],xpd=T)
    }
    if(t1[i,2]>0){
      rect(ybottom=xstart,ytop=1+ny+0.2+0.9*t1[i,2]/mm,
           xleft = ystart2[i],xright=yend2[i],col=c1[1],xpd=T)
    }
    # < 0.05 是 *，<0.01 是 **， <0.001 是 ***
    my <- max(1+ny+0.2+0.9*t1[i,]/mm)+0.05
    segments(y0=my,y1=my,x0=i-0.2+0.5,x1=i+0.2+0.5)
    if(pv[i]<=0.001)text(x=i+0.5,y=my+0.15,sprintf('P=%s***',pv[i]),cex=cexpv,xpd=T,adj=0.5)
    if(pv[i]<=0.01&pv[i]>0.001)text(x=i+0.5,y=my+0.15,sprintf('P=%s**',pv[i]),cex=cexpv,xpd=T,adj=0.5)
    if(pv[i]<=0.05&pv[i]>0.01)text(x=i+0.5,y=my+0.15,sprintf('P=%s*',pv[i]),cex=cexpv,xpd=T,adj=0.5)
    if(pv[i]>0.05)text(x=i+0.5,y=my+0.15,'ns',cex=cexpv,xpd=T,adj=0.5)
  }
  ## legend
  # legend
  segments(x0=1,x1=nx+1,y0=ny+1.2,y1=ny+1.2)
  segments(y0=ny+1.2,y1=ny+2.1,x0=1,x1=1)
  ss <- seq(ny+1.2,ny+2.1,length.out=5)
  bb <- round(seq(0,mm,length.out=5)*1000)/10
  segments(y0=ss,y1=ss,x0=1,x1=0.95)
  text(y=ss,x=0.9,sprintf('%s%s',bb,'%'),
       srt=0,adj=1,xpd=T,cex=cexaxis)
  ##
  legend(x=nx+1,y=ny+2,u1$class,
         fill=rev(c1),border=NA,bty='n',cex=cexgroup,
         xpd=T)
}