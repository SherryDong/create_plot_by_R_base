## create plot by R-base
library(RColorBrewer)
#' Create DNA
#'
#' @param DNA_length Unit length of DNA double helix.
#' @examples
#' draw.DNA(DNA_length=5)
#' @export
draw.DNA <- function(DNA_length=4){
  cc <- brewer.pal(11,'Set3') ## 2-5,1-4
  cc[5] <- brewer.pal(11,'Paired')[1]
  cc[1] <- brewer.pal(8,'Accent')[1]
  if(DNA_length%%2==0){
    x <- seq(-DNA_length*pi/2,DNA_length*pi/2,length.out=1000) ##
    y1 <- cos(x)
    y2 <- cos(x+pi)
    plot(y1~x,pch=16,type='l',xlab='',ylab='',xaxt='n',yaxt='n',main='',bty='n',lwd=6,col=brewer.pal(8,'Set1')[2])
    xx <- seq(DNA_length*pi/2,-DNA_length*pi/2,length.out = DNA_length*5+1); xx <- xx+(xx[2]-xx[1])/2
    xx <- setdiff(xx,c(xx[c(1:DNA_length)*5-2],min(xx)))
  }
  if(DNA_length%%2==1){
    x <- seq(-(DNA_length+1)*pi/2,(DNA_length-1)*pi/2,length.out=1000) ##
    y1 <- cos(x)
    y2 <- cos(x+pi)
    plot(y1~x,pch=16,type='l',xlab='',ylab='',xaxt='n',yaxt='n',main='',bty='n',lwd=6,col=brewer.pal(8,'Set1')[2])
    xx <- seq(-(DNA_length+1)*pi/2,(DNA_length-1)*pi/2,length.out = DNA_length*5+1); xx <- xx+(xx[2]-xx[1])/2
    xx <- setdiff(xx,max(xx))
  }
  up_rr <- c();
  for(i in 1:length(xx)){
    ybottom <- cos(xx[i])
    ytop    <- cos(xx[i]+pi)
    rr <- sample(1:4,1) ## ATCG
    if(rr==2){
      segments(y0=ybottom,y1=0,x0=xx[i],x1=xx[i],col=cc[2],lwd=4) ##
      segments(y0=0,y1=ytop,x0=xx[i],x1=xx[i],col=cc[5],lwd=4)
    }
    if(rr==4){
      segments(y0=ybottom,y1=0,x0=xx[i],x1=xx[i],col=cc[5],lwd=4) ##
      segments(y0=0,y1=ytop,x0=xx[i],x1=xx[i],col=cc[2],lwd=4)
    }
    if(rr==1){
      segments(y0=ybottom,y1=0,x0=xx[i],x1=xx[i],col=cc[1],lwd=4) ##
      segments(y0=0,y1=ytop,x0=xx[i],x1=xx[i],col=cc[4],lwd=4)
    }
    if(rr==3){
      segments(y0=ybottom,y1=0,x0=xx[i],x1=xx[i],col=cc[4],lwd=4) ##
      segments(y0=0,y1=ytop,x0=xx[i],x1=xx[i],col=cc[1],lwd=4)
    }
    up_rr[i] <- rr
  }
  lines(y1~x,pch=16,lwd=8,col=brewer.pal(8,'Set1')[2])
  lines(y2~x,pch=16,lwd=8,col=brewer.pal(8,'Set1')[2])
}

