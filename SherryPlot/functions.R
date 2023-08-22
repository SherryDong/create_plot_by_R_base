library(RColorBrewer)
twopi <-  2 * pi;
##########
##########
t2xy <- function(t,radius=1,init.angle=0) {
  t2p <- twopi * t + init.angle * pi/180
  list(x = radius * cos(t2p), y = radius * sin(t2p))
}
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
  dd <- x[2:length(x)]-x[1:(length(x)-1)]
  x1 <- x
  w1 <- c()
  d1 <- dd[1]
  for(i in 1:length(dd)){
    print(x[w1])
    if(d1<min_w*(length(w1)+1)){
      #print(c(names(dd)[i],d1,i,length(w1)+1))
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
  (x * (1 - bleach))
}
get_line <- function(line.start = NULL, line.end = NULL){
  P0 <- line.start
  P2 <- line.end
  t <- seq(0, 1, 0.001)
  linkX <- (1 - t)^2 * P0[1] + t^2 * P2[1]
  linkY <- (1 - t)^2 * P0[2] + t^2 * P2[2]
  return(list(pos.x = linkX, pos.y = linkY))
}





