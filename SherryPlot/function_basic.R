
plot_new <- function(xlim=c(),ylim=c(),main='',...){
  plot(1,xlim=xlim,ylim=ylim,
       xaxt='n',yaxt='n',xlab='',ylab='',bty='n',col='white',main=main,...)
}
hp_prepare <- function(){
  hp_info <- read.delim('D:/analysis_eng/SherryPlot/data/hp_info.txt',stringsAsFactors = F,header = F)
  rownames(hp_info) <- hp_info$V1
  return(hp_info)
}
hp2list <- function(dat,add_des=FALSE){
  lapply(dat,function(x){
    unique(unlist(lapply(strsplit(x,';'),function(x1){
      x2 <- gsub('(HP:\\d\\d\\d\\d\\d\\d\\d).*','\\1',x1)
      if(add_des==TRUE) x2 <- sprintf('%s(%s)',x2,hp_info[x2,2])
    })))
  })
}
draw_image <- function(mat,cc,bb,add_line=TRUE){
  image(mat,bty='n',xaxt='n',yaxt='n',col=cc,breaks=bb);
  pp <- par()$usr
  xx <- seq(pp[1],pp[2],length.out=1+dim(mat)[1]);
  xx1 <- xx[1:(length(xx)-1)]/2+xx[2:length(xx)]/2
  yy <- seq(pp[3],pp[4],length.out=1+dim(mat)[2])
  yy1 <- yy[1:(length(yy)-1)]/2+yy[2:length(yy)]/2
  if(add_line==TRUE){
    segments(x0=pp[1],x1=pp[2],y0=yy,y1=yy,lwd=0.5)
    segments(x0=xx,x1=xx,y0=pp[3],y1=pp[4],lwd=0.5)
  }
  dxx <- xx1[2]-xx1[1]
  dyy <- yy1[2]-yy1[1]
  return(list(xx=xx,yy=yy,xx1=xx1,yy1=yy1,dxx=dxx,dyy=dyy,pp=pp))
}
t2xy <- function(t,init.angle=90,radius=1) {
  t2p <- 2*pi * t + init.angle * pi/180
  list(x = -radius * cos(t2p), y = radius * sin(t2p))
}
draw_filletRect <- function(xleft, ybottom, xright, ytop,
                            fillet_len=0.1,
                            density=NULL,angle=45,border=NULL,col=NA,lty= par("lty"),
                            N=100,
                            ...){
  #
  if(fillet_len==0){
    rect(xleft, ybottom, xright, ytop, density = density, angle = angle,
         border = border, col = col, lty = lty,...)
    return()
  }
  all_circosXY <- t2xy(seq(0,1,length.out=4*N),init.angle=0,radius=fillet_len)
  # left,top,right,bottom
  x_val <- c(rep(xleft,length.out=N),
             all_circosXY$x[1:N]+xleft+fillet_len,
             seq(xleft+fillet_len,xright-fillet_len,length.out = N),
             all_circosXY$x[c(N+1):c(2*N)]+xright-fillet_len,
             rep(xright,length.out=N),
             all_circosXY$x[c(2*N+1):c(3*N)]+xright-fillet_len,
             seq(xright-fillet_len,xleft+fillet_len,length.out = N),
             all_circosXY$x[c(3*N+1):c(4*N)]+xleft+fillet_len)
  y_val <- c(seq(ybottom+fillet_len,ytop-fillet_len,length.out = N),
             all_circosXY$y[1:N]+ytop-fillet_len,
             rep(ytop,length.out=N),
             all_circosXY$y[c(N+1):c(2*N)]+ytop-fillet_len,
             seq(ytop-fillet_len,ybottom+fillet_len,length.out = N),
             all_circosXY$y[c(2*N+1):c(3*N)]+ybottom+fillet_len,
             rep(ybottom,length.out=N),
             all_circosXY$y[c(3*N+1):c(4*N)]+ybottom+fillet_len)
  polygon(x=x_val, y = y_val, 
          density = density, angle = angle,
          border = border, col = col, lty = lty,
          ...)
}
draw_pie <- function(x,enter=T,...){
  x1 <- sprintf('%s%s',round(x/sum(x)*10000)/100,'%')
  if(enter==T) pie(x,labels=sprintf('%s\n(%s)',names(x),x1),...)
  if(enter==F) pie(x,labels=sprintf('%s(%s)',names(x),x1),...)
}
#plot_new(xlim=c(0,2),ylim=c(0,2))
#draw_filletRect(0,0,1,1,fillet_len=0.2)
##
boxtext <- function(x,y,label,width_label=label,fill='white',border='black',
                    adj=0.5,srt=0,cex=1,font=1,add_per=0,col=1){
  pp <- par()$usr
  ratio_x2y <- (pp[4]-pp[3])/(pp[2]-pp[1])
  charW <- strwidth(width_label,cex=cex,font=font);
  dW <- charW*(add_per/100)
  charH <- strheight(width_label,cex=cex,font=font)*1.5*(1+add_per/100);
  if(srt==90){
    charH <- strwidth('W',cex=cex,font=font)*1.5*(1+add_per/100);
    #charW <- strwidth(width_label,cex=cex,font=font)*ratio_x2y
    charW <- strheight(width_label,cex=cex,font=font)*nchar(width_label)
    dW <- charW*(add_per/100);
  }
  an <- srt/180*pi
  w1 <- c((dW+(1-adj)*charW)*cos(an),(dW+(1-adj)*charW)*sin(an))
  w2 <- c(charH/2*sin(an),-charH/2*cos(an))
  w3 <- c((dW+charW*adj)*cos(an),(dW+charW*adj)*sin(an))
  xx <- c(x+w1[1]+w2[1],x+w1[1]-w2[1],x-w3[1]-w2[1],x-w3[1]+w2[1])
  yy <- c(y+w1[2]+w2[2],y+w1[2]-w2[2],y-w3[2]-w2[2],y-w3[2]+w2[2])
  polygon(x=xx,y=yy,col=fill,border=border,xpd=TRUE)
  text(x=x,y=y,label=label,adj=adj,srt=srt,cex=cex,font=font,col=col)
}
#plot_new(xlim=c(-1,1),ylim=c(-1,1))
#par(mar=c(1,1,1,1))
#boxtext(0,0,'test',srt=90,cex=2,adj=0.5,add_per=30)
#boxtext(0,0,'test',srt=0,cex=2,adj=0,add_per=30)
#boxtext(0,0,'test',srt=90,cex=2,adj=0.5,add_per=30)
#boxtext(0,0,'test',srt=90,cex=2,adj=0,add_per=30)
# simple functions
list2mat <- function(input_list,all_x=NULL){
  if(is.null(all_x)==TRUE) all_x <- base::unique(unlist(input_list))
  all_y <- base::unique(names(input_list))
  mat1 <- matrix(0,nrow=base::length(all_x),ncol=base::length(all_y))
  rownames(mat1) <- all_x; colnames(mat1) <- all_y;
  for(i in names(input_list)){
    mat1[input_list[[i]],i] <- 1
  }
  return(mat1)
}
# prod
list2mat_value <- function(input_list,name_col,value_col,
                           val_strategy,default_value=NA){
  all_x <- base::unique(unlist(lapply(input_list,function(x)x[,name_col])))
  all_y <- base::unique(names(input_list))
  mat1 <- matrix(default_value,nrow=base::length(all_x),ncol=base::length(all_y))
  rownames(mat1) <- all_x; colnames(mat1) <- all_y;
  for(i in names(input_list)){
    mat1[input_list[[i]][,name_col],i] <- apply(input_list[[i]][,value_col,drop=F],1,function(x){get(val_strategy)(as.numeric(x))})
  }
  return(mat1)
}

list2df_narrow <- function(x){
  x1 <- unlist(lapply(names(x),function(xx)rep(xx,length.out=length(x[[xx]]))))
  x2 <- unlist(x)
  data.frame(x1,x2)
}
vec2list <- function(input_v,sep=NULL){
  if(is.null(sep)==TRUE){
    tmp2 <- list()
    input_vn <- names(input_v)
    input_v <- clean_charVector(input_v); names(input_v) <- input_vn
    for(i in 1:base::length(input_v)){
      if(input_v[i] %in% names(tmp2)){
        tmp2[[input_v[i]]] <- c(tmp2[[input_v[i]]],names(input_v)[i])
      }else{
        tmp2[[input_v[i]]] <- names(input_v)[i]
      }
    }
  }else{
    tmp1 <- stats::aggregate(names(input_v),list(input_v),function(x)base::paste(x,collapse=sep))
    tmp2 <- tmp1$x; names(tmp2) <- tmp1$Group.1
  }
  tmp2
}
pie_mod <- function (xmid,ymid,x,radius = 0.8, edges = 200,
                     clockwise = FALSE, col=NULL,
          init.angle = if (clockwise) 90 else 0,
          density = NULL, angle = 45,
          border = NULL, lty = NULL)
{
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  pin <- par("pin")
  if (is.null(col))
    col <- if (is.null(density))
      c("white", "lightblue", "mistyrose",
        "lightcyan", "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col))
    col <- rep_len(col, nx)
  if (!is.null(border))
    border <- rep_len(border, nx)
  if (!is.null(lty))
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density))
    density <- rep_len(density, nx)
  twopi <- ifelse (clockwise,-2 * pi,2 * pi)
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x+xmid, xmid), c(P$y+ymid, ymid),
            density = density[i], angle = angle[i],
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
  }
}
clean_table <- function(x){
  x1 <- lapply(x,function(x){
    xx <- gsub('\\s+$','',x)
    xx <- gsub('^\\s+','',xx)
    xx <- gsub('&#10;','',xx)
    xx <- gsub('\\s+$','',xx)
    xx <- gsub('^\\s+','',xx)
    xx <- gsub('\n','',xx)
  })
  if(class(x)=='data.frame'){
    x1 <- as.data.frame(x1)
  }
  return(x1)
}
df2list <- function(x,sepC='',remainC=colnames(x)){
  tmp1 <- lapply(unique(x[,sepC]),function(xx){
    x[which(x[,sepC]==xx),remainC]
  })
  names(tmp1) <- unique(x[,sepC]);return(tmp1)
}
get_label_manual <- function(x){
  x1 <- sapply(x,function(x2){
    x3 <- unlist(strsplit(as.character(x2),""))
    x4 <- base::length(x3)%/%3 ## add number
    if(x4>0){
      pp <- base::length(x3)-base::seq(1,x4)*3; x3[pp] <- paste0(x3[pp],','); base::paste(x3,collapse="")
    }else{
      x2
    }
  })
  unlist(x1)
}
par.pos2inch<-function(){
  user.range <- par("usr")[c(2,4)] - par("usr")[c(1,3)]
  region.pin <- par("pin")
  return(region.pin/user.range)
}
add_m<-function(x,log='',cex=2,xadjust=TRUE){
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
draw_Fisher <- function(x1,x2,x3,p,OR){
  r1 <- 0.75
  plot_new(xlim=c(-0.5,2.5),ylim=c(-0.5,1))
  rect(xleft=-0.5,ybottom = -0.5,xright=2.5,ytop=1)
  draw.ellipse(x=0.8,y=0.5,a=0.3,b=0.25,xpd=T)
  draw.ellipse(x=1.2,y=0.5,a=0.3/r1,b=0.25/r1,xpd=T)
  text(0.5/2+(1.2-0.3/r1)/2,0.5,x2,adj=0.5,cex=0.8)
  text(1.1/2+(1.2+0.3/r1)/2,0.5,x3,adj=0.5,cex=0.8)
  text(1.1/2+(1.2-0.3/r1)/2,0.5,x1,adj=0.5,cex=0.8)
  text(1,-0.15,sprintf('P=%s\nOdds ratio=%s',signif(p,3),signif(OR,3)),adj=0.5,cex=0.8)
  text(0.1,0.5,'Inferred\ntargets',cex=0.8)
  text(1.4+0.3/r1,0.5,'ChIP-seq',cex=0.8,adj=0)
}

