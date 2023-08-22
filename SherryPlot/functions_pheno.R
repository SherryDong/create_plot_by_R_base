draw_phenoDist <- function(data,Name_col='System',Total_col='Total',Sub_col='Diagnosed',cex.label = 0.9,cex.legend=1.2,
                           color_name='Paired'){
  p2 <- data[,Total_col]
  names(p2) <- data[,Name_col]
  p3 <- sort(p2)
  p4 <- cumsum(p3)
  t2xy <- function(t,init.angle=0,radius=0.8) {
    twopi <- pi*2
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  tt  <- t2xy(c(1:c(max(p4)+max(p3)/2)/c(max(p3)/2+max(p4))))
  tt1 <- t2xy(c(1:c(max(p4)+max(p3)/2)/c(max(p3)/2+max(p4))),radius = 1)
  tt2 <- t2xy(c(1:c(max(p4)+max(p3)/2)/c(max(p3)/2+max(p4))),radius = 0.9)
  ptt <- cbind(tt$x[p4],tt$y[p4])
  ptt1 <- cbind(tt1$x[p4],tt1$y[p4]);
  ptt2 <- cbind(tt2$x[p4],tt2$y[p4]);
  rownames(ptt)<-names(p4);rownames(ptt1)<-names(p4);rownames(ptt2)<-names(p4)
  cc1 <- brewer.pal(8,'Pastel2')
  cc2 <- colorRampPalette(brewer.pal(8,'Reds'))(100)
  par(mar=c(3,3,3,3))
  plot(1,col='white',xaxt='n',yaxt='n',bty='n',xlab='',ylab='',xlim=c(-1,1),ylim=c(-1,1))
  w1 <- which(p3/max(p3)<0.3)
  w2 <- which(p3/max(p3)>=0.3)
  draw.ellipse(0,0,a=0.8,b=0.8,border=cc1[1])
  segments(x0=ptt2[w1,1],y0=ptt2[w1,2],x1=ptt[w1,1],y1=ptt[w1,2],col='grey')
  draw.ellipse(ptt,a=0.3*p3/max(p3),b=0.3*p3/max(p3),
               xpd=TRUE,col='grey',border='white')
  ## add pie
  rownames(data) <- data[,Name_col]
  if(length(Sub_col)>0){
    if(length(Sub_col)<=8) cc1 <- brewer.pal(8,color_name)[1:length(Sub_col)]
    if(length(Sub_col)>8) cc1 <- colorRampPalette(brewer.pal(8,color_name))(length(Sub_col))
    for(i in rownames(data)){
      xx_mid <- ptt[i,1]
      yy_mid <- ptt[i,2]
      ra <- 0.3*p3[i]/max(p3); rb <- 0.3*p3[i]/max(p3)
      s1 <- as.numeric(data[i,Sub_col])/data[i,Total_col]; s2 <- c(0,cumsum(s1));#print(s2)
      for(j in 1:length(Sub_col)){
        ss <- seq.int(s2[j],s2[j+1],length.out=100)
        xx <- ra*cos(2*pi*ss); yy <- rb*sin(2*pi*ss);
        polygon(c(xx_mid,xx+xx_mid,xx_mid),c(yy_mid,yy+yy_mid,yy_mid),col=cc1[j],border=NA,xpd=TRUE)
      }
    }
    legend(0,0,fill=cc1,Sub_col,xpd=T,border=NA,bty='n',cex=cex.legend,xjust = 0.5,yjust = 0.5)
  }
  ##
  text(ptt[w2,],names(p4)[w2],cex=cex.label,xpd=T)
  text(ptt1[w1,],names(p4)[w1],cex=cex.label,xpd=T)
  ##
}

draw_phenoDist_Two <- function(data,Name_col='System_name',Total_num=4214,Total_col='System',
                               Total_legend_text='Total',Sub_legend_text='Diagnosed',
                               Sub_num=294,legend_text='Patients with involvement of this organ',
                               Sub_col='System_Diagnosed',cex.label = 0.9,cex.legend=1.2,
                           color_name='Paired'){
  p2 <- rep(1,length.out=length(data[,Total_col]))
  names(p2) <- data[,Name_col]
  p3 <- sort(p2)
  p4 <- cumsum(p3)
  t2xy <- function(t,init.angle=0,radius=0.8) {
    twopi <- pi*2
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  tt  <- t2xy(c(0:(length(p2)-1))/length(p2))
  tt1  <- t2xy(c(0:(length(p2)-1))/length(p2),radius=1)
  ptt <- cbind(tt$x,tt$y);ptt1 <- cbind(tt1$x,tt1$y);
  rownames(ptt)<-names(p2);rownames(ptt1)<-names(p2);
  cc1 <- brewer.pal(8,'Pastel2')
  cc2 <- colorRampPalette(brewer.pal(8,'Reds'))(100)
  plot(1,col='white',xaxt='n',yaxt='n',bty='n',xlab='',ylab='',xlim=c(-1,1),ylim=c(-1,1))
  #w1 <- which(p3/max(p3)<0.3)
  #w2 <- which(p3/max(p3)>=0.3)
  draw.ellipse(0,0,a=0.8,b=0.8,border=cc1[1])
  draw.ellipse(ptt,a=0.15,b=0.15,
               xpd=TRUE,col='light grey',border='white')

  rownames(data) <- data[,Name_col]
  ## add pie
  cc1 <- brewer.pal(8,color_name)[1:2]
  ##
  if(length(Total_col)>0){
    for(i in rownames(data)){
      xx_mid <- ptt[i,1]
      yy_mid <- ptt[i,2]
      ra <- 0.15*p3[i]/max(p3); rb <- 0.15*p3[i]/max(p3)
      s1 <- as.numeric(data[i,Total_col])/Total_num; s2 <- c(0,cumsum(s1));#print(s2)
      for(j in 1:length(Total_col)){
        ss <- seq.int(s2[j],s2[j+1],length.out=100)
        xx <- ra*cos(2*pi*ss); yy <- rb*sin(2*pi*ss);
        polygon(c(xx_mid,xx+xx_mid,xx_mid),c(yy_mid,yy+yy_mid,yy_mid),col=cc1[j],border=NA,xpd=TRUE)
      }
    }
  }
  ##
  draw.ellipse(ptt,a=0.1,b=0.1,
               xpd=TRUE,col='white',border='white')
  draw.ellipse(ptt,a=0.07,b=0.07,
               xpd=TRUE,col='light grey',border='white')
  ##
  if(length(Sub_col)>0){
    for(i in rownames(data)){
      xx_mid <- ptt[i,1]
      yy_mid <- ptt[i,2]
      ra <- 0.07*p3[i]/max(p3); rb <- 0.07*p3[i]/max(p3)
      s1 <- as.numeric(data[i,Sub_col])/Sub_num; s2 <- c(0,cumsum(s1));#print(s2)
      for(j in 1:length(Sub_col)){
        ss <- seq.int(s2[j],s2[j+1],length.out=100)
        xx <- ra*cos(2*pi*ss); yy <- rb*sin(2*pi*ss);
        polygon(c(xx_mid,xx+xx_mid,xx_mid),c(yy_mid,yy+yy_mid,yy_mid),col=cc1[j],border=NA,xpd=TRUE)
      }
    }
  }
  ## legend
  draw.ellipse(0,0,a=0.15,b=0.15,xpd=TRUE,col='light grey',border='white')
  text(0,0.15,Total_legend_text,xpd=TRUE,cex=cex.legend*0.9)
  draw.ellipse(0,0,a=0.1,b=0.1,xpd=TRUE,col='white',border='white')
  draw.ellipse(0,0,a=0.07,b=0.07,xpd=TRUE,col='light grey',border='white')
  text(0,0,Sub_legend_text,xpd=TRUE,cex=cex.legend*0.9)
  legend(0,-0.25,fill=cc1,legend_text,xpd=T,border=NA,bty='n',cex=cex.legend,xjust = 0.5,yjust = 0.5)
  ##
  for(i in 1:nrow(ptt1)){
    #text(x=ptt1[i,1],y=ptt1[i,2],names(p4)[i],cex=cex.label,xpd=T,srt=180*atan(ptt1[i,2]/ptt1[i,1])/pi,adj=ifelse(ptt1[i,1]>0,0,1))
    text(x=ptt1[i,1],y=ptt1[i,2],names(p4)[i],cex=cex.label,xpd=T,srt=0,adj=ifelse(ptt1[i,1]>0,0,1))
  }
  ##
}
