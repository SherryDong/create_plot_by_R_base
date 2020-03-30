library(RColorBrewer)
pdf('web.pdf',width=3,height=3)
par(mar=c(1,1,1,1))
t2xy <- function(t,init.angle=0,radius=sqrt(2)) {
  t2p <- pi*2 * t + init.angle * pi/180
  list(x = radius * cos(t2p), y = radius * sin(t2p))
}
cc <- brewer.pal(8,'Paired')
plot(1,xlim=c(-1,1),ylim=c(-1,1),
     bty='n',xlab='',ylab='',xaxt='n',yaxt='n',col='white',xaxs='i',yaxs='i')
draw.ellipse(x=0,y=0,a=1,b=1,lwd=14,border=cc[2],col='white',xpd=T)
segments(x0=-1,x1=1,y0=0,y1=0,lwd=8,col=cc[2])
segments(x0=0,x1=0,y0=-1,y1=1,lwd=8,col=cc[2])
xx <- 0.75
p1 <- t2xy(seq((pi-atan(1/xx))/(2*pi),(pi+atan(1/xx))/(2*pi),length.out = 1000),
           radius=sqrt(1+xx*xx))
lines(p1$x+xx,p1$y,xpd=F,lwd=8,col=cc[2])
lines(-p1$x-xx,p1$y,xpd=F,lwd=8,col=cc[2])
##
yy <- 1.75
p1 <- t2xy(seq((1.5*pi-asin(1/yy))/(2*pi),(1.5*pi+asin(1/yy))/(2*pi),length.out = 1000),
           radius=sqrt(yy*yy-1))
lines(p1$x,p1$y+yy,xpd=F,lwd=8,col=cc[2])
lines(p1$x,-p1$y-yy,xpd=F,lwd=8,col=cc[2])
dev.off()

