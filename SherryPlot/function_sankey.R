curve_sankey <- function (x0, x1, y0, y1,nsteps = 1000){
  xx <- seq(-pi/2, pi/2, length.out = nsteps)
  yy <- y0 + (y1 - y0) * (sin(xx) + 1)/2
  xx <- seq(x0, x1, length.out = nsteps)
  list(x=xx,y=yy)
}
sankey_polygon <- function(Fx0, Fx1, Fy0, Fy1,
                           Sx0, Sx1, Sy0, Sy1,
                           nsteps = 100,col='red',border=NA,...){
  curve_pos1 <- curve_sankey(Fx0,Fx1,Fy0,Fy1)
  curve_pos2 <- curve_sankey(Sx0,Sx1,Sy0,Sy1)
  polygon(x=c(curve_pos1$x,rev(curve_pos2$x)),
          y=c(curve_pos1$y,rev(curve_pos2$y)),
          col=col,border=border,xpd=TRUE,...)
}
#
plot_trace <- function(all_mat,top_n=3,each_width=0.25,
                       remove_char='.',
                       count_link_thre=0.1,part_thre=0.1,
                       group_info=NULL,
                       only_use_group=TRUE,
                       group_name=NULL){
  n <- ncol(all_mat)
  cc <- colorRampPalette(brewer.pal(8,'Spectral'))(top_n+1)
  if(is.null(group_info)==FALSE){
    par(mar=c(5,3,ifelse(is.null(group_name)==TRUE,2,5),5))
  }else{
    par(mar=c(5,3,ifelse(is.null(group_name)==TRUE,2,5),2))
  }
  plot_new(xlim=c(1.1,n+0.25),ylim=c(0,1))
  hap_start <- list(); hap_end <- list(); all_top <- list();
  all_valid <- list()
  for(i in 1:n){
    w1 <- which(all_mat[,i]!=remove_char)
    all_valid[[i]] <- w1
    tmp1 <- sort(table(all_mat[w1,i]),decreasing = T)
    tmp1 <- tmp1/sum(tmp1)
    if(is.null(group_info)==FALSE){
      w1 <- unlist(lapply(group_info,function(x)x[[i]]))
      if(only_use_group==TRUE){
        tmp1 <- tmp1[c(w1)]
      }else{
        tmp1 <- tmp1[c(w1,setdiff(names(tmp1),w1))]
      }
    }
    if(length(tmp1)>top_n){
      tmp1 <- c(tmp1[1:top_n],other=sum(tmp1)-sum(tmp1[1:top_n]))
      top_k <- top_n+1
    }else{
      top_k <- length(tmp1)
    }
    tmp2 <- cumsum(rev(tmp1))
    tmp31 <- c(0,tmp2[1:(top_k-1)]); tmp32 <- c(tmp2[1:top_k])
    hap_start[[i]] <- rev(tmp31);
    hap_end[[i]]   <- rev(tmp32);
    all_top[[i]] <- names(tmp1)
    for(j in 1:top_k){
      rect(xleft=i,xright = i+each_width,
           ybottom = tmp31[j],ytop=tmp32[j],col=cc[top_k+1-j],xpd=TRUE)
    }
    legend(x=i,y=-0.05,names(tmp1),rev(cc[top_k+1-c(1:top_k)]),
           cex=0.5,xpd=T,border = NA,bty='n',yjust = 1)
  }
  if(is.null(group_name)==FALSE){
    text(1:n+each_width/2,1.05,group_name,srt=90,adj = 0,xpd=TRUE,cex=0.7)
  }
  mat1 <- as.matrix(all_mat)
  for(i in 1:(n-1)){
    start_left <- rep(0,length(all_top[[i]]))
    end_right <- rep(0,length(all_top[[i+1]]))
    for(j1 in 1:length(all_top[[i]])){
      for(j2 in 1:length(all_top[[i+1]])){
        inter_count <- length(which(mat1[,i]==all_top[[i]][j1] & 
                                      mat1[,i+1]==all_top[[i+1]][j2]))
        count_link <- inter_count/length(intersect(all_valid[[i]],all_valid[[i+1]]))
        per1_part <- inter_count/length(which(mat1[,i]==all_top[[i]][j1]))
        per1 <- inter_count/length(all_valid[[i]])
        per2 <- inter_count/length(all_valid[[i+1]])
        if(count_link > count_link_thre & per1_part>part_thre){
          Fy0 <- hap_end[[i]][j1]-start_left[j1]
          Fy1 <- hap_end[[i+1]][j2]-end_right[j2]
          Sy0 <- Fy0-per1
          Sy1 <- Fy1-per2
          sankey_polygon(i+each_width,i+1,Fy0,Fy1,i+each_width,i+1,Sy0,Sy1,
                         col = adjustcolor(cc[j1],0.3),
                         border = 1,lwd=0.4)
          start_left[j1] <- start_left[j1]+per1
          end_right[j2] <- end_right[j2]+per2
        }
      }
    }
  }
  if(is.null(group_info)==FALSE){
    legend(x=n+0.5,y=0.5,names(group_info),
           fill=cc[1:length(names(group_info))],xpd=T,cex=0.75,
           border = NA,bty='n')
  }
  axis(side=2,at=c(0:10)/10,labels = sprintf('%s%s',c(0:10)*10,'%'),
       las=2,cex.axis=0.7,cex=0.7)
}