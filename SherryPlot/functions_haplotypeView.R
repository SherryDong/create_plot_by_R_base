source('D:/analysis_eng/SherryPlot/function_basic.R')
# Family ID ('FID')
# Within-family ID (sample ID) ('IID'; cannot be '0')
# Within-family ID of father (Paternal ID)('0' if father isn't in dataset)
# Within-family ID of mother (Maternal ID)('0' if mother isn't in dataset)
# Sex code ('1' = male, '2' = female, '0' = unknown)
# Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
# 第7th 8th 是第一个SNP的alleles. 9th, 10th 是第二个SNP的alleles. 以此类推。。。 0 0代表missing data
read_ped <- function(ped_file,map_file=NULL,info_file=NULL){
  d1 <- read.delim(ped_file,sep=' ',stringsAsFactors = F,header=F,colClasses = 'character')
  if(ncol(d1)==1){
    d1 <- read.delim(ped_file,sep='\t',stringsAsFactors = F,header=F,colClasses = 'character')
  }
  ind <- d1$V2 ## individual id
  d2 <- d1[,7:ncol(d1)]
  x1 <- seq(1,ncol(d2)-1,by=2)
  x2 <- seq(2,ncol(d2),by=2)
  d21 <- d2[,x1]
  d22 <- d2[,x2]
  rownames(d21) <- sprintf("%s_allele1",ind)
  rownames(d22) <- sprintf("%s_allele2",ind)
  if(is.null(map_file)==F) f1 <- read.delim(map_file,stringsAsFactors = F,header=F)[,c(2,4)]
  if(is.null(info_file)==F) f1 <- read.delim(info_file,stringsAsFactors = F,header=F)
  #f1 <- read.delim(info_file,sep='\t',stringsAsFactors = F,header=F,colClasses = 'character')
  colnames(d21) <- colnames(d22) <- as.character(f1[,2])
  d2 <- rbind(d21,d22)
  return(list(allele1=d21,allele2=d22,all=d2,info=f1))
}
read_tped <- function(tped_file,tfam_file=NULL){
  d1 <- read.delim(tped_file,stringsAsFactors = F,header = F,sep=' ')
  s1 <- read.delim(tfam_file,stringsAsFactors = F,header = F,sep=' ')
  ind <- s1$V2
  f1 <- d1[,c(2,4)]
  d2 <- d1[,5:ncol(d1)]
  x1 <- seq(1,ncol(d2)-1,by=2)
  x2 <- seq(2,ncol(d2),by=2)
  d21 <- d2[,x1]
  d22 <- d2[,x2]
  d21 <- data.frame(t(d21),stringsAsFactors=F); 
  d22 <- data.frame(t(d22),stringsAsFactors=F)
  rownames(d21) <- sprintf("%s_allele1",ind)
  rownames(d22) <- sprintf("%s_allele2",ind)
  colnames(d21) <- colnames(d22) <- as.character(f1[,2])
  d2 <- rbind(d21,d22)
  return(list(allele1=d21,allele2=d22,all=d2,info=f1))
}
## read block --> tab SNP, involve pos
read_block <- function(block_file,map_file=NULL,info_file=NULL){
  b1 <- read.table(block_file,stringsAsFactors = F,header=T)
  if(is.null(map_file)==F) m1 <- read.delim(map_file,stringsAsFactors = F,header=F)[,c(2,4)]
  if(is.null(info_file)==F) m1 <- read.delim(info_file,stringsAsFactors = F,header=F)
  b1_list <- apply(b1,1,function(x){
    w1 <- which(m1[,2]>=x['BP1'] & m1[,2]<=x['BP2'])
    x1 <- m1[w1,2]
    r1 <- unlist(strsplit(x['SNPS'],'\\|'))
    x2 <- m1[which(m1[,1] %in% r1),2]
    list(all_pos=x1,mark_pos=x2)
  })
  names(b1_list) <- sprintf('Block%s',1:length(b1_list))
  return(b1_list)
}
## strategy: overlap/cover
get_share_block <- function(block_list,start_pos,end_pos,
                            strategy='overlap',
                            min_size=1000){
  # get block margin
  block_margin <- lapply(block_list,function(x){
    do.call(rbind,lapply(x,function(xx)c(min(xx$all_pos),max(xx$all_pos))))
  })
  names(block_margin) <- names(block_list)
  # filter by start end
  if(strategy == 'overlap'){
    block_use <- lapply(block_margin,function(x){
      x1 <- which(x[,2]>=start_pos)
      x2 <- which(x[,1]<=end_pos)
      x[intersect(x1,x2),]
    }) 
  }
  if(strategy == 'cover'){
    block_use <- lapply(block_margin,function(x){
      x1 <- which(x[,1]>=start_pos)
      x2 <- which(x[,2]<=end_pos)
      x[intersect(x1,x2),]
    }) 
  }
  # filter by size
  block_use <- lapply(block_use,function(x){
    x[which(x[,2]-x[,1]>=min_size),,drop=F]
  })
  names(block_use) <- names(block_list)
  print(unlist(lapply(block_use,nrow)))
  # get block info
  use_block <- lapply(names(block_list),function(x){
    block_list[[x]][rownames(block_use[[x]])]
  })
  names(use_block) <- names(block_list)
  # get shared block
  block_use2pos <- lapply(block_use,function(x){
    unique(apply(x,1,function(xx)xx[1]:xx[2]))
  })
  N = length(block_list)
  x <- table(unlist(block_use2pos))
  x1 <- as.numeric(names(which(x==N)))
  x11 <- x1[1:(length(x1)-1)]
  x12 <- x1[2:(length(x1))]
  x2 <- x12-x11
  w1 <- which(x2!=1)
  w11 <- c(1,w1+1) # start
  w12 <- c(w1,length(x1)) # end
  share_block <- cbind(x1[w11],x1[w12])
  share_block <- share_block[which(share_block[,2]-share_block[,1]>=min_size),,drop=F]
  share_block <- lapply(1:nrow(share_block),function(x)share_block[x,])
  names(share_block) <- sprintf('Block%s',1:length(share_block))
  block_use_share <- lapply(names(share_block),function(x){
    x1<-share_block[[x]];x2<-x1[1]:x1[2]
    x3 <- lapply(names(block_use),function(xx){
      xx1 <- rownames(block_use[[xx]])
      xx2 <- unique(unlist(lapply(block_list[[xx]][xx1],function(xxx)xxx$mark_pos)))
      intersect(xx2,x2)
    })
    x4 <- unlist(x3);x5<-as.numeric(names(which(table(x4)==length(x3))));
    list(all_pos=x5,mark_pos=x5)
  })
  names(block_use_share) <- names(share_block)
  return(list(share_block=share_block,block_use=block_use,block_use_share=block_use_share))
}
plot_share_block <- function(share_block,cex.block=1){
  tmp1 <- as.numeric(unlist(share_block))
  min_pos <- min(tmp1); max_pos <- max(tmp1); 
  N=length(share_block$block_use)
  par(mar=c(2,4,2,2))
  plot_new(xlim=c(min_pos,max_pos),ylim=c(0,N+1))
  cc <- get.class.color(names(share_block$block_use))
  for(i in 1:length(names(share_block$block_use))){
    for(j in rownames(share_block$block_use[[i]])){
      rect(xleft=share_block$block_use[[i]][j,1],
           xright=share_block$block_use[[i]][j,2],
           ybottom = i,ytop=i+0.5,border = cc[i],
           lwd=1,col=adjustcolor(cc[i],0.4));
    }
  }
  text(x=min_pos,y=c(0:3)+0.25,c('Shared',names(cc)),xpd=T,pos=2)
  for(i in names(share_block$share_block)){
    abline(v=share_block$share_block[[i]][1],lty=2,col='grey',lwd=0.5)
    abline(v=share_block$share_block[[i]][2],lty=2,col='grey',lwd=0.5)
    rect(xleft=share_block$share_block[[i]][1],
         xright=share_block$share_block[[i]][2],
         ybottom = 0,ytop=0.5,border = 'grey',lwd=1,col=adjustcolor('grey',0.4))
    nsnp <- share_block$block_use_share[[i]]$mark_pos
    text(x=share_block$share_block[[i]][1]/2+share_block$share_block[[i]][2]/2,
         y=0,sprintf('%s\n(SNP=%s)',i,length(nsnp)),pos=1,xpd=T,
         cex=cex.block)
    segments(x0=nsnp,x1=nsnp,
             y0=0.15,y1=0.35,lwd=0.25,cex=cex.block)
  }
  ss <- seq(min_pos,max_pos,length.out=11);
  rr <- floor(min(log10(ss)))-5
  ss1 <- round(ss/(10^rr))*10^rr
  axis(side=3,at=ss,labels=get_label_manual(ss1),xpd=TRUE,cex.axis=0.75)
}
##
add_miss <- function(x,p){
  x1 <- data.frame(matrix(0,nrow=nrow(x),ncol=length(p)))
  rownames(x1) <- rownames(x)
  colnames(x1) <- p
  for(i in 1:ncol(x)) x1[,colnames(x)[i]] <- x[,i]
  x1
}
ped2eset <- function(ped_dat){
  m1 <- t(as.matrix(ped_dat$all))
  rownames(m1) <- rownames(ped_dat$info)
  f1 <- ped_dat$info
  info <- data.frame(mark=as.character(f1$V4),
                     pos=as.numeric(f1$V4),
                     end_pos=as.numeric(f1$V4),
                     rs_mark=f1$V2,stringsAsFactors = F)
  ## 
  w1 <- which(table(info$mark)>=2)
  if(length(w1)==0){
    rownames(info) <- info$mark
    rownames(m1) <- rownames(info)
  }else{
    w11 <- which(info$mark %in% names(w1))
    info1 <- info[w11,]
    info2 <- info[setdiff(1:nrow(info),w11),]
    tmp1 <- lapply(names(w1),function(x){
      x1 <- which(info$mark==x)
      m11 <- m1[x1,]
      apply(m11,2,function(x)paste(unique(x),collapse=''))
    })
    tmp1 <- do.call(rbind,tmp1)
    m2 <- rbind(m1[setdiff(1:nrow(info),w11),],tmp1)
    tmp2 <- lapply(names(w1),function(x){
      x1 <- which(info$mark==x)
      m11 <- info[x1,]
      apply(m11,2,function(x)paste(unique(x),collapse=''))
    })
    tmp2 <- do.call(rbind,tmp2)
    info2 <- rbind(info[setdiff(1:nrow(info),w11),],tmp2)
    od1 <- order(as.numeric(info2$pos))
    info <- info2[od1,]
    m1 <- m2[od1,]
    rownames(info) <- info$mark
    rownames(m1) <- rownames(info)
  }
  eset <- generate.eset(exp_mat=m1, ## pos * sample
                        feature_info = info)
  return(eset) 
}
merge_ped <- function(ped_list){
  if(class(ped_list[[1]])=='data.frame') ped_list <- list(each=ped_list)
  pn <- names(ped_list)
  all_pos <- do.call(rbind,lapply(ped_list,function(x)x$info))
  tmp1 <- aggregate(all_pos[,1],list(all_pos[,2]),function(x)paste(x,collapse='|'))
  p2 <- tmp1$Group.1; p1 <- tmp1$x
  info <- data.frame(mark=as.character(p2),pos=as.numeric(p2),end_pos=as.numeric(p2),rs_mark=p1,stringsAsFactors = F)
  rownames(info) <- info$mark
  new_ped_list <- lapply(ped_list,function(x){
    list(allele1=add_miss(x$allele1,p2),allele2=add_miss(x$allele2,p2),all=add_miss(x$all,p2),info=info)
  })
  new_ped_list_combine <- do.call(rbind,lapply(new_ped_list,function(x)x$all))
  new_ped_list_sample_info <- do.call(rbind,lapply(pn,function(x)cbind(rownames(new_ped_list[[x]]$all),x)))
  rownames(new_ped_list_combine) <- new_ped_list_sample_info[,1]
  new_ped_list_sample_info <- as.data.frame(new_ped_list_sample_info,stringsAsFactors=FALSE)
  colnames(new_ped_list_sample_info) <- c('sample_name','group')
  new_ped_list_sample_info$ori_sample_name <- gsub('(.*)_allele.*','\\1',new_ped_list_sample_info$sample_name)
  new_ped_list_pos_info <- info
  rownames(new_ped_list_pos_info) <- new_ped_list_pos_info$mark ## position info --> like feature
  rownames(new_ped_list_sample_info) <- new_ped_list_sample_info$sample_name ## sample info --> like pheno
  eset <- generate.eset(exp_mat=t(new_ped_list_combine), ## pos * sample
                        phenotype_info = new_ped_list_sample_info[rownames(new_ped_list_combine),],
                        feature_info = new_ped_list_pos_info[colnames(new_ped_list_combine),])
  return(eset) 
}
change_sampleName <- function(eset,transfer_tab,from_col,to_col){
  p1 <- pData(eset)
  p1$old_ori_sample_name <- p1$ori_sample_name
  p1$old_sample_name <- p1$sample_name
  rownames(transfer_tab) <- transfer_tab[,from_col]
  p1$ori_sample_name <- transfer_tab[p1$old_ori_sample_name,to_col]
  tmp1 <- gsub('(.*)_(allele.*)','\\2',p1$old_sample_name)
  p1$sample_name <- sprintf('%s_%s',p1$ori_sample_name,tmp1)
  rownames(p1) <- p1$sample_name
  m1 <- exprs(eset); 
  colnames(m1) <- p1$sample_name
  eset <- generate.eset(exp_mat=m1, ## pos * sample
                        phenotype_info = p1,
                        feature_info =fData(eset))
  return(eset) 
}
##
combine_pos_ped <- function(x){
  mat <- exprs(x)
  pos <- fData(x)
  w1 <- apply(mat,1,function(x)length(unique(x)))
  w2 <- which(w1==1)
  cc  <- list();mark <- pos$mark[1]
  if(length(w2)==0) return(x) ## no need to combine
  for(i in 2:length(w2)){
    x0 = w2[i-1]
    x1 = w2[i]
    if(x1==x0+1){
      cc[[mark]] <- c(cc[[mark]],x1)
    }else{
      mark <- pos$mark[i]
      cc[[mark]] <- c(x1)
    }
  }
  cc_use <- cc[which(unlist(lapply(cc,length))>1)] ## combine
  if(length(cc_use)==0) return(x) ## no need to combine
  w3 <- setdiff(1:nrow(pos),unlist(cc_use)) ## keep
  pos_keep <- pos[w3,]
  mat_keep <- mat[w3,]
  mat_new <- list(); pos_new <- list()
  for(i in 1:length(cc_use)){
    x1 <- cc_use[[i]]
    x2 <- pos[x1,]
    pos1 <- min(x2$pos); pos2 <- max(x2$pos);
    mark <- paste(x2$rs_mark,collapse='|')
    nt <- paste(apply(mat[x1,],1,function(x)unique(x)),collapse='')
    pos_new[[i]] <- c(paste(pos1,pos2,sep='_'),pos1,pos2,mark)
    mat_new[[i]] <- rep(nt,length.out=ncol(mat))
  }
  pos_new <- do.call(rbind,pos_new)
  mat_new <- do.call(rbind,mat_new)
  pos_new <- data.frame(pos_new,stringsAsFactors=F);
  colnames(pos_new) <- colnames(pos)
  pos_new$pos <- as.numeric(pos_new$pos)
  pos_new$end_pos <- as.numeric(pos_new$end_pos)
  rownames(mat_new) <- pos_new$mark
  colnames(mat_new) <- colnames(mat)
  mat_combine <- rbind(mat_keep,mat_new)
  pos_combine <- rbind(pos_keep,pos_new)
  pos_combine <- pos_combine[order(pos_combine$pos),];
  rownames(pos_combine) <- pos_combine$mark
  mat_combine <- mat_combine[rownames(pos_combine),]
  eset <- generate.eset(exp_mat=mat_combine, ## pos * sample
                        phenotype_info = pData(x),
                        feature_info = pos_combine)
  return(eset) 
}
##
random_select_eset <- function(eset,n,use_sample){
  mat <- exprs(eset)
  phe <- pData(eset)
  if(is.null(use_sample)==FALSE){
    w0 <- which(colnames(mat) %in% use_sample)
    w1 <- sample(setdiff(1:ncol(mat),w0),n-length(w0))
    w1 <- sort(c(w0,w1))
  }else{
    w1 <- sample(1:ncol(mat),n)  
  }
  eset <- generate.eset(exp_mat=mat[,w1], ## pos * sample
                        phenotype_info = pData(eset)[w1,],
                        feature_info = fData(eset))
  return(eset) 
}
## filter positions with too missing value
filter_pos <- function(eset,thre=1,use_pos=NULL,remove_identical=FALSE){
  mat <- exprs(eset)
  phe <- pData(eset)
  pos <- fData(eset)
  w1 <- apply(mat,1,function(x)length(which(x=='0')))
  w2 <- w1/ncol(mat)
  w3 <- which(w2<=thre)
  pos <- pos[w3,]
  if(is.null(use_pos)==FALSE){
    pos <- pos[intersect(as.character(use_pos),rownames(pos)),]
  }
  if(remove_identical==TRUE){
    x <- mat[rownames(pos),,drop=F]
    w1 <- apply(x,1,function(x)length(unique(x)))
    w2 <- which(w1>1)
    x <- x[w2,,drop=F]
    pos <- pos[rownames(x),,drop=F]
  }
  eset <- generate.eset(exp_mat=mat[rownames(pos),,drop=T], ## pos * sample
                        phenotype_info = pData(eset),
                        feature_info = pos)
  return(eset) 
}
filter_sample <- function(eset,thre=1,use_sample=NULL){
  mat <- exprs(eset)
  phe <- pData(eset)
  pos <- fData(eset)
  w1 <- apply(mat,2,function(x)length(which(x=='0')))
  w2 <- w1/nrow(mat)
  w3 <- which(w2<=thre)
  phe <- phe[w3,]
  if(is.null(use_sample)==FALSE){
    phe <- phe[which(phe$ori_sample_name %in% use_sample),]
  }
  eset <- generate.eset(exp_mat=mat[,rownames(phe),drop=T], ## pos * sample
                        phenotype_info = phe,
                        feature_info = pos)
  return(eset) 
}
filter_allele <- function(eset,thre=1,use_allele=NULL){
  mat <- exprs(eset)
  phe <- pData(eset)
  pos <- fData(eset)
  w1 <- apply(mat,2,function(x)length(which(x=='0')))
  w2 <- w1/nrow(mat)
  w3 <- which(w2<=thre)
  phe <- phe[w3,,drop=F]
  if(is.null(use_allele)==FALSE){
    phe <- phe[which(rownames(phe) %in% use_allele),,drop=F]
  }
  eset <- generate.eset(exp_mat=mat[,rownames(phe),drop=T], ## pos * sample
                        phenotype_info = phe,
                        feature_info = pos)
  return(eset) 
}

find_block_eset <- function(eset,target_pos,inner_text=FALSE){
  x <- exprs(eset); f1 <- fData(eset)
  w1 <- apply(x,1,function(x)length(unique(x)))
  w2 <- which(w1>1) ## not identical
  ww <- which(fData(eset)[,2]==target_pos)
  if(length(ww)==0) return(NULL)
  wwd <- w2-ww
  wwd1 <- max(which(wwd<0))
  wwd2 <- min(which(wwd>0))
  u1 <- c(w2[wwd1],w2[wwd1]+1,ww,ww+1,w2[wwd2]) # start
  u2 <- c(w2[wwd1],ww-1,ww,w2[wwd2]-1,w2[wwd2]) # end
  s1 <- paste0(x[(w2[wwd1]+1):(ww-1),1],collapse = '')
  s2 <- paste0(x[(ww+1):(w2[wwd2]-1),1],collapse = '')
  if(inner_text==FALSE){s1='0';s2='0';}
  new_mat <- rbind(x[w2[wwd1],],s1,x[ww,],s2,x[w2[wwd2],])
  rownames(new_mat) <- c(f1[w2[wwd1],1],
                         sprintf('%s-%s',f1[w2[wwd1]+1,1],f1[ww-1,1]),
                         f1[ww,1],
                         sprintf('%s-%s',f1[ww+1,1],f1[w2[wwd2]-1,1]),
                         f1[w2[wwd2],1])
  eset <- generate.eset(new_mat)
}
alter_rowname <- function(eset,new_col=1){
  m1 <- exprs(eset); f1 <- fData(eset);
  rownames(m1) <- f1[,new_col]; rownames(f1) <- f1[,new_col]
  eset <- generate.eset(m1,feature_info = f1)
}
##
col_ATCG <- c(brewer.pal(8,'Accent')[1],brewer.pal(11,'Set3')[4],brewer.pal(11,'Set3')[2],
              brewer.pal(11,'Paired')[1],'grey')
names(col_ATCG) <- c('A','T','C','G','0') ## # A-green, T-red, C-yellow, G-blue 
##
draw_pos2sample_pdf <- function(pdf_file,eset,...){
  ww <- dim(exprs(eset))[1]/15+6;hh <- dim(exprs(eset))[2]/10+2;
  pdf(pdf_file,height=hh,width=ww)
  draw_pos2sample(eset,...)
  dev.off()
}
draw_pos2sample <- function(eset,point_strategy='color',order_strategy='data',
                            remove_identical=FALSE,show_sample_name=TRUE,
                            show_pos=TRUE,
                            mark_pos=NULL,mark_block=NULL,
                            top_hap_number=0,top_hap_number_thre=0.05,
                            only_mark_tag_SNP=TRUE,
                            intersect_thre=0.1,cex.point=0.6){
  x <- exprs(eset)
  if(remove_identical==TRUE){
    w1 <- apply(x,1,function(x)length(unique(x)))
    w2 <- which(w1>1)
    x <- x[w2,]
  }
  #p1 <- 3; p2 <- 3; p3 <- 4; p4 <- 3;
  #if(show_sample_name==TRUE) p2 <- 8;
  #if(top_hap_number>0 & is.null(mark_block)==FALSE) p1 <- 3+top_hap_number/3
  #par(mar=c(p1,p2,p3,p4))
  x1 <- ifelse(x=='A',1,ifelse(x=='T',2,ifelse(x=='C',3,ifelse(x=='G',4,ifelse(x=='0',5,10)))))
  if(order_strategy=='data'){o1 <- hclust(as.dist(1-cor(x1)))$order;x1 <- x1[,o1]}
  if(order_strategy=='name'){x1 <- x1[,sort(colnames(x1))]}
  x <- x[rownames(x1),colnames(x1)]
  plot(1,xlim=c(1,nrow(x1)),ylim=c(1,ncol(x1)),xaxt='n',yaxt='n',col='white',bty='n',xlab='',ylab='')
  pp <- par()$usr
  ##
  if(is.null(mark_block)==FALSE){
    p1 <- fData(eset)[rownames(x1),]
    for(ni in names(mark_block)){
      i <- mark_block[[ni]]
      all_pos_min <- min(i$all_pos)
      all_pos_max <- max(i$all_pos)
      w1 <- which(p1$pos>=all_pos_min)
      w2 <- which(p1$end_pos<=all_pos_max)
      w3 <- intersect(w1,w2)
      if(length(w3)==0) next
      rect(xleft=min(w3)-0.25,xright=max(w3)+0.25,ybottom =-2,ytop=-1,col='light grey',border=NA,xpd=TRUE)
      w4 <- which(p1$pos %in% i$mark_pos)
      segments(x0=w4,x1=w4,y0=0.5,y1=0.5+ncol(x1),lwd=0.5,lty=3,col='grey')
      points(x=w4,y=rep(-1,length.out=length(w4)),pch=17,cex=0.4,xpd=T,col='red',xpd=T)
      text(mean(w3),-1.5,ni,cex=0.4,xpd=T)
    }
    if(top_hap_number>0){
      h2p <- list();s2p <- list()
      for(ni in names(mark_block)){
        i <- mark_block[[ni]]
        all_pos_min <- min(i$all_pos)
        all_pos_max <- max(i$all_pos)
        w1 <- which(p1$pos>=all_pos_min)
        w2 <- which(p1$end_pos<=all_pos_max)
        w3 <- intersect(w1,w2)
        if(length(w3)==0) next
        x31 <- min(w3)
        x32 <- max(w3)
        cont1 <- x[x31:x32,,drop=F]
        w4 <- which(p1$pos %in% i$mark_pos)
        if(only_mark_tag_SNP==TRUE) cont1 <- x[w4,,drop=F]
        cont2 <- apply(cont1,2,function(x)paste(x,collapse='_'))
        t1 <- sort(table(cont2),decreasing = T)
        t1 <- t1[1:min(top_hap_number,length(t1))]
        t1p <- t1/ncol(x1); 
        w1 <- which(t1p>=top_hap_number_thre) ## 
        t1 <- t1[w1]
        t2 <- names(t1)
        c1 <- colorRampPalette(brewer.pal(8,'Spectral'))(length(t1))
        c1 <- adjustcolor(c1,0.5)
        names(c1) <- names(t1)
        for(j in names(t1)){
          w4 <- which(cont2==j)
          segments(x0=x31,x1=x32,y0=w4,y1=w4,col=c1[j],lwd=2)
        }
        ss <- seq(-2,by=-1,length.out=length(t1));names(ss) <- names(t1);
        h2p[[ni]] <- data.frame(y=ss-1.5,x1=x31,x2=x32,row.names = names(t1))
        s2p[[ni]] <- cont2
        rect(xleft=min(w3)-0.25,xright=max(w3)+0.25,ybottom =ss-1.7,ytop=ss-1.2,col=c1,border=NA,xpd=TRUE,xpd=T)
        text(x=x31/2+x32/2,y=ss-1.5,sprintf("%s (n=%s,%s)",gsub("_","",names(t1)),
                                            t1,signif(t1/ncol(x1),2)),xpd=T,cex=0.3)
      }
      ## interaction
      if(length(h2p)>=2){
        for(i2 in 2:length(h2p)){
          i1 <- i2-1;
          b1 <- h2p[[i1]];b2 <- h2p[[i2]];
          for(each_b1 in rownames(b1)){
            for(each_b2 in rownames(b2)){
              #count1 <- length(which(s2p[[i1]]==each_b1 & s2p[[i2]]==each_b2))/ncol(x1)
              count1 <- length(which(s2p[[i1]]==each_b1 & s2p[[i2]]==each_b2))/length(which(s2p[[i1]]==each_b1))
              if(count1>intersect_thre){
                segments(x0=b1[each_b1,'x2'],x1=b2[each_b2,'x1'],
                         y0=b1[each_b1,'y'],y1=b2[each_b2,'y'],lwd=count1*2,
                         xpd=T,col=adjustcolor('black',0.4))
                ll <- length(rownames(b1))
                text(b1[each_b1,'x2'],(2*ll-1)*b1[each_b1,'y']/(2*ll)+b2[each_b2,'y']/(2*ll),
                     signif(count1,2),cex=0.15,adj=1,xpd=T)
              }
            }
          }
        }
      }
    }
  }
  if(is.null(mark_pos)==FALSE){
    p1 <- fData(eset)[rownames(x1),]
    for(i in mark_pos){
      m1 <- i
      w1 <- which(m1>=p1$pos)
      w2 <- which(m1<=p1$end_pos)
      w3 <- intersect(w1,w2)
      if(length(w3)==0) w3 <- max(w1)+0.5 else w3 <- min(w3)
      segments(x0=w3,x1=w3,y0=0.5,y1=0.5+ncol(x1),lwd=0.5)
      text(x=w3,y=0,i,cex=0.3,xpd=T)
    }
  }
  for(i in 1:nrow(x1)){ ## for each position
    v1 <- unique(x1[i,])
    if(length(v1)==1){c1 = 'dark grey';p1 <- 1;}else{c1 <- col_ATCG[x[i,]];c1[which(is.na(c1)==TRUE)]<-'purple';p1<-16}
    if(point_strategy=='color')points(rep(i,length.out=ncol(x1)),1:ncol(x1),pch=p1,col=c1,xpd=T,cex=cex.point) ## 
    if(point_strategy=='text')text(rep(i,length.out=ncol(x1)),1:ncol(x1),x[i,],cex=0.3,font=2,xpd=T) ## 
  }
  if(show_sample_name==TRUE) text(0,1:ncol(x1),colnames(x1),adj=1,xpd=T,cex=0.5)
  if(show_pos==TRUE) text(1:nrow(x1),ncol(x1)+1,rownames(x1),adj=0,srt=60,cex = 0.5,xpd=T)
  if(point_strategy=='color'){legend(1+nrow(x1),pp[3]/2+pp[4]/2,
                                     c('A','T','C','G','missing','complex','identical'),
                                     col=c(col_ATCG,'purple','dark grey'),
         pch=c(16,16,16,16,16,16,1),xpd=T,bty='n',yjust = 0.5,cex=0.6,xjust=0)
  }
}
##
get_block_hap <- function(eset,mark_block,remove_miss=TRUE,only_mark_tag_SNP=FALSE){
  x <- exprs(eset)
  if(is.null(mark_block)==FALSE){
    p1 <- fData(eset)
    s2p <- list()
    for(ni in names(mark_block)){
       i <- mark_block[[ni]]
       all_pos_min <- min(i$all_pos)
       all_pos_max <- max(i$all_pos)
       w1 <- which(p1$pos>=all_pos_min)
       w2 <- which(p1$end_pos<=all_pos_max)
       w3 <- intersect(w1,w2)
       if(length(w3)==0) next
       x31 <- min(w3)
       x32 <- max(w3)
       cont1 <- x[x31:x32,,drop=F]
       w4 <- which(p1$pos %in% i$mark_pos)
       if(only_mark_tag_SNP==TRUE) cont1 <- x[w4,,drop=F]
       cont2 <- apply(cont1,2,function(x)paste(x,collapse='_'))
       #t1 <- sort(table(cont2),decreasing = T)
       #t2 <- names(t1)
       if(remove_miss==TRUE) cont2[grep('0',cont2)] <- NA 
       s2p[[ni]] <- cont2
      }
  }
  x1 <- do.call(rbind,s2p)
  eset <- generate.eset(exp_mat=x1, ## pos * sample
                        phenotype_info = pData(eset),
                        feature_info = NULL)
  return(eset) 
}
test_hap_assoc <- function(eset_hap,group='group',pv_thre=0.01){
  m1 <- exprs(eset_hap)
  s1 <- pData(eset_hap)[,group]
  names(s1) <- rownames(pData(eset_hap))
  pv <- apply(m1,1,function(x){
    t2 <- table(list(x,s1))
    print(colSums(t2))
    chisq.test(t2)$p.value
  })
  w1 <- which(pv<pv_thre)
  if(length(w1)==0) return(list(pv=pv))
  sig_block <- list()
  for(i in w1){
    x=m1[i,]
    t2 <- table(list(x,s1))
    pv1 <- apply(t2,1,function(x){
      chisq.test(cbind(x,colSums(t2)))$p.value
    })
    pv1 <- p.adjust(pv1)
    w2 <- which(pv1<pv_thre)
    t3 <- do.call(rbind,lapply(w2,function(w3)c(t2[w3,]/colSums(t2),colSums(t2))))
    sig_block[[rownames(m1)[i]]] <- cbind(t2[w2,,drop=F],pv=pv1[w2],t3)
  }
  return(list(pv=pv,sig_block=sig_block))
}
get_mark <- function(use_pos){
  x1 <- hg19_pos[which(hg19_pos$V2==use_pos),'V4']
  x2 <- gsub('(.*):(.*)','\\2',x1)
  x2
}
test_similarity <- function(eset){
  m1 <- exprs(eset)
  all1 <- colnames(m1)
  res <- list()
  for(i in all1){
    for(j in all1){
      if(i!=j){
        v1 <- length(which(m1[,i]==m1[,j]))
        res[[sprintf('%s_%s',i,j)]] <- c(i,j,v1)
      }
    }
  }
  res <- do.call(rbind,res)
  res <- as.data.frame(res,stringsAsFactors=F);colnames(res) <- c('SA1','SA2','Overlap')
  res$Overlap <- as.numeric(res$Overlap)
  res$Percentage <- res$Overlap/nrow(m1)
  res <- res[order(res$Percentage,decreasing = T),]
  return(res)
}
test_similarity_pair <- function(eset){
  m1 <- exprs(eset)
  all1 <- colnames(m1)
  res <- matrix(nrow(m1),nrow=length(all1),ncol=length(all1))
  rownames(res) <- colnames(res) <- all1
  for(i in all1){
    for(j in all1){
      if(i!=j){
        v1 <- length(which(m1[,i]==m1[,j]))
        res[i,j] <- v1
      }
    }
  }
  return(res)
}

##
merge_eset_multi <- function(eset_list){
  x1<-lapply(eset_list,exprs)
  x2<-lapply(eset_list,pData)
  x2<-lapply(names(x2),function(x)data.frame(group=x,x2[[x]]))
  x1_tmp <- unique(unlist(lapply(x1,rownames)))
  x11 <- do.call(cbind,lapply(x1,function(x)as.data.frame(x)[x1_tmp,]))
  x22 <- do.call(rbind,x2);rownames(x22) <- colnames(x11)
  generate.eset(x11,phenotype_info = x22)
}
##
get_allele_from <- function(child_allele,parent_allele,main='',per_thre=0.95,n2pos=NULL,
                            child_allele_AF=NULL,AFdiff_thre=0.1){
  # 0,0.5,1
  if(is.null(child_allele_AF)==FALSE){
    child_allele_AF_mod <- sapply(child_allele_AF,function(x)min(abs(x-0.5),x,1-x))
  }
  ss <- round(seq(1,length(n2pos),length.out=10))
  r1 <- lapply(1:length(child_allele),function(x){
    colnames(parent_allele)[which(parent_allele[x,]==child_allele[x])]
  })  
  r2 <- table(unlist(r1))/length(r1)
  max_r <- names(which.max(r2))
  r3 <- lapply(r1,function(x){
    if(max_r %in% x)return(max_r) else return(x) ## if max_r exist, remove others
  })
  r4 <- table(unlist(r3))/length(r3)
  r5 <- do.call(cbind,lapply(names(r4),function(x){
    unlist(lapply(r3,function(xx)ifelse(x %in% xx,1,0)))
  }))
  colnames(r5) <- names(r4)
  ##
  par(mar=c(4,4,8,4))
  r6 <- draw_image(r5,cc=c('grey',2),add_line = F)
  text(0,r6$yy1,colnames(r5),xpd=T,pos=2)
  text(x=0.5,y=2,sprintf('%s from %s\n%s%s',main,max_r,signif(r2[max_r]*100,2),'%'),xpd=T)
  text(x=r6$xx1[ss],y=-0.25,n2pos[ss],xpd=T,srt=90,adj=1,cex=0.5)
  if(is.null(child_allele_AF)==FALSE){
    text(x=r6$xx1,y=r6$pp[4]+0.01+child_allele_AF/5,'*',xpd=T,
         cex=ifelse(child_allele_AF_mod<AFdiff_thre,0.3,0.5),
         col=ifelse(child_allele_AF_mod<AFdiff_thre,'light grey','black'))
  }
  ## if only one mix
  if(r4[max_r]<per_thre){
    for(i in setdiff(colnames(r5),max_r)){
    t1 <- sum(r5[,max_r]+r5[,i])/nrow(r5)
    if(t1>=per_thre){
      r7 <- lapply(r3,function(x){
        if(max_r %in% x){return(max_r)}else{
          if(i %in% x) return(i) else return(x)
        }
      })
      r5 <- do.call(cbind,lapply(names(r4),function(x){
        unlist(lapply(r7,function(xx)ifelse(x %in% xx,1,0)))
      }))
      colnames(r5) <- names(r4)
      r6 <- draw_image(r5,cc=c('grey',2),add_line = F)
      text(0,r6$yy1,colnames(r5),xpd=T,pos=2)
      text(x=0.5,y=1.6,sprintf('%s from %s>%s \n%s%s',
                               main,max_r,i,signif(t1*100,2),'%'),xpd=T)
      text(x=r6$xx1[ss],y=-0.25,n2pos[ss],xpd=T,srt=90,adj=1,cex=0.5)
      if(is.null(child_allele_AF)==FALSE){
        text(x=r6$xx1,y=r6$pp[4]+0.01+child_allele_AF/5,'*',xpd=T,
             cex=ifelse(child_allele_AF_mod<AFdiff_thre,0.3,0.5),
             col=ifelse(child_allele_AF_mod<AFdiff_thre,'light grey','black'))
      }
    }
  }
  }
  ##
}
##
get_allele_from_inputOrder <- function(child_allele,parent_allele,main='',per_thre=0.95,n2pos=NULL,allele_order=NULL){
  ss <- round(seq(1,length(n2pos),length.out=10))
  r1 <- lapply(1:length(child_allele),function(x){
    colnames(parent_allele)[which(parent_allele[x,]==child_allele[x])]
  })
  r2 <- table(unlist(r1))/length(r1)
  for(i in 1:length(allele_order)){
    r3 <- lapply(r1,function(x){
      if(allele_order[i] %in% x) return(allele_order[i]) else return(x) ## if max_r exist, remove others
    })
    r4 <- table(unlist(r3))/length(r3)
    r5 <- do.call(cbind,lapply(names(r4),function(x){
      unlist(lapply(r3,function(xx)ifelse(x %in% xx,1,0)))
    }))
    colnames(r5) <- names(r4)
    ##
    r6 <- draw_image(r5,cc=c('grey',2),add_line = F)
    text(0,r6$yy1,colnames(r5),xpd=T,pos=2)
    ii <- paste(allele_order[1:i],collapse='>')
    text(x=0.5,y=1.4,sprintf('%s from %s\n%s%s',main,ii,signif(sum(r2[allele_order[1:i]])*100,3),'%'),xpd=T)
    text(x=r6$xx1[ss],y=-0.25,n2pos[ss],xpd=T,srt=90,adj=1,cex=0.5)
    r1 <- r3
    r2 <- table(unlist(r1))/length(r1);
  }
}
##
FM_percentage <- function(Father_GT,Mother_GT,Child_AF,sep_char='\\|'){
  split1 <- function(x,sep_char){
    x1 <- as.data.frame(do.call(rbind,strsplit(x,sep_char)))
    for(i in 1:ncol(x1)) x1[,i] <- as.numeric(x1[,i])
    return(x1)
  }
  c_af <- Child_AF
  f_gt <- split1(Father_GT,sep_char)
  m_gt <- split1(Mother_GT,sep_char)
  w_f <- which(f_gt[,1]!=f_gt[,2]) ## L/R
  w_m <- which(m_gt[,1]!=m_gt[,2]) ## L/R
  f1_int <- f_gt[,1]; f1_int[w_f] <- 1 ## first choose alt
  m1_int <- m_gt[,1]; m1_int[w_m] <- 1 ## first choose alt
  # Y=pF+(1-P)M --> Y-M=p(F-M)
  fm_diff <- f1_int-m1_int
  fit1 <- lm(c_af-m1_int~0+fm_diff)
  p <- coef(fit1);pdiff <- sum(abs(c_af-m1_int-predict(fit1)));
  pdiff_old <- pdiff;print(p);print(pdiff_old)
  #
  need_re <- union(w_f,w_m)
  for(r in 1:10){
    for(i in need_re){
      yLL=p*f_gt[i,1]+(1-p)*m_gt[i,1]
      yLR=p*f_gt[i,1]+(1-p)*m_gt[i,2]
      yRL=p*f_gt[i,2]+(1-p)*m_gt[i,1]
      yRR=p*f_gt[i,2]+(1-p)*m_gt[i,2]
      ydiff <- abs(c(yLL,yLR,yRL,yRR)-c_af[i])
      wy <- which.min(ydiff)
      if(wy==1){ f1_int[i] <- f_gt[i,1]; m1_int[i] <- m_gt[i,1];}
      if(wy==2){ f1_int[i] <- f_gt[i,1]; m1_int[i] <- m_gt[i,2];}
      if(wy==3){ f1_int[i] <- f_gt[i,2]; m1_int[i] <- m_gt[i,1];}
      if(wy==4){ f1_int[i] <- f_gt[i,2]; m1_int[i] <- m_gt[i,2];}
    }
    fm_diff <- f1_int-m1_int
    fit1 <- lm(c_af-m1_int~0+fm_diff)
    p <- coef(fit1);pdiff <- sum(abs(c_af-m1_int-predict(fit1)))
    if(abs(pdiff-pdiff_old)<0.1) break
    pdiff_old <- pdiff
    print(p);print(pdiff)
  }
  # output m1_int/f1_int, p
  print(p)
  print(str(m1_int))
  print(str(f1_int))
  return(list(F1_percentage=p,F1_GT=f1_int,M1_GT=m1_int))
}

