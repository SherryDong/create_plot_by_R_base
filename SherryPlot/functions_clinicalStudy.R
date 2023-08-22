library(NetBID2)
library(readxl)
library(gridExtra)
library(coxme)
library(ggsci)
library("scales")
library(merTools)
## match
library(table1)
library(ccoptimalmatch)
#library(Matching)
library(MatchIt)
cc <- pal_jama("default")(7)
source('D:/analysis_eng/SherryPlot/function_basic.R')
##
check_for_dupid <- function(x,id_col){
  w1 <- names(which(table(x[,id_col])>1));
  print(sprintf('Dup ID(Name+ID):%s',w1));
  x1 <- x[which(x[,id_col] %in% w1),]
  return(x1)
}
check_diff_for_dupid <- function(x,id_col,diff_col){
  tmp1 <- unique(x[,c(id_col,diff_col)])
  w1 <- names(which(table(tmp1[,id_col])>1))
  print(sprintf('Diff %s for same %s:%s',diff_col,id_col,w1));
  x1 <- x[which(x[,id_col] %in% w1),]
  return(x1)
}
NA2char <- function(x){
  x[which(is.na(x)==TRUE)] <- NA;
  x[which(x=='NA')] <- NA;
  x
}
small2NA <- function(x,min=5){
  x1 <- table(x);x2<-names(x1)[which(x1<min)]
  x[which(x%in%x2)]<-NA
  x
}
ToYesNo <- function(x){
  ifelse(x==1,'YES','NO')
}
df2list <- function(x,id_col='',use_id){
  x1 <- lapply(sort(unique(x[,id_col])),function(xx){
    xx1 <- unique(x[which(x[,id_col]==xx),use_id]);xx1[which(is.na(xx1)==F)]
  })
  names(x1) <- sort(unique(x[,id_col]));
  x11 <- unlist(lapply(x1,length));x1[which(x11>0)]
}
df2vec <- function(x,id_col='',use_id){
  x1 <- unique(x[,c(id_col,use_id)]);
  x1 <- x1[which(is.na(x1[,1])==F & is.na(x1[,2])==F),]
  x2 <- x1[,id_col];names(x2)<-x1[,use_id];x2
}
detect_NA <- function(dat,check_col){
  x1 <- apply(dat[,check_col,drop=F],2,function(x){
    which(is.na(x)==T)
  })
  setdiff(1:nrow(dat),unique(unlist(x1)))
}
removeNA_level <- function(x){
  if(class(x)!='factor') x<-factor(x)
  x1 <- as.character(x)
  x2 <- factor(x1,levels=intersect(levels(x),unique(x1)))
  x2
}
findDup <- function(dat,use_col){
  w1 <- table(dat[,use_col]); w2 <- names(w1)[which(w1>1)]; w2
}
require(survival);library(survminer)
#survfit,survdiff,coxph,cox.zph,lss
plot_Surv<-function(dat,time_col='TestNegative-Onset(day)',
                    group_col='group',
                    censor_col='Censor',main=''){
  dat1 <- data.frame(time=dat[,time_col],status=dat[,censor_col],group=dat[,group_col])
  fit1 <- survfit(Surv(time,status)~group,data=dat1)
  logrank.test <- survdiff(Surv(time,status)~group,data=dat1);
  plot1=ggsurvplot(fit1,data=dat1,risk.table = TRUE,pval=TRUE,conf.int=TRUE,palette='jama')
  plot1$plot$labels$y <- 'Test Positive probability';
  plot1$plot$labels$x <- 'Time to Onset (Day)';
  plot1$plot$labels$title <- sprintf('%s, N=%s',main,nrow(dat1))
  pv <- 1-pchisq(logrank.test$chisq,df=length(logrank.test$n)-1)
  return(list(plot=plot1,test=pv))
}
plot_Box <- function(dat,group_col='group',ylab='TestNegative-Onset(day)',val_col='TestNegative-Onset(day)',xlab='Group',main=''){
  dat[,group_col] <- factor(as.character(dat[,group_col]),levels=intersect(levels(dat[,group_col]),unique(dat[,group_col])))
  fit1 <- aov(dat[,val_col]~dat[,group_col])
  pv <- summary(fit1)[[1]][1,5]
  boxplot(dat[,val_col]~dat[,group_col],xlab=xlab,ylab=ylab,main=sprintf('%s\n(Pv=%s)',main,signif(pv,3)),cex.lab=1.2)
  return(pv)
}

plot_BarPre <- function(dat,group_col=NULL,use_col=c('Fever','Cough','Diarrhea','Vomiting','Fatigue','Sore throat','Body aches','Running nose','Febrile seizures')){
  if(is.null(group_col)==FALSE){
    t1 <- do.call(rbind,lapply(use_col,function(x)table(dat[,x],dat[,group_col])['YES',]));
    rownames(t1) <- use_col
    t2 <- t(t1)/as.numeric(table(dat[,group_col]))
  }else{
    t1 <- do.call(rbind,lapply(use_col,function(x)table(dat[,x])['YES']));rownames(t1) <- use_col;t2<-t(t1)
  }
  t2 <- t2[which(is.nan(rowSums(t2))==FALSE),,drop=F]
  return(t2)
}
plot_BarPre_simple <- function(dat,group_col=NULL,use_col=c('Fever')){
  t1 <- table(dat[,use_col],dat[,group_col]);
  t2 <- t(t1)/as.numeric(table(dat[,group_col]))
  return(t2)
}
plot_Bar <- function(dat,cc=pal_jama("default")(ncol(dat)),mar=c(8,5,6,3),ylab='Percentage',legend_pos='topright',xtext_srt=60,xtext_adj=1,main='',bty='n',...){
  if(ncol(dat)>7){cc <- colorRampPalette(pal_jama("default")(7))(ncol(dat))}
  names(cc) <- colnames(dat);
  mm <- max(dat);mmr <- ceiling(mm*10)/10; mmr1 <- 10*mmr
  par(mar=mar);a <- barplot(t(dat),beside = T,col=cc,yaxt='n',xaxt='n',ylim=c(0,mmr),ylab=ylab,
                            cex.axis = 1.2,cex.names = 1.4,main=main,cex.lab=1.4);
  pp<-par()$usr
  axis(side=2,at=c(0:mmr1)/10,labels=sprintf('%s%s',c(0:mmr1)*10,'%'),las=2,xpd=TRUE)
  if(ncol(a)>1) text(x=colMeans(a),pp[3]-(pp[4]-pp[3])/30,rownames(dat),srt=xtext_srt,xpd=T,adj=xtext_adj,cex=1.2)
  if(nrow(a)>1) legend(legend_pos,fill=cc,legend=names(cc),bty=bty,...)
}
get_PvCoxme <- function(x, rcoef=FALSE, digits=options()$digits, ...) {
  cat("Cox mixed-effects model fit by maximum likelihood\n")
  if (!is.null(x$call$data)) 
    cat("  Data:", deparse(x$call$data))
  if(!is.null(x$call$subset)) {
    cat(";  Subset:", deparse(x$call$subset), sep="\n")
  }
  else cat("\n")
  
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail<- nrow(x$var) - nvar
  
  omit <- x$na.action
  cat("  events, n = ", x$n[1], ', ', x$n[2], sep='')
  if(length(omit))
    cat(" (", naprint(omit), ")", sep = "")
  loglik <- x$loglik + c(0,0, x$penalty)
  temp <- matrix(loglik, nrow=1)
  cat("\n  Iterations=", x$iter, "\n")
  dimnames(temp) <- list("Log-likelihood", 
                         c("NULL", "Integrated", "Fitted"))
  print(temp)
  cat("\n")
  chi1 <- 2*diff(x$loglik[c(1,2)]) 
  
  
  chi1 <- 2*diff(loglik[1:2]) 
  chi2 <- 2*diff(loglik[c(1,3)])
  temp <- rbind(c(round(chi1,2), round(x$df[1],2),
                  signif(1- pchisq(chi1,x$df[1]),5),
                  round(chi1- 2*x$df[1],2),
                  round(chi1- log(x$n[1])*x$df[1],2)),
                c(round(chi2,2), round(x$df[2],2),
                  signif(1- pchisq(chi2,x$df[2]),5),
                  round(chi2- 2*x$df[2],2),
                  round(chi2- log(x$n[1])*x$df[2],2)))
  dimnames(temp) <- list(c("Integrated loglik", " Penalized loglik"),
                         c("Chisq", "df", "p", "AIC", "BIC"))
  print(temp, quote=F, digits=digits)
  
  cat ("\nModel: ", deparse(x$call$formula), "\n")
  
  if (nvar > 0)  { # Not a ~1 model
    se <- sqrt(diag(x$var)[nfrail+1:nvar])
    tmp <- cbind(beta, exp(beta), se, round(beta/se,2),
                 signif(1 - pchisq((beta/ se)^2, 1), 2))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                         "se(coef)", "z", "p"))
  }
  if (rcoef) { # print the random coefs
    #next line unlists, trying to give good names to the coefs
    coef <- unlist(lapply(ranef(x), function(y) {
      if (is.matrix(y)) {
        z <- c(y)
        dd <- dimnames(y)
        names(z) <- c(outer(dd[[1]], dd[[2]], paste,sep=':'))
        z}
      else y
    }))
    
    se <- sqrt(diag(x$var)[1:nfrail])
    rtmp <- cbind(coef, exp(coef), se)
    dimnames(rtmp) <- list(names(coef), c("coef", "exp(coef)",
                                          "Penalized se"))
  }
  
  if (nvar>0 && rcoef) {
    cat("Fixed and penalized coefficients\n")
    print(rbind(tmp, cbind(rtmp,NA,NA)), na.print='', digits=digits)
  }
  else if (rcoef) {
    cat("Penalized coefficients\n")
    print(rtmp, digits=digits)
  }
  else if (nvar>0) {
    cat("Fixed coefficients\n")
    print(tmp, digits=digits)
    pv=tmp
  }
  
  cat("\nRandom effects\n")
  
  random <- VarCorr(x)
  nrow <-  sapply(random, 
                  function(x) if (is.matrix(x)) nrow(x) else length(x))
  maxcol <-max(sapply(random,
                      function(x) if (is.matrix(x)) 1+ncol(x) else 2))
  temp1 <- matrix(NA, nrow=sum(nrow), ncol=maxcol)
  indx <- 0
  for (term in  random) {
    if (is.matrix(term)) {
      k <- nrow(term)
      nc <- ncol(term)  #assume nc > nr (only cases I know so far)
      for (j in 1:k) {
        temp1[j+indx, 1] <- sqrt(term[j,j])
        temp1[j+indx, 2] <- term[j,j]
        if (nc>j) {
          indx2 <- (j+1):nc
          temp1[j+indx, 1+ indx2] <- term[j, indx2]
        }
      }
    }
    else {
      k <- length(term)
      temp1[1:k + indx,1] <- sqrt(term)
      temp1[1:k + indx,2] <- term
    }
    indx <- indx + k
  }
  
  indx <- cumsum(c(1, nrow))   # starting row of each effect
  temp3 <- rep("", nrow(temp1))
  temp3[indx[-length(indx)]] <- names(random)
  xname <- unlist(lapply(random, 
                         function(x) if (is.matrix(x)) dimnames(x)[[1]] else names(x)))
  temp <- cbind(temp3, xname, ifelse(is.na(temp1), "", 
                                     format(temp1, digits=digits)))
  if (maxcol == 2)
    temp4 <- c("Group", "Variable", "Std Dev", "Variance")
  else 
    temp4 <- c("Group","Variable", "Std Dev", "Variance", "Corr", 
               rep("", maxcol-3))
  dimnames(temp) <- list(rep("", nrow(temp)), temp4)
  
  print(temp, quote=F)
  invisible(x)
  return(pv)
}
plot_SurvwithEffect<-function(dat,time_col='TestNegative-Onset(day)',
                              group_col='group',
                              censor_col='Censor',main='',effect_col='family_id',plot=TRUE){
  dat1 <- data.frame(time=dat[,time_col],status=dat[,censor_col],group=dat[,group_col],effect=dat[,effect_col])
  dat1 <- dat1[which(is.na(dat1$effect)==FALSE),]
  fit0 <- survfit(Surv(time,status)~group,data=dat1)
  fit1 <- coxph(Surv(time,status)~group,data=dat1)
  fit2 <- coxme(Surv(time,status)~group+(1|effect),data=dat1)
  pv_effect <- anova(fit1,fit2)[2,4];#print(anova(fit1,fit2));print(anova(fit2))
  logrank.test <- survdiff(Surv(time,status)~group,data=dat1);
  pv <- 1-pchisq(logrank.test$chisq,df=length(logrank.test$n)-1)
  get_PvCoxme(fit2)->tmp1;
  pv_fixEffect <- anova(fit2)[2,4];pv_fixEffect_c <- signif(pv_fixEffect,3); 
  if(pv_fixEffect==0) pv_fixEffect_c <- '<2.2e-16'
  if(plot==TRUE){
    plot1=ggsurvplot(fit0,data=dat1,risk.table = TRUE,pval=pv_fixEffect_c,conf.int=TRUE,palette='jama')
    plot1$plot$labels$y <- 'Test Positive probability';
    plot1$plot$labels$x <- 'Time to Onset (Day)';
    plot1$plot$labels$title <- sprintf('%s, N=%s',main,nrow(dat1))
  }else{
    plot1=NA;
  }
  return(list(plot=plot1,test=pv,test_effect=pv_fixEffect,effect_pv=pv_effect,group_pv=tmp1[,5]))
}
plot_SurvwithEffectTWO<-function(dat,time_col='TestNegative-Onset(day)',
                                 group_col='group',
                                 censor_col='Censor',main='',effect_col='family_id',effect_col2='family_id',plot=TRUE){
  dat1 <- data.frame(time=dat[,time_col],status=dat[,censor_col],group=dat[,group_col],
                     effect=dat[,effect_col],effect2=dat[,effect_col2])
  dat1 <- dat1[which(is.na(dat1$effect)==FALSE & is.na(dat1$effect2)==FALSE),]
  fit0 <- survfit(Surv(time,status)~group,data=dat1)
  fit1 <- coxph(Surv(time,status)~group,data=dat1)
  fit2 <- coxme(Surv(time,status)~group+(1|effect)+(1|effect2),data=dat1)
  pv_effect <- anova(fit1,fit2)[2,4];#print(anova(fit1,fit2));print(anova(fit2))
  logrank.test <- survdiff(Surv(time,status)~group,data=dat1);
  pv <- 1-pchisq(logrank.test$chisq,df=length(logrank.test$n)-1)
  get_PvCoxme(fit2)->tmp1;
  pv_fixEffect <- anova(fit2)[2,4];pv_fixEffect_c <- signif(pv_fixEffect,3); 
  if(pv_fixEffect==0) pv_fixEffect_c <- '<2.2e-16'
  if(plot==TRUE){
    plot1=ggsurvplot(fit0,data=dat1,risk.table = TRUE,pval=pv_fixEffect_c,conf.int=TRUE,palette='jama')
    plot1$plot$labels$y <- 'Test Positive probability';
    plot1$plot$labels$x <- 'Time to Onset (Day)';
    plot1$plot$labels$title <- sprintf('%s, N=%s',main,nrow(dat1))
  }else{
    plot1=NA;
  }
  return(list(plot=plot1,test=pv,test_effect=pv_fixEffect,effect_pv=pv_effect,group_pv=tmp1[,5]))
}
ind2family <- function(x,family_col='family_id',use_col=colnames(x)){
  x_family <- lapply(unique(x[,family_col]),function(xx){
    x[which(x[,family_col]==xx),use_col,drop=F]
  })
  names(x_family) <- unique(x[,family_col])
  return(x_family)
}
## 家长孩子满足某个变量相等
filterbyfamily <- function(dat,judge_col,use_col=c('TestNegative-Onset(day)','disease_group','Parent_or_Child','family_id','Vaccines','individualID')){
  dat_family <- ind2family(dat,use_col=use_col)
  tmp1 <- lapply(dat_family,function(x){
    x1 <- x[which(x$Parent_or_Child=='Parent'),];x2 <- x[which(x$Parent_or_Child=='Child'),]
    x3 <- intersect(x1[,judge_col],x2[,judge_col])
    if(length(x3)==0) return(NULL)
    x11 <- x1[which(x1[,judge_col] %in% x3),];x22 <- x2[which(x2[,judge_col] %in% x3),]
    rbind(x11,x22)
  })
  tmp11 <- do.call(rbind,tmp1);rownames(tmp11) <- tmp11$individualID
  message(sprintf('%s individuals passed the inclusion criteria with total=%s',nrow(tmp11),nrow(dat)))
  return(tmp11)
}
filterbyFamilyID <- function(x,use_id){
  x[which(x$family_id %in% use_id),]
}
plot_Pie <- function(x,...){
  t1 <- table(x); p1<-round(1000*t1/sum(t1))/10
  t2 <- sprintf('%s (N=%s,%s%s)',names(t1),t1,p1,'%')
  pie(t1,labels=t2,col=cc,...)
}
Q5basic1 <- function(dat_child,group_col = 'Age_category'){
  t1 <- plot_BarPre(dat_child,group_col = group_col)
  plot_Bar(t(t1),main='Symptoms for Children')
  plot_Bar(t1,main='Symptoms for Children',bty='n',xtext_srt = 0,xtext_adj = 0.5)
  t1 <- plot_BarPre(dat_child,group_col = group_col,use_col=use_col2);colnames(t1) <- use_col1
  plot_Bar(t(t1),main='Symptoms in Course for Children')
  plot_Bar(t1,main='Symptoms in Course for Children',bty='n',xtext_srt = 0,xtext_adj = 0.5)
  t1 <- plot_BarPre_simple(dat_child,group_col = group_col,use_col = 'WithBasicDisease')
  plot_Bar(t(t1),main='WithBasicDisease for Children',xtext_srt = 0,xtext_adj = 0.5)
  plot_Bar(t1,main='WithBasicDisease for Children',xtext_srt = 0,xtext_adj = 0.5,bty='n')
  t1 <- plot_BarPre_simple(dat_child,group_col = group_col,use_col = 'disease_group')
  plot_Bar(t(t1),main='disease_group for Children',xtext_srt = 0,xtext_adj = 0.5)
  plot_Bar(t1,main='disease_group for Children',xtext_srt = 0,xtext_adj = 0.5,bty='n')
  t1 <- plot_BarPre_simple(dat_child,group_col = group_col,use_col = 'Vaccines')
  plot_Bar(t(t1),main='Vaccines for Children',xtext_srt = 0,xtext_adj = 0.5)
  plot_Bar(t1,main='Vaccines for Children',xtext_srt = 0,xtext_adj = 0.5,bty='n')
  t1 <- plot_BarPre_simple(dat_child,group_col = group_col,use_col = 'Vaccines_YN')
  plot_Bar(t(t1),main='Vaccines for Children',xtext_srt = 0,xtext_adj = 0.5)
  plot_Bar(t1,main='Vaccines for Children',xtext_srt = 0,xtext_adj = 0.5,bty='n')
  t1 <- plot_BarPre_simple(dat_child,group_col = group_col,use_col = 'TestNegative-Onset(day)')
  cc1 <- colorRampPalette(brewer.pal(8,'Spectral'))(ncol(t1)); names(cc1) <- colnames(t1)
  plot_Bar(t1,main='Test Negative to onset days for Children',xtext_srt = 0,xtext_adj = 0.5,ncol=4,cc=cc1)
  cc1 <- colorRampPalette(brewer.pal(8,'Spectral'))(nrow(t1)); names(cc1) <- rownames(t1)
  plot_Bar(t(t1),main='Test Negative to onset days for Children',xtext_srt = 0,xtext_adj = 0.5,ncol=4,cc=cc1)
  boxplot(dat_child$`TestNegative-Onset(day)`~dat_child[,group_col],xlab=group_col,ylab='Test Negative to onset days for Children',
          cex.lab=1.4,cex.axis=1.2)
}
Q5basic2 <- function(dat_child,group_col = 'Age_category',effect_col=c('Age_category','Gender','WithBasicDisease','Vaccines','Vaccines_YN','disease_group',use_col1,use_col2)){
  #
  all_pv <- list(); res_plot<-list()
  tmp1=plot_Surv(dat_child,group_col =group_col,main=group_col)
  #print(tmp1$plot,newpage=FALSE);
  res_plot[[sprintf('%s_survival',group_col)]]=tmp1$plot
  all_pv[[sprintf('%s_survival',group_col)]] <- c(test_pv=tmp1$test,effect_pv=NA,group_pv=NA);
  #
  for(each_effect_col in effect_col){
    if(each_effect_col==group_col) next
    dat_child[,each_effect_col] <- removeNA_level(small2NA(dat_child[,each_effect_col]))
    if(min(table(dat_child[,each_effect_col]))<5) next
    print(each_effect_col)
    tmp1=plot_SurvwithEffect(dat_child,group_col = group_col,main=sprintf('%s with %s',group_col,each_effect_col),effect_col = each_effect_col,plot=FALSE)
    #print(tmp1$plot,newpage=TRUE);
    # plot for each
    u1 <- dat_child[,each_effect_col];tmp2<-c()
    for(j in unique(u1[which(is.na(u1)==F)])){
      print(j);w1 <- which(u1==j)
      if(min(table(dat_child[w1,group_col]))<5 | length(table(dat_child[w1,group_col]))<2) next
      tmp11=plot_Surv(dat_child[which(u1==j),],group_col = group_col,main=sprintf('%s,%s=%s',group_col,each_effect_col,j))
      tmp2[j]<-tmp11$test;res_plot[[sprintf('%s_%s_%s_survival',group_col,each_effect_col,j)]]=tmp11$plot;
    }
    #
    all_pv[[sprintf('%s_survival_with%sEffect',group_col,each_effect_col)]] <- c(test_pv=tmp1$test_effect,effect_pv=tmp1$effect_pv,group_pv=tmp1$group_pv,sep_pv=tmp2);
  }
  #
  return(list(pv=all_pv,plot=res_plot))
}
Q5basic3 <- function(dat_all,group_col = 'group',effect_col='family_id',effect_col2=c('Age_category','Gender','WithBasicDisease','Vaccines','Vaccines_YN','disease_group')){
  #
  all_pv <- list(); res_plot<-list()
  tmp1=plot_SurvwithEffect(dat_all,group_col =group_col,main=group_col,effect_col=effect_col)
  #print(tmp1$plot,newpage=FALSE);
  res_plot[[sprintf('%s_survival',group_col)]]=tmp1$plot
  all_pv[[sprintf('%s_survival',group_col)]] <- c(test_pv=tmp1$test_effect,effect_pv=tmp1$effect_pv,group_pv=tmp1$group_pv);
  #
  for(each_effect_col in effect_col2){
    if(each_effect_col==group_col) next
    dat_all[,each_effect_col] <- removeNA_level(small2NA(dat_all[,each_effect_col]))
    if(min(table(dat_all[,each_effect_col]))<5) next
    print(each_effect_col)
    tmp1=plot_SurvwithEffectTWO(dat_all,group_col = group_col,main=sprintf('%s with %s',group_col,each_effect_col),effect_col=effect_col,effect_col2 = each_effect_col,plot=FALSE)
    #print(tmp1$plot,newpage=TRUE);
    # plot for each
    u1 <- dat_all[,each_effect_col];tmp2<-c()
    for(j in unique(u1[which(is.na(u1)==F)])){
      print(j);w1 <- which(u1==j)
      if(min(table(dat_all[w1,group_col]))<5 | length(table(dat_all[w1,group_col]))<2) next
      tmp11=plot_SurvwithEffect(dat_all[which(u1==j),],group_col = group_col,main=sprintf('%s,%s=%s',group_col,each_effect_col,j),effect_col=effect_col)
      tmp2[j]<-tmp11$test;res_plot[[sprintf('%s_%s_%s_survival',group_col,each_effect_col,j)]]=tmp11$plot;
    }
    #
    all_pv[[sprintf('%s_survival_with%sEffect',group_col,each_effect_col)]] <- c(test_pv=tmp1$test_effect,effect_pv=tmp1$effect_pv,group_pv=tmp1$group_pv,sep_pv=tmp2);
  }
  #
  return(list(pv=all_pv,plot=res_plot))
}


table1_getPv <- function(out_tab,K=2,ori_dat,group_col='group'){
  out2 <- as.data.frame(out_tab)
  # add p-value
  bg_N <- as.numeric(gsub('\\(N=(\\d+)\\)','\\1',out2[1,]));all_pv <- list()
  # prop test
  w1 <- grep('%',out2[,2]);
  for(k in c(1:K)*2){
    pv<-c()
    for(i in w1){
      p1 <- as.numeric(gsub('(\\d+) \\(.*','\\1',out2[i,k:(k+1)]));
      pv1 <- fisher.test(cbind(p1,bg_N[k:(k+1)]-p1))$p.value;pv[i] <- signif(pv1,3)
    }
    all_pv[[k/2]] <- pv
  }
  # t.test for mean
  bg_C <- colnames(out2)
  w1 <- grep('Mean',out2[,1]);
  for(k in c(1:K)*2){
    for(i in w1){
      v1 <- out2[i-1,1];w2 <- which(ori_dat[,group_col] %in% bg_C[k:(k+1)])
      r1 <- get_datG(ori_dat[w2,],group_col = group_col,day_col = v1);r2 <- unlist(lapply(r1,length))
      if(min(r2)<3){
        pv1<-1
      }else{
        pv1 <- t.test(ori_dat[w2,v1]~ori_dat[w2,group_col])$p.value
      }
      all_pv[[k/2]][i] <- pv1
    }
  }
  # wilcox for median
  bg_C <- colnames(out2)
  w1 <- grep('Median',out2[,1]);
  for(k in c(1:K)*2){
    for(i in w1){
      v1 <- out2[i-2,1];w2 <- which(ori_dat[,group_col] %in% bg_C[k:(k+1)])
      pv1 <- wilcox.test(ori_dat[w2,v1]~ori_dat[w2,group_col])$p.value
      all_pv[[k/2]][i] <- pv1
    }
  }
  #for(i in 1:K){all_pv[[i]]<-p.adjust(all_pv[[i]])}
  out3 <- c('class'=out2[,1])
  for(i in 1:K){
    out3 <- cbind(out3,out2[,c(2*i):c(2*i+1)],`P-value`=all_pv[[i]])
  }
  out3 <- cbind(out3,total=out2[,ncol(out2)]);colnames(out3)[1]<-'Class'
  return(out3)
}

get_datG <- function(dat=use_dat_pair1,day_col='day',group_col='group',effect_col='family_id'){
  tmp1 <- aggregate(dat[,day_col],list(dat[,group_col]),c,simplify=F);x<-tmp1$x;names(x)<-tmp1$Group.1;
  x
}
get_pPair <- function(x){
  pv<-c();nn<-names(x);
  for(i in 1:(length(x)-1)){
    for(j in (i+1):length(x)){
      if(nn[i]=='Child' & nn[j] == 'Parent'){
        pv[sprintf('%s_%s',nn[i],nn[j])]<-t.test(x[[i]],x[[j]],paired=TRUE)$p.value
      }else{
        pv[sprintf('%s_%s',nn[i],nn[j])]<-t.test(x[[i]],x[[j]])$p.value
      }
    }
  }
  return(pv)
}
get_summ <- function(x){
  data.frame(names(x),N=sapply(x,length),mean=sapply(x,mean),sd=sapply(x,sd),median=sapply(x,median),min=sapply(x,min),max=sapply(x,max))
}

get_CI <- function(dat,y_col,x_col,use_col=cc[4],use_sd=1.64,do_test=FALSE){
  v1 <- lowess(y=dat[,y_col],x=dat[,x_col]);
  lines(v1, col=use_col,lwd=2);
  tmp1 <- aggregate(dat[,y_col],list(dat[,x_col]),c)
  tmp2 <- tmp1$x; names(tmp2) <- tmp1$Group.1
  tmp2m <- unlist(lapply(tmp2,function(x)mean(x,na.rm=TRUE)))
  tmp2s <- unlist(lapply(tmp2,function(x)sd(x,na.rm=TRUE)));tmp2s[which(is.na(tmp2s)==TRUE)] <- 0
  v12 <- lowess(y=tmp2m+use_sd*tmp2s,x=tmp1$Group.1);
  v13 <- lowess(y=tmp2m-use_sd*tmp2s,x=tmp1$Group.1);
  polygon(x=c(v12$x,rev(v13$x)),y=c(v12$y,rev(v13$y)),col=adjustcolor(use_col,0.3),border=NA)
  if(do_test==TRUE){
    pv <- unlist(lapply(tmp2,function(x){
      if(length(x)>3) t.test(x)$p.value else 1
    }))
    pv1 <- p.adjust(pv)
    return(cbind(mean=tmp2m,sd=tmp2s,pvalue=pv,padjust=pv1))
  }
}
get_CIMod <- function(dat,y_col,x_col,use_col=cc[4],use_sd=1.64,do_test=FALSE,plot=TRUE){
  v1 <- lowess(y=dat[,y_col],x=dat[,x_col]);
  if(plot==TRUE) lines(v1, col=use_col,lwd=2);
  tmp1 <- aggregate(dat[,y_col],list(dat[,x_col]),c)
  tmp2 <- tmp1$x; names(tmp2) <- tmp1$Group.1
  pv <- lapply(tmp2,function(x){
    if(length(x)>3) t.test(x) else NULL
  })
  pv_conf <- lapply(pv,function(x)x$conf.int)
  tmp2 <- do.call(rbind,pv_conf);tmp2 <- tmp2[which(is.nan(tmp2[,1])==F),]
  #tmp2m <- unlist(lapply(tmp2,function(x)mean(x,na.rm=TRUE)))
  #tmp2s <- unlist(lapply(tmp2,function(x)sd(x,na.rm=TRUE)));tmp2s[which(is.na(tmp2s)==TRUE)] <- 0
  #v12 <- lowess(y=tmp2m+use_sd*tmp2s,x=tmp1$Group.1);
  #v13 <- lowess(y=tmp2m-use_sd*tmp2s,x=tmp1$Group.1);
  v12 <- lowess(y=tmp2[,1],x=rownames(tmp2))
  v13 <- lowess(y=tmp2[,2],x=rownames(tmp2))
  if(plot==TRUE) polygon(x=c(v12$x,rev(v13$x)),y=c(v12$y,rev(v13$y)),col=adjustcolor(use_col,0.3),border=NA)
  if(do_test==TRUE){
    pv1 <- unlist(lapply(pv[as.character(rownames(tmp2))],function(x){x$p.value}))
    pv2 <- p.adjust(pv1)
    return(cbind(CI_low=tmp2[,1],CI_high=tmp2[,2],pvalue=pv1,padjust=pv2))
  }
}




