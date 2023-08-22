## X: binary, category, numeric
## Y: binary, category, numeric
#source('D:/analysis_eng/SherryPlot/function_ML.R')
test_feature <- function(input_dataset,
                         input_featureType=NULL,
                         X_columns=NULL,
                         Y_columns=NULL,
                         p.adjust.method='fdr',
                         plot=FALSE,remove_character=NULL,min_count=3){
  ##
  if(is.null(input_featureType)==TRUE){
    input_featureType <- auto_dataType(input_dataset=input_dataset,
                                       ID_column=NULL,
                                       factorize=FALSE,
                                       remove_identical=FALSE,
                                       remove_unique=FALSE,
                                       remain_column=NULL,
                                       return_dataset=FALSE,
                                       return_featureType=TRUE)
  }
  ##
  use_featureType <- c('binary','category','numeric','binary_withNA','category_withNA','numeric_withNA')
  all_columns <- names(input_featureType)[which(input_featureType %in% use_featureType)]
  ##
  X_columns <- intersect(X_columns,all_columns)
  Y_columns <- intersect(Y_columns,all_columns)
  if(length(Y_columns)==0 | length(X_columns)==0){
    stop('Please check input columns!')
  }
  ## test-statistics, p-value, value for plot
  all_res <- list()
  for(each_Y in Y_columns){
    ori_Y <- input_dataset[,each_Y]
    for(each_X in X_columns){
      print(each_X)
      ori_X <- input_dataset[,each_X]
      if(is.null(remove_character)==FALSE){
        u1 <- which(!ori_X %in% remove_character)
        X <- ori_X[u1];Y<-ori_Y[u1]
        tmpX <- table(X);
        u2 <- which(X %in% names(tmpX)[which(tmpX>=min_count)])
        X <- X[u2];Y<-Y[u2]
        if(length(u2)<=3) next
        if(length(unique(X))<=1) next
        if(length(unique(Y))<=1) next
      }else{
        X <- ori_X;Y<-ori_Y
      }
      ## binary, category, numeric
      XF <- input_featureType[each_X]
      YF <- input_featureType[each_Y]
      if(grepl('numeric',XF) & grepl('numeric',YF)){
          # correlation
		X <- as.numeric(X)
		Y <- as.numeric(Y)
        res1 <- cor.test(X,Y,method='pearson');
        each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='pearson',p.value=res1$p.value,statistic=as.numeric(res1$statistic))
        #res1 <- cor.test(X,Y,method='spearman')
        #each_res <- list(method='spearman',p.value=res1$p.value,statistic=res1$statistic)
      }else{
        if(grepl('numeric',XF) & !grepl('numeric',YF)){
		  X <- as.numeric(X)
          # t-test/anova
          if(grepl('binary',YF)){
		  	tryres <- try({res1 <- t.test(X~Y)})
            if (class(tryres)!='try-error'){
          		each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='t-test',p.value=res1$p.value,statistic=as.numeric(res1$statistic))
			}else{
          		each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='t-test',p.value=NA,statistic=NA)
			}
		  }else{
            tryres <- try({res1 <- aov(X~Y)})
            if (class(tryres)!='try-error'){
          		each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='ANOVA',p.value=summary(res1)[[1]][1,5],statistic=summary(res1)[[1]][1,4])
			}else{
          		each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='ANOVA',p.value=NA,statistic=NA)
			}
		  }
        }
        if(!grepl('numeric',XF) & grepl('numeric',YF)){
		  Y <- as.numeric(Y)
          # t-test/anova
          if(grepl('binary',XF)){
		  	tryres <- try({res1 <- t.test(Y~X)})
            if (class(tryres)!='try-error'){
            	each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='t-test',p.value=res1$p.value,statistic=as.numeric(res1$statistic))
			}else{
          		each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='t-test',p.value=NA,statistic=NA)
			}
		  }else{
            tryres <- try({res1 <- aov(Y~X)})
			if (class(tryres)!='try-error'){
            	each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='ANOVA',p.value=summary(res1)[[1]][1,5],statistic=summary(res1)[[1]][1,4])
			}else{
				each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='ANOVA',p.value=NA,statistic=NA)
			}
		  }
        }
        if(!grepl('numeric',XF) & !grepl('numeric',YF)){
          # chisq
          t1 <- table(list(X,Y))
          tryres <- try({res1 <- chisq.test(t1)})
		  if(class(tryres) != 'try-error'){
              each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='chisq-test',p.value=res1$p.value,statistic=as.numeric(res1$statistic))
		  }else{
              each_res <- list(X_var=each_X,Y_var=each_Y,X=X,Y=Y,method='chisq-test',p.value=NA,statistic=NA)
		  }
        }
        all_res[[sprintf('%s-%s',each_X,each_Y)]] <- each_res
      }
    }
  }
  ##
  all_X_var <- unlist(lapply(all_res,function(x)x$X_var))
  all_Y_var <- unlist(lapply(all_res,function(x)x$Y_var))
  all_pv <- unlist(lapply(all_res,function(x)x$p.value))
  all_sta <- unlist(lapply(all_res,function(x)x$statistic))
  tab_res <- data.frame(X_var=all_X_var,Y_var=all_Y_var,p.value=all_pv,statistics=all_sta)
  tab_res$adjust.p <- p.adjust(tab_res$p.value,method=p.adjust.method)
  tab_res <- tab_res[order(tab_res$p.value),]
  return(list(table_result=tab_res,all_result=all_res))
}
