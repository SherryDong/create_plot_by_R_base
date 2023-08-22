### functions
library(caret)
library(ordinal)
library(pROC)
library(gbm)
library(glmnet)
library(ROCR)
library(RColorBrewer)
library(MLmetrics)
##
# unique,identical, binary, category, numeric
# unique_withNA,identical_withNA, binary_withNA, category_withNA, numeric_withNA
## load demo dataset
load_demoDataset <- function(dataset='metabric'){
  # demo dataset:https://www.cbioportal.org/study/summary?id=brca_metabric
  if(dataset == 'metabric'){
    d1 <- read.delim('data/brca_metabric_clinical_data.tsv',stringsAsFactors = F)
  }
  return(d1)
}
## auto-define data type
auto_dataType <- function(input_dataset,
                          ID_column=NULL,
                          factorize=TRUE,
                          remove_identical=TRUE,
                          remove_unique=TRUE,
                          remain_column=NULL,
                          return_dataset=TRUE,
                          return_featureType=TRUE){
  # unique,identical, binary, category, numeric
  # unique_withNA,identical_withNA, binary_withNA, category_withNA, numeric_withNA
  input_dataset <- unique(input_dataset);
  all_feature <- colnames(input_dataset)
  all_featureType <- c();N <- nrow(input_dataset)
  for(i in all_feature){
    x1 <- input_dataset[,i]; x1 <- as.character(x1)
    x1_noNA <- x1[which(is.na(x1)==F)]
    x1_uniq <- unique(x1);x1_noNA_uniq <- unique(x1_noNA)
    if(length(x1_uniq)==1){
      all_featureType[i]<-'identical' ## only one identical value
      next
    }
    if(length(x1_uniq)==N){
      all_featureType[i]<-'unique' ## each with unique value
      next
    }
    if(length(x1_noNA_uniq)==1){
      all_featureType[i]<-'identical_withNA' ## only one identical value with NA
      next
    }
    if(length(x1_noNA_uniq)==N){
      all_featureType[i]<-'unique_withNA' ## each with unique value with NA
      next
    }
    if(length(x1_uniq)==2){
      all_featureType[i]<-'binary' ## binary
      next
    }
    if(length(x1_noNA_uniq)==2){
      all_featureType[i]<-'binary_withNA' ## binary
      next
    }
    # numeric,category,NA
    x1_tryNum <- as.numeric(x1);
    x1_tryNum_noNA <- x1_tryNum[which(is.na(x1_tryNum)==F)]
    x1_tryNum_uniq <- unique(x1_tryNum);x1_tryNum_noNA_uniq <- unique(x1_tryNum_noNA)
    if(length(x1_tryNum_noNA_uniq)<length(x1_noNA_uniq)){ ## not all numeric
      if(length(x1_uniq)==length(x1_noNA_uniq)){
        all_featureType[i]<-'category' 
      }else{
        all_featureType[i]<-'category_withNA'
      }
      if(length(x1_noNA_uniq)>N/10){
        message(sprintf('%s column has %d unique non-NA values, please check!',i,
                        length(x1_noNA_uniq)))
      }
      next
    }else{
      if(length(x1_tryNum_uniq)==length(x1_tryNum_noNA_uniq)){
        all_featureType[i]<-'numeric' 
      }else{
        all_featureType[i]<-'numeric_withNA'
      }
    }
  }
  #
  output_dataset <- input_dataset; 
  # ID_column=NULL
  candidate_ID_column <- names(all_featureType)[which(all_featureType=='unique')]
  if(length(candidate_ID_column)==0){
    message('No unique column, will set fake rownames for the dataset!')
    rownames(output_dataset) <- sprintf('Sample:%s',1:nrow(output_dataset))
    use_ID_column <- NULL
  }else{
    if(is.null(ID_column)==TRUE){
      use_ID_column <- candidate_ID_column[1]
    }else{
      if(ID_column %in% candidate_ID_column){
        use_ID_column <- ID_column
      }else{
        use_ID_column <- candidate_ID_column[1]
        message(sprintf('%s is not unique, will use %s instead !',ID_column,use_ID_column))
      }
    }
    message(sprintf('Use %s as ID column !',use_ID_column))
    rownames(output_dataset) <- as.character(input_dataset[,use_ID_column])
  }
  # factorize=FALSE
  if(factorize==TRUE){
    use_featureType <- c('binary','binary_withNA','category','category_withNA')
    use_feature <- names(all_featureType)[which(all_featureType %in% use_featureType)]
    for(i in use_feature){
      output_dataset[,i] <- factor(output_dataset[,i])
    }
  }
  # # unique,identical, binary, category, numeric
  # unique_withNA,identical_withNA, binary_withNA, category_withNA, numeric_withNA
  # remove_identical=TRUE
  if(remove_identical==TRUE){
    use_featureType <- c('identical_withNA','identical')
    use_feature <- names(all_featureType)[which(!all_featureType %in% use_featureType)]
    rc <- intersect(use_feature,colnames(output_dataset))
    rc <- unique(rc,remain_column)
    output_dataset <- output_dataset[,rc,drop=F]
  }
  # remove_unique=TRUE
  if(remove_unique==TRUE){
    use_featureType <- c('unique_withNA','unique')
    use_feature <- names(all_featureType)[which(!all_featureType %in% use_featureType)]
    rc <- intersect(use_feature,colnames(output_dataset))
    rc <- unique(c(rc,remain_column))
    rc <- unique(c(use_ID_column,rc))
    output_dataset <- output_dataset[,rc,drop=F]
  }
  #output
  output_featureType <- all_featureType[colnames(output_dataset)]
  if(return_dataset==TRUE & return_featureType==TRUE){
    return(list(dataset=output_dataset,featureType=output_featureType))
  }
  if(return_dataset==TRUE & return_featureType==FALSE){
    return(output_dataset)
  }
  if(return_dataset==FALSE & return_featureType==TRUE){
    return(output_featureType)
  }
}
## auto transfer numeric (with small unique value to category)
auto_num2Category <- function(input_dataset,
                              check_column=NULL,
                              min_N=3,
                              force_column=NULL,factorize=TRUE){
  input_dataset <- unique(input_dataset);
  all_feature <- colnames(input_dataset);
  output_dataset <- input_dataset; 
  if(is.null(check_column)==FALSE){
    use_column <- intersect(all_feature,check_column)
  }else{
    check_column <- all_feature
  }
  if(length(use_column)==0){error('No remained column!');}
  for(i in use_column){
    x1 <- input_dataset[,i]; x1_noNA <- x1[which(is.na(x1)==F)]
    x1_uniq <- unique(x1);x1_noNA_uniq <- unique(x1_noNA)
    x1_tryNum <- as.numeric(x1);
    x1_tryNum_noNA <- x1_tryNum[which(is.na(x1_tryNum)==F)]
    x1_tryNum_uniq <- unique(x1_tryNum);
    x1_tryNum_noNA_uniq <- unique(x1_tryNum_noNA)
    if(length(x1_tryNum_noNA_uniq)<=min_N | i %in% force_column){
      message(sprintf('%s is transformed into category!',i))
      output_dataset[,i] <- sprintf('Class%s',output_dataset[,i])
      if(factorize==TRUE){
        output_dataset[,i] <- factor(output_dataset[,i])
      }
    }
  }
  return(output_dataset)
  #
}
## pre-processing
#基于临床研究等实际业务的处理，如合并、填补、分组等
merge_feature <- function(input_dataset,
                          old_feature=NULL,
                          new_feature=NULL,
                          remove_oldFeature=TRUE,
                          sep='-',
                          factorize=TRUE){
  tmp_d <- apply(input_dataset[,old_feature],1,function(x)paste(x,collapse = sep))
  tmp_NA <- paste(rep('NA',length.out=length(old_feature)),collapse = sep)
  tmp_d[which(tmp_d==tmp_NA)] <- NA
  if(is.null(new_feature)==TRUE){
    tmp_f <- paste(old_feature,collapse = sep)
  }else{
    tmp_f <- new_feature
  }
  output_dataset <- input_dataset
  if(remove_oldFeature==TRUE){
    rc <- setdiff(colnames(output_dataset),old_feature)
    output_dataset <- output_dataset[,rc,drop=F]
  }
  if(factorize==TRUE){tmp_d <- factor(tmp_d)}
  output_dataset <- cbind(output_dataset,tmp_d);
  colnames(output_dataset)[ncol(output_dataset)] <- tmp_f
  return(output_dataset)
}
#连续变量进行 Box-Cox 变换、scaling 等变换无量纲化并转为近似标准高斯分布缺失处理。
# 例如
  #缺失率过高：变量删除（也包括 almost-constant 变量）
  #缺失率过低：均值/众数填补或其它填补
  #离散变量：缺失作为单独level
  #连续变量：均值/随机填补+新建 missing indicator 二值变量；mock数据生成后再合并回来。
## NA_remove_threPct: percentage to remove the column
## binary_NA_strategy:choose from removeByPct,imputeByMode,createNewCate
## category_NA_strategy:choose from removeByPct,imputeByMode,createNewCate
## order: imputeByMode>createNewCate
## numeric_NA_strategy: choose from removeByPct,imputeByMean>imputeByMedian
impute_data <- function(input_dataset,input_featureType=NULL,
                  binary_NA_strategy=c('removeByPct','imputeByMode'), 
                  category_NA_strategy=c('removeByPct','createNewCate'),
                  numeric_NA_strategy=c('removeByPct','imputeByMean'),
                  NA_remove_threPct=0.25){
  input_dataset <- unique(input_dataset);
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
  output_dataset <- input_dataset
  ## remove:unique_withNA,identical_withNA
  use_featureType <- c('unique_withNA','identical_withNA')
  use_feature <- names(input_featureType)[which(!input_featureType %in% use_featureType)]
  rc <- intersect(use_feature,colnames(output_dataset))
  output_dataset <- output_dataset[,rc,drop=F]
  ## binary_withNA/category_withNA/numeric_withNA (removeByPct,imputeByMode,createNewCate)
  # removeByPct
  N <- nrow(output_dataset)
  use_featureType <- c('numeric_withNA','binary_withNA','category_withNA')
  use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
  rc <- intersect(use_feature,colnames(output_dataset))
  rc_NA_count <- unlist(lapply(rc,function(x){length(which(is.na(output_dataset[,x])==TRUE))}))
  rc_NA_count_pct <- rc_NA_count/N;w1 <- which(rc_NA_count_pct>=NA_remove_threPct)
  if('removeByPct' %in% binary_NA_strategy){
    use_featureType <- c('binary_withNA')
    use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
    remove_rc <- intersect(rc[w1],use_feature)
    if(length(remove_rc)>0){
      message(sprintf('%s will be removed due to too NA values (%s >%s)',
                      paste(remove_rc,collapse=';'),
                      paste(signif(rc_NA_count_pct[w1],4),collapse=';'),NA_remove_threPct))
      output_dataset <- output_dataset[,setdiff(colnames(output_dataset),remove_rc),drop=F]
    }
  }
  if('removeByPct' %in% category_NA_strategy){
    use_featureType <- c('category_withNA')
    use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
    remove_rc <- intersect(rc[w1],use_feature)
    if(length(remove_rc)>0){
      message(sprintf('%s will be removed due to too NA values (%s >%s)',
                      paste(remove_rc,collapse=';'),
                      paste(signif(rc_NA_count_pct[w1],4),collapse=';'),NA_remove_threPct))
      output_dataset <- output_dataset[,setdiff(colnames(output_dataset),remove_rc),drop=F]
    }
  }
  if('removeByPct' %in% numeric_NA_strategy){
    use_featureType <- c('numeric_withNA')
    use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
    remove_rc <- intersect(rc[w1],use_feature)
    if(length(remove_rc)>0){
      message(sprintf('%s will be removed due to too NA values (%s >%s)',
                      paste(remove_rc,collapse=';'),
                      paste(signif(rc_NA_count_pct[w1],4),collapse=';'),NA_remove_threPct))
      output_dataset <- output_dataset[,setdiff(colnames(output_dataset),remove_rc),drop=F]
    }
  }
  # imputeByMode
  use_featureType <- c('binary_withNA','category_withNA')
  use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
  rc <- intersect(use_feature,colnames(output_dataset))
  rc_mode <- unlist(lapply(rc,function(x){names(sort(table(output_dataset[,x]),decreasing = T))[1]}))
  if(length(rc_mode)>0) names(rc_mode) <- rc
  if('imputeByMode' %in% binary_NA_strategy){
    use_featureType <- c('binary_withNA')
    use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
    use_feature <- intersect(use_feature,colnames(output_dataset))
    if(length(use_feature)>0){
      for(i in use_feature){
        message(sprintf('%s will be imputed by %s',i,rc_mode[i]))
        x1 <- output_dataset[,i];w1 <- which(is.na(x1)==TRUE)
        x1[w1] <- rc_mode[i];output_dataset[,i] <- x1
        input_featureType[i] <- 'binary'
      }
    }
  }
  if('imputeByMode' %in% category_NA_strategy){
    use_featureType <- c('category_withNA')
    use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
    use_feature <- intersect(use_feature,colnames(output_dataset))
    if(length(use_feature)>0){
      for(i in use_feature){
        message(sprintf('%s will be imputed by %s',i,rc_mode[i]))
        x1 <- output_dataset[,i];w1 <- which(is.na(x1)==TRUE)
        x1[w1] <- rc_mode[i];output_dataset[,i] <- x1
        input_featureType[i] <- 'category'
      }
    }
  }
  # createNewCate
  use_featureType <- c('binary_withNA','category_withNA')
  use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
  if('createNewCate' %in% binary_NA_strategy){
    use_featureType <- c('binary_withNA')
    use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
    use_feature <- intersect(use_feature,colnames(output_dataset))
    if(length(use_feature)>0){
      for(i in use_feature){
        message(sprintf('%s will add one missing class',i))
        x1 <- as.character(output_dataset[,i]);w1 <- which(is.na(x1)==TRUE)
        x1[w1] <- 'missing';output_dataset[,i] <- factor(x1)
        input_featureType[i] <- 'binary'
      }
    }
  }
  if('createNewCate' %in% category_NA_strategy){
    use_featureType <- c('category_withNA')
    use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
    use_feature <- intersect(use_feature,colnames(output_dataset))
    if(length(use_feature)>0){
      for(i in use_feature){
        message(sprintf('%s will add one missing class',i))
        x1 <- as.character(output_dataset[,i]);w1 <- which(is.na(x1)==TRUE)
        x1[w1] <- 'missing';output_dataset[,i] <- factor(x1)
        input_featureType[i] <- 'category'
      }
    }
  }
  ## numeric_withNA, 
  use_featureType <- c('numeric_withNA')
  use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
  use_feature <- intersect(use_feature,colnames(output_dataset))
  if('imputeByMean' %in% numeric_NA_strategy){
    if(length(use_feature)>0){
      for(i in use_feature){
        message(sprintf('%s will be imputed by mean',i))
        x1 <- output_dataset[,i];w1 <- which(is.na(x1)==TRUE)
        x1[w1] <- mean(x1,na.rm=T);output_dataset[,i] <- x1
        input_featureType[i] <- 'numeric'
      }
    }
  }
  if('imputeByMedian' %in% numeric_NA_strategy){
    if(length(use_feature)>0){
      for(i in use_feature){
        message(sprintf('%s will be imputed by median',i))
        x1 <- output_dataset[,i];w1 <- which(is.na(x1)==TRUE)
        x1[w1] <- median(x1,na.rm=T);output_dataset[,i] <- x1
        input_featureType[i] <- 'numeric'
      }
    }
  }
  ##
  return(output_dataset)
}

#######################################################################################
#2.3. 预测建模
#假设源数据有K个变量。对于每个变量，以其它变量作为自变量，该变量为预测目标构建预测模型。
#输出K个学习器 Lk, k=1,…,K
#具体实现方式：对于每个变量 x[k]，针对不同数据类型采取不同的模型
#binary: Logistic Regression，输出阳性概率
#categorical: Multi-class Logistic Regression，输出每个类别的概率
#numerical: Linear Regression，输出预测均值和方差
#说明
#使用线性模型的原因是，输出概率需要模型是 well-calibrated。
#也可以考虑ANN等最后一层为线性层的预测模型，或者tree-based model + calibrating 后处理。
#建模细节因素
#使用 CV/holdout 等方式调参
#自变量需要 one-hot encoding 等预处理
#使用单变量/Lasso/ElasticNet等常规方式进行特征选择
# unique,identical, binary, category, numeric
# ... other parameters into select_featureBywrapper
generate_predModel <- function(input_dataset,input_featureType=NULL,...){
  input_dataset <- unique(input_dataset);
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
  use_featureType <- c('binary','category','numeric')
  use_feature <- names(input_featureType)[which(input_featureType %in% use_featureType)]
  rc <- intersect(use_feature,colnames(input_dataset))
  use_featureType1 <- c('binary_withNA','category_withNA','numeric_withNA')
  use_feature1 <- names(input_featureType)[which(input_featureType %in% use_featureType1)]
  rc1 <- intersect(use_feature1,colnames(input_dataset))
  message(sprintf('%s features (%s) in binary, category, numeric category will be processed; %s features in binary_withNA, category_withNA, numeric_withNA (%s) will pass !',length(rc),paste(rc,collapse=';'),length(rc1),paste(rc1,collapse=';')))
  all_use_fit <- list()
  for(pred_feature in rc){
    print(sprintf('Process for %s',pred_feature))
    res1 <- get_XY(input_dataset=input_dataset,
                   input_featureType=input_featureType,
                   pred_feature=pred_feature)
    X <- res1$X; Y <- res1$Y; Y_class <- res1$Y_class
    err <- try(res1_p <- select_featureBywrapper(X=X,Y=Y,Y_class=Y_class,...))
    if(class(err)=='try-error'){
      print(sprintf('%s cannot be processed',pred_feature))
      next
    }
    all_use_fit[[pred_feature]] <- process_wrapperResult(result_list=res1_p,draw=FALSE)
  }
  return(all_use_fit)
}
# get X,Y
get_XY <- function(input_dataset,input_featureType=NULL,pred_feature=''){
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
  Y <- input_dataset[,pred_feature]
  # remove features cannot be used as predictors
  use_featureType <- c('unique','identical','unique_withNA','identical_withNA')
  use_feature <- names(input_featureType)[which(!input_featureType %in% use_featureType)]
  rc <- setdiff(intersect(use_feature,colnames(input_dataset)),pred_feature)
  X <- input_dataset[,rc,drop=F];names(Y) <- rownames(X)
  return(list(X=X,Y=Y,Y_class=input_featureType[pred_feature]))
}

# feature selection (binary/category/numeric)
# method_name (gbm/glmnet)
select_featureBywrapper <- function(X,Y,Y_class='binary',
                                    method_name='glmnet',
                                    cv_number=10){
  ## transform Y class
  use_Y <- Y; ## original Y
  if(Y_class=='numeric'){
    class_transform <- NA;
  }else{
    levels(use_Y) <- sprintf('class%d',1:length(levels(Y)));
    class_transform <- levels(Y); names(class_transform) <- levels(use_Y);
    sprintf('%s to %s',paste(levels(Y),collapse=';'),paste(levels(use_Y),collapse=';'))
  }
  use_X <- X; ## original X
  new_Xcolnames <- lapply(names(use_X),function(x)paste(x,levels(use_X[,x]),sep=''));
  names(new_Xcolnames) <- names(use_X);
  ##
  if(Y_class=='binary'){
    ctrl <- trainControl(method="repeatedcv",summaryFunction=twoClassSummary,
                         classProbs=T,savePredictions = T, number=cv_number) # number 10 means 10-fold cv
  }
  if(Y_class=='category'){
    ctrl <- trainControl(method="repeatedcv",summaryFunction=multiClassSummary,
                         classProbs=T,savePredictions = T, number=cv_number) # number 10 means 10-fold cv
  }
  if(Y_class=='numeric'){
    ctrl <- trainControl(method="repeatedcv",
                         classProbs=T,savePredictions = T, number=cv_number) # number 10 means 10-fold cv
  }
  all_p <- list()
  i=0;
  while(ncol(use_X)>1){
    i=i+1;
    print(dim(use_X))
    dat <- data.frame(use_X,Y=use_Y)
    if(Y_class=='numeric'){
      err <- try(fit1 <- train(Y ~ ., data=dat, method=method_name, trControl=ctrl))
      if(class(err[1])=='try-error'){
        stop('Model cannot be generated!')
      }
      #fit1 <- train(Y ~ ., data=dat, method=method_name, trControl=ctrl)
      predict_CA <- predict(fit1);
      f1 <- cor(use_Y,predict_CA);
      f2 <- cor(fit1$pred$obs,fit1$pred$pred)
    }else{
      err <- try(fit1 <- train(Y ~ ., data=dat, method=method_name, trControl=ctrl, metric = 'Accuracy'))
      if(class(err[1])=='try-error'){
        stop('Model cannot be generated!')
      }
      #fit1 <- train(Y ~ ., data=dat, method=method_name, trControl=ctrl, metric = 'Accuracy')
      predict_CA <- predict(fit1,type='prob');
      f1 <- multiclass.roc(use_Y,predict_CA)$auc[[1]];
      f2 <- multiclass.roc(fit1$pred$obs,fit1$pred[,levels(use_Y)])$auc[[1]];
    }
    feat_imp <- varImp(fit1, scale = TRUE)$importance;
    feat_imp <- feat_imp[order(feat_imp$Overall,decreasing = T),,drop=F]
    w1 <- which(feat_imp==min(feat_imp)) ##
    w2 <- which(feat_imp>min(feat_imp)) ##
    w3 <- lapply(new_Xcolnames,function(x)intersect(x,rownames(feat_imp)[w2]))
    w4 <- unlist(lapply(w3,length))
    w5 <- names(w4)[which(w4>0)]
    # next feature use: 
    new_feature <- w5; old_feature <- colnames(use_X)
    all_p[[i]] <- list(AUROC=f1,AUROC_cv=f2,fit=fit1,feature_importance=feat_imp[w2,,drop=F],
                       use_feature=old_feature)
    # next feature use: 
    if(length(setdiff(new_feature,old_feature))==0 & length(setdiff(old_feature,new_feature))==0){
      break
    }else{
      use_X <- use_X[,new_feature,drop=F]
    }
  }
  ## return
  return(list(all_result=all_p,class_transform=class_transform,method_name=method_name,Y_class=Y_class))
}

## deal with result from select_featureBywrapper
#list(AUROC=f1,AUROC_cv=f2,fit=fit1,feature_importance=feat_imp[w2,,drop=F],use_feature=w5)
#return(list(all_result=all_p,class_transform=class_transform))
process_wrapperResult <- function(result_list=NULL,draw=FALSE,main=''){
  each_value_CV <- unlist(lapply(result_list$all_result,function(x){x$AUROC_cv}))
  each_value_MG <- unlist(lapply(result_list$all_result,function(x){x$AUROC}))
  each_feature_num <- unlist(lapply(result_list$all_result,function(x){length(x$use_feature)}))
  ww <- max(which.max(each_value_CV))
  #
  if(draw==TRUE){
    layout(t(matrix(1:2)));par(mar=c(5,4,4,3))
    draw_each <- function(use_value,use_feature_num,main=main,ww=NULL,ylab=''){
      cc1 <- brewer.pal(8,'Set1')[c(1,5,3,4)];
      plot(y=use_value,x= -use_feature_num,pch=16,col=adjustcolor(cc1[1],0.6),
           xlab='Remained feature number',ylab=ylab,xaxt='n',
           main=main,cex=0.5)
      lines(y=use_value,x= -use_feature_num,col=adjustcolor('dark grey',0.8));
      points(y=use_value,x= -use_feature_num,col=adjustcolor(cc1[1],0.6),cex=0.5)
      axis(side=1,at= -seq(min(use_feature_num),max(use_feature_num),length.out=5),
           labels=round(seq(min(use_feature_num),max(use_feature_num),length.out=5)))
      if(is.null(ww)==TRUE) ww <- max(which.max(use_value))
      abline(h=use_value[ww],lty=2,col='purple')
      text(par()$usr[2],use_value[ww],signif(use_value[ww],3),pos=4,xpd=T)
      abline(v=-use_feature_num[ww],lty=2,col='purple')
      text(y=par()$usr[3],x=use_feature_num[ww],use_feature_num[ww],pos=3,xpd=T)
    }
    draw_each(each_value_CV,each_feature_num,main=sprintf('%s\nCV results',main),
              ww=ww)
    draw_each(each_value_MG,each_feature_num,main=sprintf('%s\nfinal model',main),
              ww=ww)    
  }
  use_fit <- result_list$all_result[[ww]]
  return_res <- list()
  return_res$best_fit <- use_fit
  return_res$class_transform <- result_list$class_transform
  return_res$method_name <- result_list$method_name
  return_res$Y_class <- result_list$Y_class
  return(return_res)
}

## max_T: maximum iteration time
generate_mockData <- function(all_fitModel,input_dataset,rand_seed=329,max_T=100,verbose=TRUE){
  set.seed(rand_seed)
  # process final model
  all_use_fit_final <- all_fitModel
  model_feature <- names(all_use_fit_final)
  model_feature_class <- unlist(lapply(all_use_fit_final,function(x)x$Y_class))
  model_feature_transform <- lapply(all_use_fit_final,function(x)x$class_transform)
  model_feature_fit <- lapply(all_use_fit_final,function(x)x$best_fit$fit)
  model_feature_useFeat <- lapply(all_use_fit_final,function(x)x$best_fit$use_feature)
  names(model_feature_class) <- names(model_feature_transform) <- 
    names(model_feature_fit) <- names(model_feature_useFeat) <- model_feature
  wn <- names(model_feature_class)[which(model_feature_class=='numeric')]
  num_model_sd <- unlist(lapply(model_feature_fit[wn],function(x){sd(predict(x))}))
  # "best_fit"                "class_transform" "method_name"     "Y_class" 
  # "AUROC"              "AUROC_cv"           "fit"      "feature_importance" "use_feature"  
  d3 <- input_dataset[,model_feature,drop=F] ## original
  K <- length(model_feature); ## feature number
  w1 <- sample(1:nrow(d3),size=max_T,replace = T);
  old_X <- d3[w1,,drop=F];new_X <- old_X;
  # 随机初始化t=0 时的一组向量 x0=(x0[1],…,x0[K])
  w1 <- sample(1:nrow(d3),size=1,replace = T);X0 <- d3[w1,,drop=F];
  for(each_t in 1:max_T){
    for(i in 1:K){
      pred_feature <- model_feature[i]
      # get the newest feature
      if(each_t == 1 & i == 1){
        new_pX <- X0
      }
      if(each_t == 1 & i>1){
        new_pX <- X0; new_pX[1:(i-1)] <- new_X[each_t,c(1:(i-1))]
      }
      if(each_t > 1 & i == 1){
        new_pX <- new_X[each_t-1,]
      }
      if(each_t > 1 & i > 1){
        new_pX <- new_X[each_t-1,]; new_pX[1:(i-1)] <- new_X[each_t,c(1:(i-1))]
      }
      # predict for X[each_t,i]
      pred_fit <- model_feature_fit[[pred_feature]]
      c1 <- model_feature_class[pred_feature]
      if(c1 == 'numeric'){
        res <- predict(pred_fit,newdata = new_pX)
        # 连续变量是预测值为均值，模型误差为方差，取一个 rnorm
        rs <- rnorm(n=1,mean=res,sd=num_model_sd[pred_feature])
      }else{
        res <- predict(pred_fit,newdata = new_pX,type='prob')
        # binary/category: use probability to select one res and transfer to class
        t1 <- model_feature_transform[[pred_feature]][names(res)]
        rs <- sample(t1,size=1,prob=res)
      }
      new_X[each_t,i] <- rs;
      if(verbose==TRUE) message(sprintf('Time:%s for feature:%s, result is %s',each_t,pred_feature,rs))
    }
  }
  return(new_X)
}
get_sta <- function(TP,TN,FP,FN){
  ACC <- (TP + TN) / (TP + TN + FP + FN)
  TPR <- TP / (TP + FN)
  TNR <- TN / (TN + FP)
  FPR <- FP / (FP + TN)
  PPV <- TP / (TP + FP) # precision
  NPV <- TN / (TN + FN) # recall
  F1_score <- 2 * PPV * TPR / (PPV + TPR)
  c(ACC,TPR,TNR,FPR,PPV,NPV,F1_score)
}
