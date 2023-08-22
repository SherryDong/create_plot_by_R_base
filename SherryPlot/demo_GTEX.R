library(NetBID2)
#################
load('data/GTEx_dds.RData')
load('data/GTEx_feature_info.RData')
#
dds <- DESeq(dds)
#save(dds,file='data/GTEx_dds_normalize.RData')
vsd <- DESeq2::vst(dds)
mat <- SummarizedExperiment::assay(vsd)
tmp1 <- colData(dds)@listData
tmp1 <- as.data.frame(tmp1)
eset <- generate.eset(exp_mat=mat, 
                      phenotype_info = tmp1, 
                      feature_info = feature_info)
#save(eset,file='data/GTEx_eset_normalize.RData')
eset_symbol <- update_eset.feature(eset,use_feature_info = feature_info,
                                   from_feature = 'Name',
                                   to_feature = 'Description')
save(eset_symbol,file='data/GTEx_esetSymbol_normalize.RData')
#
tmp1 <- exprs(eset_symbol)
tmp2 <- pData(eset_symbol)
tmp3 <- aggregate(t(tmp1),list(tmp2$design),median)
mat_median <- t(tmp3[,-1])
colnames(mat_median) <- tmp3[,1]
save(mat_median,file='data/GTEx_mat_median.RData')

