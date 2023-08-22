library(deconvSeq)
library(DESeq2)
library(NetBID2)
library(scMINER)
#dependencies
library(e1071)
library(parallel)
#source('/home/dongxinran/project/SingleCell/soft/preprocessCore/R/normalize.quantiles.R')
library(preprocessCore)
source('D:/analysis_eng/SherryPlot//CIBERSORT.R')
## demo run
# input: exp_mat, cell_type, group
# https://rosedu1.github.io/deconvSeq_vignette.html
# devtools::install_github("rosedu1/deconvSeq", dependencies=TRUE)
get_signatureMatrix <- function(exp_mat,cell_type){
  use_gene  <- rownames(exp_mat)
  names(cell_type) <- colnames(exp_mat)
  # filter genes and samples
  IQR1 <- aggregate(t(exp_mat),list(cell_type),sd)
  w1 <- unlist(lapply(IQR1[,-1],min));w1 <- which(w1>0)
  IQR2 <- apply(exp_mat,2,sd);w2 <- which(IQR2>0)
  exp_mat <- exp_mat[w1,w2]; cell_type <- cell_type[colnames(exp_mat)]
  ## signature matrix
  group <- factor(cell_type);design <- model.matrix(~0 + group)
  colnames(design) <- gsub(' ','_',levels(group));
  rownames(design) <- colnames(exp_mat)
  exp_mat <- as.matrix(exp_mat)
  design.singlecell <- design 
  dge.celltypes = getdge(exp_mat, design.singlecell,ncpm.min=1, nsamp.min=4, method="bin.loess") ## input original count
  b0 = getb0.rnaseq(dge.celltypes, design.singlecell, ncpm.min=1, nsamp.min=4) # model
  X <- b0$b0; X <- X[which(is.na(X[,1])==F),]
  return(X)
}
deconvolution.Activity <- function(target_mat,ref_mat,ref_celltype){
  #假设我们有一个RNASeq数据矩阵X，每行代表一个基因，每列代表一个样本。我们还有一个参考基因表达矩阵G，每行代表一个基因，每列代表一个细胞类型。我们的目标是从X中估计不同细胞类型的相对丰度。
  #假设X可以表示为G的线性组合，即X = G*S + E，其中S是细胞类型相对丰度矩阵，E是误差矩阵。我们可以使用线性回归模型来拟合这个方程，从而估计S
  # X = G*S + E
  # S = (G'G)^(-1)G'X
  tmp1 <- aggregate(t(ref_mat),list(ref_celltype),median)
  G <- as.matrix(t(tmp1[,-1]));colnames(G) <- tmp1$Group.1
  w1 <- intersect(rownames(G),rownames(target_mat))
  G <- G[w1,];IQR <- IQR.filter(G,thre=0.95); G <- G[IQR,]
  X <- target_mat[rownames(G),]
  S <- solve(t(G) %*% G) %*% t(G) %*% X
  S_norm <- apply(S,2,function(x)(x-mean(x))/sd(x))
  r1 <- sigmoid(S_norm)/colSums(sigmoid(S_norm))
  return(r1)
}

#### demo for GSE104276 human PFC, expression
load('D:/写写文章/GD-Driver/BrainCortexDriver_project/BrainCortexDriver/DATA/GSE104276_scminerResult.RData')
eset.log2 <- output_result$eset.log2
exp_mat <- as.matrix(Biobase::exprs(eset.log2))
cell_type <- pData(eset.log2)$cell_groups
X <- get_signatureMatrix(exp_mat,cell_type)
save(X,exp_mat,cell_type,file='D:/analysis_eng/SherryPlot/data/CIBERSORT_humanPFC.RData')
####

#### demo for GSE162170 human cerebral cortex, expression
load('D:/写写文章/GD-Driver/BrainCortexDriver_project/BrainCortexDriver/DATA/GSE162170_scminerResult.RData')
eset.log2 <- output_result$eset.log2
exp_mat <- as.matrix(Biobase::exprs(eset.log2))
cell_type <- pData(eset.log2)$celltype
X <- get_signatureMatrix(exp_mat,cell_type)
save(X,exp_mat,cell_type,file='D:/analysis_eng/SherryPlot/data/CIBERSORT_GSE162170.RData')

#### demo for GSE104276 human PFC, activity
load('D:/写写文章/GD-Driver/BrainCortexDriver_project/BrainCortexDriver/DATA/GSE104276_scminerResult.RData')
eset.log2 <- output_result$eset.log2
exp_mat <- as.matrix(Biobase::exprs(eset.log2))
cell_type <- pData(eset.log2)$cell_groups
load('D:/analysis_eng/SherryPlot/data/brainspan_net_whole.RData') # network
use_ac <- cal.Activity(igraph_obj = brainspan_net_whole$igraph_obj,
                       cal_mat = exp_mat,
                       es.method = 'weightedmean')
use_ac_exp <- exp(use_ac)
X <- get_signatureMatrix(use_ac_exp,cell_type)
save(X,use_ac_exp,cell_type,file='D:/analysis_eng/SherryPlot/data/CIBERSORT_humanPFC_activityBrainSpan.RData')
####


### demo for SYTL3 project with hESC/hNPC/Neuron
#exp_mat <- read.delim('D:/工作工作/绘制图表/XiongMan/RNASeq/ref/gene_count_mat.xls')
#genes <- gsub('(.*)_(.*)','\\2',exp_mat$EnsemblGene_GeneSymbol)
#exp_mat <- apply(exp_mat[,-1],2,as.numeric);rownames(exp_mat) <- genes
#cell_type <- read.xlsx('D:/工作工作/绘制图表/XiongMan/RNASeq/SYTL3_RNASeq_phenotype.xlsx')$design
#w1 <- grep('WT',cell_type)
#exp_mat <- exp_mat[,w1]; cell_type <- cell_type[w1]
#X <- get_signatureMatrix(exp_mat,cell_type)
#save(X,file='D:/analysis_eng/SherryPlot/data/CIBERSORT_humanES-NPC-Neuron.RData')

#### demo for usage
#load('D:/analysis_eng/SherryPlot/data/CIBERSORT_mouseCortex.RData')
#load('data/GSE132044/eset_symbol.RData') ## yechang autism rnaseq (mouse)
#new_mat_symbol <- exprs(eset_symbol)
load('D:/工作工作/绘制图表/XiongMan/ALMS1/ALMS1_eset.RData') ## xiongman alms1 rnaseq (human)
new_mat_symbol <- exprs(use_eset)
new_mat_symbol <- new_mat_symbol[intersect(rownames(new_mat_symbol),rownames(X)),]
######### original
res1 <- CIBERSORT.mod(X = X,Y=new_mat_symbol)
res2 <- t(res1[colnames(new_mat_symbol),colnames(X)])
NetBID2::draw.heatmap(res2,phenotype_info = pData(use_eset),use_phe = 'group')



