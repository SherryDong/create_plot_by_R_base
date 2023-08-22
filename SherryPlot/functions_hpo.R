## functions for hpo
library(igraph)
hp_p2c <- read.delim('D:/analysis_eng/SherryPlot/data/hp_p2c.txt',stringsAsFactors = F,header = F)
hp_info <- read.delim('D:/analysis_eng/SherryPlot/data/hp_info.txt',stringsAsFactors = F,header = F)
rownames(hp_info) <- hp_info$V1
hp_net <- graph.data.frame(hp_p2c,directed = T)
##
find_nearest2root <- function(hp_root='HP:0000118',all_hp){
  hp_net_dist <- distances(hp_net,v=all_hp,to=hp_root,mode='in')
  all_d <- rowSums(hp_net_dist)
  return(hp_net_dist[order(all_d),,drop=F])
}
get_hpo2root <- function(hp_root=hp_p2c[which(hp_p2c$V1=='HP:0000118'),2]){
  hp_net_dist <- distances(hp_net,to=hp_root,mode='in')
  hpo2root_uni <- do.call(rbind,lapply(hp_root,function(x){
    cbind(names(which(hp_net_dist[,x]!=Inf)),x)
  }))
  hpo2root_uni <- as.data.frame(hpo2root_uni,stringsAsFactors=F)
  colnames(hpo2root_uni) <- c('HPO.ID','root.HPO')  
  return(hpo2root_uni)
}
impute_hpo2root <- function(x,hp_root='HP:0000118'){
  x1 <- all_simple_paths(hp_net,from=x,to=hp_root,mode='in')
  unique(unlist(lapply(x1,function(xx)names(xx))))
}
impute_hpo2root_list <- function(all_x,hp_root='HP:0000118'){
  all_x <- intersect(all_x,V(hp_net)$name)
  res <- c()
  for(x in all_x){
    res <- c(res,impute_hpo2root(x)) 
  }
  unique(res)
}
load_hpoDistbg <- function(RData_path='CCGT',sample2hpo_path='D:/写写文章/GTLC/GTLC/CCGT_result_202103/data/sample2hpo.txt'){
  if(file.exists(RData_path)){
    load(RData_path)
  }else{
    sample2hpo_bg <- read.delim(sample2hpo_path,stringsAsFactors = F,header=F)
    sample2hpo_bg_list <- lapply(unique(sample2hpo_bg$V1),function(x)unique(sample2hpo_bg[which(sample2hpo_bg$V1==x),2]))
    names(sample2hpo_bg_list) <- unique(sample2hpo_bg$V1)
    sample2hpo_bg_list <- lapply(sample2hpo_bg_list,function(x){
      impute_hpo2root_list(x)
    })
    w1 <- unlist(lapply(sample2hpo_bg_list,length));u1 <- which(w1>0);N_bg <- length(u1)
    sample2hpo_bg_list <- sample2hpo_bg_list[u1]
    hpDist_count <- sort(table(unlist(sample2hpo_bg_list)),decreasing = T);
    hpDist <- sort(hpDist_count/N_bg,decreasing = T)
    save(sample2hpo_bg_list,N_bg,hpDist_count,hpDist,file=RData_path)
  }
}
find_sampleWithHPO <- function(x='',sample2hpo){
  w1 <- lapply(sample2hpo,function(xx)intersect(xx,x))
  w2 <- unlist(lapply(w1,length))
  names(which(w2>0))
}
hpo2case_from_case2hpo <- function(case2hp_list2root){
  hp2case <- do.call(rbind,lapply(names(case2hp_list2root),function(x){cbind(x,case2hp_list2root[[x]])}))
  all_hp <- unique(hp2case[,2])
  hp2case <- lapply(all_hp,function(x){
    hp2case[which(hp2case[,2]==x),1]
  })
  names(hp2case) <- all_hp;hp2case
}
####################################### HPO enrichment analysis
preprocess=0
if(preprocess==1){
  load('D:/analysis_eng/SherryPlot/data/OMIM/phenotype2gene.RData')
  gene2hpo_list <- lapply(unique(phenotype2gene$GeneSymbol),function(x){
    unique(phenotype2gene$HPID[which(phenotype2gene$GeneSymbol==x)])
  })
  names(gene2hpo_list) <- unique(phenotype2gene$GeneSymbol)
  source('D:/analysis_eng/SherryPlot/functions_hpo.R')
  gene2hpo_list_full <- lapply(gene2hpo_list,function(x){
    impute_hpo2root_list(x)
  })
  tmp1 <- list2mat(gene2hpo_list_full)
  gr <- graph_from_incidence_matrix(tmp1)
  tmp2 <- igraph::as_data_frame(gr, what="edges")
  tmp3 <- aggregate(tmp2$to,list(tmp2$from),unique)
  hpo2gene <- tmp3[,-1]
  names(hpo2gene) <- sprintf('%s_%s',tmp3[,1],hp_info[tmp3[,1],2])
  save(hpo2gene,file='D:/analysis_eng/SherryPlot/data/OMIM/hpo2gene.RData')
  gene2hpo_list_full<-lapply(gene2hpo_list_full,function(x)sprintf('%s_%s',x,hp_info[x,2]))
  gene2hpo_list <- lapply(gene2hpo_list,function(x)sprintf('%s_%s',x,hp_info[x,2]))
  save(gene2hpo_list_full,file='D:/analysis_eng/SherryPlot/data/OMIM/gene2hpo_list_full.RData')
  save(gene2hpo_list,file='D:/analysis_eng/SherryPlot/data/OMIM/gene2hpo_list.RData')
}
load('D:/analysis_eng/SherryPlot/data/OMIM/hpo2gene.RData')

