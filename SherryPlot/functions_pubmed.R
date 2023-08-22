library(easyPubMed)
library(XML)
##
get_gene_list <- function(query_string,gene_filter=hgmd_gene){
  tmp1 <- get_pubmed_ids(query_string)
  tmp2 <- as.numeric(fetch_all_pubmed_ids(tmp1))
  tmp3 <- pubmed2gene[which(pubmed2gene$PubMed_ID %in% tmp2),] 
  tmp4 <- merge(tmp3,gene_info_1)
  res <- unique(tmp4[,c('PubMed_ID','GeneID','Symbol')])
  res <- res[which(res$Symbol %in% gene_filter),]
  tmp1 <- aggregate(res,list(res$PubMed_ID),function(x)paste(sort(unique(x)),collapse=';'))
  res <- tmp1[,-1]
}
get_gene_list_fromWeb <- function(main_dir,input_files,gene_filter=hgmd_gene){
  # 'pubmed_query/GOF_LOF_20201207/GOF_2000-2015.txt'
  res1 <- lapply(input_files,function(x){
    x1 <- read.delim(sprintf('%s/%s',main_dir,x),stringsAsFactors = F)[[1]]
    x1 <- x1[grep('[0-9]+{8}',x1)]
    x2 <- lapply(x1,function(xx){
      gsub('.*([0-9]+{8}).*','\\1',xx)
    })
    unlist(x2)
  })
  tmp2 <- as.numeric(unique(unlist(res1)))
  tmp3 <- pubmed2gene[which(pubmed2gene$PubMed_ID %in% tmp2),] 
  tmp4 <- merge(tmp3,gene_info_1)
  res <- unique(tmp4[,c('PubMed_ID','GeneID','Symbol')])
  res <- res[which(res$Symbol %in% gene_filter),]
  tmp1 <- aggregate(res,list(res$PubMed_ID),function(x)paste(sort(unique(x)),collapse=';'))
  res <- tmp1[,-1]
  return(res)
}
get_detail_info <- function(pmid){
  print(pmid)
  con <- curl::curl(sprintf("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=%s",pmid))
  x <- readLines(con, n = 1000)
  x1 <- strsplit(x,'\\t')
  title <- x1[[1]]
  abs <- x1[[2]]
  detail <- do.call(rbind,x1[3:length(x1)])
  if(ncol(detail)==0) return(NULL)
  detail_sep <- lapply(unique(detail[,5]),function(x){
    unique(detail[which(detail[,5]==x),4])
  })
  names(detail_sep) <- unique(detail[,5])
  detail_sep$title <- title
  detail_sep$abstract <- abs
  return(detail_sep)
}
get_clinvar_info <- function(clinvarid){
  #clinvarid <- '456552'
  con <- curl::curl(sprintf("https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/",clinvarid))
  x <- readLines(con, n = 10000)
  return(x)
}
#library(easyPubMed) # pubtator
get_pubmed_info <- function(x){
  fd <- get_pubmed_ids(sprintf("%s[ID]",x))
  fr <- fetch_pubmed_data(fd)
  at <- custom_grep(fr, "ArticleTitle", "char")
  ad <- custom_grep(fr, "ArticleDate", "char")
  return(list(title=at,date=ad,info=fr))
}
