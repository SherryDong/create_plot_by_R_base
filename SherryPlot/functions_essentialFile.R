get_1000g <- function(x){
  x1 <- do.call(rbind,lapply(strsplit(x,'\\|'),as.numeric))
  colnames(x1)<-c('1000g_AF','1000g_Hom')
  as.data.frame(x1)
}
get_ExAC <- function(x){
  x1 <- do.call(rbind,lapply(strsplit(x,'\\|'),as.numeric))
  colnames(x1)<-c('ExAC_AC','ExAC_Hom','ExAC_AF')
  as.data.frame(x1)
}
get_gnomAD <- function(x){
  x1 <- do.call(rbind,lapply(strsplit(x,'\\|'),as.numeric))
  colnames(x1)<-c('gnomAD_AC','gnomAD_Hom')
  as.data.frame(x1)
}
get_ACAF_public <- function(x_1000g,x_ExAC,x_gnomAD){
  x1 <- cbind(get_1000g(x_1000g),get_ExAC(x_ExAC),get_gnomAD(x_gnomAD))
  x1$PublicHom <- x1$`1000g_Hom`+x1$ExAC_Hom+x1$gnomAD_Hom
  x1$PublicMaxAF <- apply(x1[,c('1000g_AF','ExAC_AF')],1,max)
  x1
}
get_mutationType <- function(x_annovar,x_vep){
  d_mt <- read.delim('D:/analysis_eng/SherryPlot/data/mutation_type.txt')
  d_mt$VEP.annotated.variant <- gsub(' ','_',d_mt$VEP.annotated.variant)
  u1 <- unique(x_annovar)
  print(setdiff(u1,d_mt$ANNOVAR.annotated.variant))
  u2 <- unique(x_vep)
  print(setdiff(u2,d_mt$VEP.annotated.variant))
  d_mt$ID <- sprintf('%s-%s',d_mt$ANNOVAR.annotated.variant,d_mt$VEP.annotated.variant)
  rownames(d_mt) <- d_mt$ID
  x_both <- sprintf('%s-%s',x_annovar,x_vep)
  x_both2impact <- d_mt[x_both,c('Impact')]
  x_both2group <- d_mt[x_both,c('Mutation.Type.Group')]
  return(list(impact=x_both2impact,group=x_both2group))
}
# Hom | Het | Lht | Hemi
# samples and families
# cannot work for chrX/Y
get_ACAN_InternalComments <- function(x){
  if(x=='.'){
    return(rep(NA,8))
  }
  x1 <- strsplit(x,';')[[1]]
  x2 <- as.numeric(strsplit(x1[1],":|\\||\\)")[[1]])
  x2 <- x2[which(is.na(x2)==F)]
  x3 <- as.numeric(strsplit(x1[2],' |\\(')[[1]])
  x3 <- x3[which(is.na(x3)==F)]
  AC <- x2[1]*2+x2[2]*1+x2[3]*1+x2[4]*1
  AN <- x3[1]*2 # not for chrX
  names(x2) <- c('Internal_Hom','Internal_Het','Internal_Lht','Internal_Hemi')
  r1 <- c(x2,Internal_total_sample=x3[1],Internal_AC=AC,Internal_AN=AN,Internal_AF=AC/AN)
  return(r1)
}
get_ACAN_InternalComments_list <- function(x,tag=''){
  x1 <- as.data.frame(do.call(rbind,lapply(x,get_ACAN_InternalComments)))
  if(tag!='') colnames(x1) <- sprintf('%s_%s',tag,colnames(x1))
  x1
}

burden_filter <- function(input_score_mat,
                          min_public_hom=1,
                          min_public_AF=0.01,
                          min_local_hom=3,
                          min_local_AF=0.01,
                          min_CADD=25,
                          min_REVEL=0.1){
  m1 <- input_score_mat
  ##
  m1 <- m1[which(m1$FILTER=='PASS'),]
  cn <- colnames(m1)
  
  # filter mutation type: Mutation_type_(VEP) Mutation_type_(annovar) --> only missense/PTR
  if('Mutation_type_(annovar)' %in% cn){
    x_annovar <- m1$`Mutation_type_(annovar)`  
  }else{
    x_annovar <- m1$Mutation_type_.annovar.
  }
  if('Mutation_type_(VEP)' %in% cn){
    x_vep <- m1$`Mutation_type_(VEP)`  
  }else{
    x_vep <- m1$Mutation_type_.VEP.
  }
  
  r1 <- get_mutationType(x_annovar,x_vep)
  m1$CADD <- as.numeric(m1$CADD)
  m1$REVEL <- as.numeric(m1$REVEL)
  if(!'ExAC(AC|AF)' %in% cn) m1$`ExAC(AC|AF)` <- m1$ExAC.AC.AF.
  if(!'GNOMAD(AC|Hom)' %in% cn) m1$`GNOMAD(AC|Hom)` <- m1$GNOMAD.AC.Hom.
  if(!'1000_genome' %in% cn) m1$`1000_genome` <- m1$X1000_genome
  m1_control <- m1[which(r1$impact %in% c('Low')),] # use as control
  m1_case_PTV <- m1[which(r1$group=='PTV'),] # use as case PTV
  m1_case_MIS <- m1[which(r1$group=='MIS'),] # use as case MIS
  all_m1 <- list(PTV=m1_case_PTV,MIS=m1_case_MIS,control=m1_control)
  ##
  # filter AC/AF public: "1000_genome" "ExAC(AC|AF)" GNOMAD(AC|AF)
  r_pub <- lapply(all_m1,function(x){
    get_ACAF_public(x$`1000_genome`,x$`ExAC(AC|AF)`,x$`GNOMAD(AC|Hom)`)
  })
  # filter AC/AF internal: InternalComments* 
  r_internal <- lapply(all_m1,function(x){
    get_ACAN_InternalComments_list(x$InternalComments,'WES')
  })
  r_internal_p <- lapply(all_m1,function(x){
    get_ACAN_InternalComments_list(x$InternalComments_Panel,'panel')
  })
  r_internal_w <- lapply(all_m1,function(x){
    get_ACAN_InternalComments_list(x$InternalComments_WGS,'WGS')
  })
  #
  r_Hom <- lapply(names(all_m1),function(x){
    apply(cbind(r_internal[[x]]$WES_Internal_Hom,r_internal_p[[x]]$panel_Internal_Hom,r_internal_w[[x]]$WGS_Internal_Hom),1,function(x)sum(x,na.rm=T))
  })
  r_Hemi <- lapply(names(all_m1),function(x){
    apply(cbind(r_internal[[x]]$WES_Internal_Hemi,r_internal_p[[x]]$panel_Internal_Hemi,r_internal_w[[x]]$WGS_Internal_Hemi),1,function(x)sum(x,na.rm=T))
  })
  r_AF <- lapply(names(all_m1),function(x){
    apply(cbind(r_internal[[x]]$WES_Internal_AF,r_internal_p[[x]]$panel_Internal_AF,r_internal_w[[x]]$WGS_Internal_AF),1,function(x)max(x,na.rm=T))
  })
  names(r_Hom) <- names(r_Hemi) <- names(r_AF) <- names(all_m1)
  r_Hom <- lapply(names(all_m1),function(x){
    chrX_PTV <- which(all_m1[[x]]$`#CHROM`=='X')
    r_Hom[[x]][chrX_PTV] <- r_Hom[[x]][chrX_PTV]+r_Hemi[[x]][chrX_PTV]  
    r_Hom[[x]]
  })
  names(r_Hom) <- names(r_Hemi) <- names(r_AF) <- names(all_m1)
  r_combine <- lapply(names(all_m1),function(x){
    cbind(r_pub[[x]],r_internal[[x]],r_internal_p[[x]],r_internal_w[[x]],
          InternalHom=r_Hom[[x]],InternalMaxAF=r_AF[[x]])
  })
  names(r_combine) <- names(all_m1)
  
  ##################################
  #   all_m1 <- list(PTV=m1_case_PTV,MIS=m1_case_MIS,control=m1_control)
  all_m1 <- lapply(all_m1,function(x){x$filter_tag<-'';x})
  n1 <- unlist(lapply(all_m1,function(x)nrow(x[which(x$filter_tag==''),])))
  message(sprintf('%s for PTV, %s for MIS and %s for control',n1[1],n1[2],n1[3]))
  
  # Hom=0 AF<=0.01
  for(i in names(all_m1)){
    x <- all_m1[[i]]
    r <- r_combine[[i]]
    w1 <- which(r$PublicHom>min_public_hom);
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('PublicHomCount:%s',r$PublicHom[w1]))
    w1 <- which(r$PublicMaxAF>min_public_AF);
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('PublicAF:%.5f',r$PublicMaxAF[w1]))
    all_m1[[i]] <- x
  }
  n1 <- unlist(lapply(all_m1,function(x)nrow(x[which(x$filter_tag==''),])))
  message(sprintf('%s for PTV, %s for MIS and %s for control',n1[1],n1[2],n1[3]))
  
  # filter by local Hom/AF
  for(i in names(all_m1)){
    x <- all_m1[[i]]
    r <- r_combine[[i]]
    w1 <- which(r$InternalHom>min_local_hom);
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('LocalHomCount:%s',r$InternalHom[w1]))
    w1 <- which(r$InternalMaxAF>min_local_AF);
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('LocalAF:%.5f',r$InternalMaxAF[w1]))
    all_m1[[i]] <- x
  }
  n1 <- unlist(lapply(all_m1,function(x)nrow(x[which(x$filter_tag==''),])))
  message(sprintf('%s for PTV, %s for MIS and %s for control',n1[1],n1[2],n1[3]))
  
  # filter function: nsdb_SIFT_score	nsdb_SIFT_prediction	nsdb_Polyphen2_score	nsdb_polyphen2_prediction	nsdb_MutationTaster_score	nsdb_MutationTaster_prediction CADD	SPIDEX	dbscSNV_RF_SCORE	REVEL
  for(i in names(all_m1)){
    x <- all_m1[[i]]
    r <- r_combine[[i]]
    w1 <- which(x$CADD<min_CADD & x$REVEL<min_REVEL | x$CADD<min_CADD & is.na(x$REVEL)==TRUE | is.na(x$CADD)==TRUE & x$REVEL<min_REVEL)
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('CADD:%s;REVEL:%s',x$CADD[w1],x$REVEL[w1]))
    all_m1[[i]] <- x
  }
  n1 <- unlist(lapply(all_m1,function(x)nrow(x[which(x$filter_tag==''),])))
  message(sprintf('%s for PTV, %s for MIS and %s for control',n1[1],n1[2],n1[3]))
  
  all_m1 <- lapply(names(all_m1),function(x){
    cbind(all_m1[[x]],r_combine[[x]])
  })
  names(all_m1) <- names(r_combine)
  return(all_m1)
}

burden_filter_2 <- function(input_list,Lht_thre=0.2,DP_thre=10,
                            min_CADD=25,
                            min_REVEL=0.1,
                            min_REFALT_len=2){
  all_m1 <- input_list
  for(i in names(all_m1)){
    x <- all_m1[[i]]
    x1 <- as.numeric(gsub('(.*):(.*),(.*):(.*)','\\3',x$FORMAT))
    x2 <- as.numeric(gsub('(.*):(.*),(.*):(.*)','\\4',x$FORMAT))
    x12 <- x1/x2
    w1 <- which(x12<=Lht_thre)
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('Lht:%.5f',x12[w1]))
    w1 <- which(x2<=DP_thre)
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('DP:%s',x2[w1]))
    all_m1[[i]] <- x
  }
  n1 <- unlist(lapply(all_m1,function(x)nrow(x[which(x$filter_tag==''),])))
  message(sprintf('%s for PTV, %s for MIS and %s for control',n1[1],n1[2],n1[3]))
  
  ##nsdb_SIFT_prediction/nsdb_polyphen2_prediction/nsdb_MutationTaster_prediction
  
  for(i in names(all_m1)){
    x <- all_m1[[i]]
    w1 <- which(x$nsdb_SIFT_prediction %in% c('B','P','T','A','N'))
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('SIFT:%s',x$nsdb_SIFT_prediction[w1]))
    w1 <- which(x$nsdb_polyphen2_prediction %in% c('B','P','T','A','N'))
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('polyphen2:%s',x$nsdb_polyphen2_prediction[w1]))
    w1 <- which(x$nsdb_MutationTaster_prediction %in% c('B','P','T','A','N'))
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('MutationTaster:%s',x$nsdb_MutationTaster_prediction[w1]))
    all_m1[[i]] <- x
  }
  n1 <- unlist(lapply(all_m1,function(x)nrow(x[which(x$filter_tag==''),])))
  message(sprintf('%s for PTV, %s for MIS and %s for control',n1[1],n1[2],n1[3]))
  
  ## remove complicated indel
  for(i in names(all_m1)){
    x <- all_m1[[i]]
    ref_len <- nchar(x$REF)
    alt_len <- nchar(x$ALT)
    w1 <- which(ref_len>min_REFALT_len | alt_len>min_REFALT_len)
    x$filter_tag <- add_tag(x$filter_tag,w1,sprintf('REF:%s;ALT:%s',ref_len[w1],alt_len[w1]))
    all_m1[[i]] <- x
  }
  n1 <- unlist(lapply(all_m1,function(x)nrow(x[which(x$filter_tag==''),])))
  message(sprintf('%s for PTV, %s for MIS and %s for control',n1[1],n1[2],n1[3]))
  
  return(all_m1)
}
add_tag <- function(x,w,tag){
  x[w] <- sprintf('%s;%s',x[w],tag)
  x <- gsub('^;','',x)
  print(sort(table(x),decreasing = T)[1:10])
  x
}  