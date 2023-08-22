draw.heatmap.2 <- function(mat=NULL,use_genes=rownames(mat),use_gene_label=use_genes,use_samples=colnames(mat),use_sample_label=use_samples,
         phenotype_info=NULL,use_phe=NULL,main="",scale='none',pdf_file=NULL,
         cluster_rows=TRUE,cluster_columns=TRUE,
         show_row_names=TRUE,show_column_names=TRUE,
         clustering_distance_rows='pearson',clustering_distance_columns='pearson',
         use_color=NULL,pre_define=NULL,
         ...){
  #
  all_input_para <- c('mat','use_genes','use_gene_label','use_samples','use_sample_label','main','scale',
                      'cluster_rows','cluster_columns','show_row_names','show_column_names','clustering_distance_rows','clustering_distance_columns')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('cluster_rows',c(TRUE,FALSE),envir=environment()),
                 check_option('cluster_columns',c(TRUE,FALSE),envir=environment()),
                 check_option('show_row_names',c(TRUE,FALSE),envir=environment()),
                 check_option('show_column_names',c(TRUE,FALSE),envir=environment()),
                 check_option('scale',c("none", "row",'column'),envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  names(use_gene_label) <- use_genes
  names(use_sample_label) <- use_samples
  if(is.null(rownames(phenotype_info))==FALSE){
    ori_phenotype_info <- phenotype_info
    phenotype_info <- as.data.frame(phenotype_info[colnames(mat),],stringsAsFactors=FALSE)
    colnames(phenotype_info) <- colnames(ori_phenotype_info)
  }
  for(i in colnames(phenotype_info)){
    phenotype_info[,i] <- clean_charVector(phenotype_info[,i])
  }
  if(exists('row_names_gp')==FALSE) row_names_gp <- gpar(fontsize = 12)
  if(exists('column_names_gp')==FALSE) column_names_gp <- gpar(fontsize = 12)
  use_genes <- base::intersect(use_genes,rownames(mat))
  use_samples <- base::intersect(use_samples,colnames(mat))
  use_mat <- mat[use_genes,use_samples]
  rownames(use_mat) <- use_gene_label[rownames(use_mat)]
  colnames(use_mat) <- use_sample_label[colnames(use_mat)]
  row_names_max_width <- base::max(strwidthMod(rownames(use_mat),'inches',cex=row_names_gp[[1]]/7))
  row_names_max_width <- unit(row_names_max_width,'inches')
  column_names_max_height <- base::max(strwidthMod(colnames(use_mat),'inches',cex=column_names_gp[[1]]/7))
  column_names_max_height <- unit(column_names_max_height,'inches')
  if(scale=='row'){use_mat <- t(apply(use_mat,1,do.std))}
  if(scale=='column'){use_mat <- apply(use_mat,2,do.std)}
  if(base::length(use_phe)==0){
    if(scale=='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main,name='Raw value',
                                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                                     show_row_names=show_row_names,show_column_names=show_column_names,
                                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,
                                     ...)
    }
    if(scale!='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main, name='Z value',
                                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                                     show_row_names=show_row_names,show_column_names=show_column_names,
                                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
  }else{
    if(base::length(use_phe)==1){
      use_phe_info <- as.data.frame(phenotype_info[,use_phe],stringsAsFactors=FALSE)
      rownames(use_phe_info) <- rownames(phenotype_info)
      colnames(use_phe_info) <- gsub(' ','.',use_phe)
    }else{
      use_phe_info <- phenotype_info[,use_phe]
      colnames(use_phe_info) <- gsub(' ','.',use_phe)
    }
    use_phe <- colnames(use_phe_info)
    l2c <- get.class.color(base::unique(as.character(as.matrix(use_phe_info))),use_color=use_color,pre_define=pre_define)
    use_col <- lapply(use_phe,function(x)l2c[base::unique(use_phe_info[,x])])
    names(use_col) <- use_phe
    ha_column <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(use_phe_info),col = use_col)
    if(scale=='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main, top_annotation = ha_column,name='Raw value',
                                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                                     show_row_names=show_row_names,show_column_names=show_column_names,
                                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
    if(scale!='none'){
      ht1 <- ComplexHeatmap::Heatmap(use_mat, column_title = main, top_annotation = ha_column,name='Z value',
                                     cluster_rows=cluster_rows,cluster_columns=cluster_columns,
                                     clustering_distance_rows=clustering_distance_rows,clustering_distance_columns=clustering_distance_columns,
                                     show_row_names=show_row_names,show_column_names=show_column_names,
                                     row_names_max_width=row_names_max_width,column_names_max_height=column_names_max_height,...)
    }
  }
  ht_list <- ht1
  if(is.null(pdf_file)==FALSE){
    ww <- 1.25*column_names_gp[[1]]/72*ncol(use_mat)+base::max(strwidthMod(rownames(use_mat),'inches',ori=TRUE))+5
    hh <- 1.25*row_names_gp[[1]]/72*nrow(use_mat)+base::max(strwidthMod(colnames(use_mat),'inches',ori=TRUE))+3
    pdf(pdf_file,width=ww,height=hh)
  }
  ComplexHeatmap::draw(ht_list,heatmap_legend_side='left',annotation_legend_side='right')
  if(is.null(pdf_file)==FALSE) {while (!is.null(dev.list()))  dev.off();}
  return(TRUE)
}