source('D:/analysis_eng/SherryPlot/functions_mutationGene.R')
## options
use_gene <- 'CYP27A1'
use_uniprot <- 'Q02318'
use_ensg <- 'ENSG00000135929'
use_enst <- 'ENST00000258415'
inner_size <- 100
main_table_file <- 'data/CYP27A1-ZP20201110.xlsx'
add_region_file <- 'data/CYP27A1-ligand.xlsx'

##
use_gene <- 'SMARCA4'
use_uniprot <- 'Q02318'
use_ensg <- 'ENSG00000135929'
use_enst <- 'ENST00000450717'
inner_size <- 0
main_table_file <- 'data/SMARCA4-210318.xlsx'
add_region_file <- 'data/SMARCA4-ligand.xlsx'

######################################################## input
variant_info <- read.xlsx(main_table_file)
## add interest site
tmp1 <- read.xlsx(add_region_file)
use_pf <- do.call(rbind,lapply(colnames(tmp1),function(i){
  x <- tmp1[,i]
  x <- x[which(is.na(x)==F)]
  x2 <- do.call(rbind,lapply(x,function(x1){
    if(grepl('-',x1)){
      x2 <- as.numeric(unlist(strsplit(x1,'-')))
      x2 <- data.frame(Source='Curate',Domain=i,Start=min(x2),End=max(x2),
                       stringsAsFactors = F)
    }else{
      x2 <- data.frame(Source='Curate',Domain=i,Start=x1,End=x1,
                       stringsAsFactors = F)  
    }
  }))
  x2
}))
################################
exon_tab <- get_exon(use_enst)
cds_tab <- get_cds(use_enst)
transcript_tab <- get_transcript(use_enst)
use_cds_tab <- get_use_cdstab(transcript_tab,cds_tab)
### draw
# mod_variant_info() could modify column
## deal with splicing site
if(inner_size>0){
  variant_info$CDS_pos <- Pos2RNA(pos=variant_info$CDS_pos,
                                  cds_start=use_cds_tab$start,
                                  cds_end=use_cds_tab$end,
                                  first_pos=0,inner_size=inner_size)
  w1 <- grep("\\+|\\-",variant_info$Variant)
  w2 <- as.numeric(gsub('c.([0-9]+)([\\+\\-])([0-9]+)[A-Z].*','\\2\\3',variant_info$Variant))[w1]
  variant_info$Variant[w1] <- variant_info$Variant[w1]
  variant_info$CDS_pos[w1] <- variant_info$CDS_pos[w1]+inner_size/2*sign(w2)

}
################################
pdf('plot/mutationGene_f1.pdf',width=8,height=4)
draw_domain_BySample(transcript_tab,exon_tab,cds_tab,use_pf,
  variant_info,mark_outcome='liver failure',unit=0.35,unit_s=0.4,
  inner_size=inner_size)
dev.off()
pdf('plot/mutationGene_f2.pdf',width=10,height=4)
draw_domain_BySite(transcript_tab,exon_tab,cds_tab,use_pf,
                   variant_info,
                   mark_outcome='liver failure',mark_type='variant',
                   unit=0.35,unit_s=0.4,inner_size=inner_size)
dev.off()

