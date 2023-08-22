## markers
panglodb <- read.delim('D:/写写文章/GD-Driver/BrainCortexDriver_project/data/PanglaoDB_markers_27_Mar_2020.tsv')
panglodb_brain <- panglodb[which(panglodb$organ=='Brain'&panglodb$gene.type=='protein-coding gene'&panglodb$species!='Mm'&panglodb$official.gene.symbol%in%feature.data$geneSymbol),]
genes_of_interest_marker_panglodb <- data.frame(celltype=panglodb_brain$cell.type,
                                       markers=panglodb_brain$official.gene.symbol,
                                       weight=1)
#### from paper PMID34390642
genes_of_interest_list <- list(
  # area
  Cortex=c('FOXG1'),Cerebellum=c('ZIC2'),allocortex=c('NRP1'),
  SubplateMarkers=c('NR4A2','CRYM','NEFL','SERPINI1'), # SP
  SubplateNeurons=c('NR4A2','CRYM','ST18','CDH18'), # 
  MedialGanglionicEminence=c('LHX6','SST'), # MGE (interneurons)
  CaudalGanglionicEminence=c('SP8','NR2F2'), # CGE (interneurons)
  PallialSubpallialBoundary=c('MEIS2','ETV1','PAX6'), # PSB (interneurons)
  # cell type
  RadialGlia=c('SOX9','HES1','SOX2','PAX6','GFAP','VIM','NES','ATP1A2'), # RG
  OPC_Oligodendrocyte=c('SOX10','NKX2-2','MBP'), #OPC_Oligo
  ExcitatoryNeurons=c('NEUROD2','NEUROD6','RBFOX1'), # eN
  MicroGlia=c('AIF1','CCL3','C1QC','CX3CR1','PTPRC'), # MG
  IntermediateProgenitorCell=c('EOMES','PPP1R17','NEUROG1'), # IPC
  InterNeurons=c('DLX2','GAD2','GAD1'), # IN (cortical interneurons)
  EndotheliaCells=c('CLDN5','PECAM1'), # EC
  Astrocyte=c('AQP4','APOE','AGT'), # 
  #
  NeuronalIntermediateProgenitorCell=c('EOMES','PPP1R17','PENK','NEUROG1','NEUROD2'), # nIPC
  Pericytes=c('FOXC2','PDGFRB'), # Peric
  LeptomeningealCells=c('COL1A1','LUM'), # VLMC
  RedBloodCells=c('HEMGN'), # RBC
  CiliatedEpendymalCells=c('FOXJ1'), # 
  Astroglia=c('GFAP','HOPX','EGFR','ASCL1','AQP4'), 
  CorticoGenesis=c('SOX9','EOMES','NEUROD2','DLX2'),
  GlutamatergicNeurons=c('NEUROD2','TBR1','BCL11B','CTIP2','SATB2','SLC17A7','VGLUT1'), # GluN
  CyclingCells=c('TOP2A','MKI67','CLSPN','AURKA'),
  VentricularRadialGlia=c('FBXO32','CTGF','HMGA2'), # vRG
  OuterRadialGlia=c('MOXD1','HOPX','FAM107A','MT3'), # oRG
  EarlyRadialGlia=c('NPY','FGFR3'),  # EarlyRG
  LateRadialGlia=c('CD9','GPX3','TNC'), # LateRG
  TruncatedRadialGlia=c('CRYAB','NR4A1','FOXJ1'), # tRG
  MultipotentGlialProgenitorCells=c('ASCL1','OLIG1','PDGFRA','EGFR'), # mGPC
  GABAergicNeurons=c('DLX2'), #
  ##
  StemCell=c("OLIG2","SOX10","NKX2-2","NKX6-2","PAX7","DBX2","EMX1"),
  OligodendrocyteProgenitorCells=c("OLIG2","OLIG1","SOX10","SOX9","CSPG4","CNP"), # OPC
  MatureOligodendrocyte=c("CNP","CLDN11","MAG","MAL","PLP1","SMARCA4","GEMIN2","CD9","MYT1"),
  Myelin=c("CNP","MBP","GALC","CD9","MAG","MOBP","MOG","MAL","PLP1","MYT1"),
  Vasculogenesis=c("VEGFA"),
  Angiogenesis=c("VEGFA","VEGFB","VEGFC","NRP1","NRP2"),
  NeuronStemstage=c("FOXG1","OTX2","DLX2","PAX6","EMX1","HES5","LHX2","EMX2","NKX2-1","SIX3","GSX2","EOMES","MEIS2","DLX1","ISL1","ASCL1"),
  Neuronalprecursor=c("SOX2","SOX9","NHLH1","EBF2","NEUROG1","NEUROD4","DCX"),
  PostMitoticNeuronalMarker=c("TUBB3","MAP2","MAPT","CUX1","CUX2","SATB2","CDH10","DKK3","TBR1","FEZF2","SST","NPY","PROX1"),
  # GSE104276
  PFC_Microglia=c('PTPRC','P2RY12'),
  PFC_NPCs=c('PAX6','SFPR1'),
  PFC_OPCs=c('OLIG1','PDGFRA','COL20A1','PMP2'),
  PFC_ExcitatoryNeurons=c('NEUROD2','RBFOX1'),
  PFC_InterNeurons=c('GAD1','PDE4DIP'),
  PFC_Astrocytes=c('GFAP','AQP4','SLCO1C1'),
  # Organoid single-cell genomic atlas uncovers human-specific features of brain development
  Organoid_CorticalEN=c('MAP2','ENO2','FOXG1','NEUROD6','NEUROD2','SLC17A7'),
  Organoid_MGECGEIN=c('MAP2','ENO2','FOXG1','DLX5','DLX2','DLX1','GAD1','GAD2','SLC32A1'),
  Organoid_LGEIN=c('MAP2','ENO2','FOXG1','DLX5','DLX2','GAD2','SLC32A1','ISL1','EBF1'),
  Organoid_DiencephalonEN=c('MAP2','ENO2','EBF1','NHLH2','SLC17A6','LHX9','GBX2','SHOX2'),
  Organoid_MesencephalonEN=c('MAP2','ENO2','EBF1','NHLH2','SLC17A6','LHX9','TFAP2B'),
  Organoid_MesencephalonIN=c('MAP2','ENO2','GAD1','GAD2','SLC32A1','GATA3','OTX2','SOX14'),
  # target
  CTDNEP1='CTDNEP1',DYRK1A='DYRK1A',ALMS1='ALMS1')
genes_of_interest <- unique(unlist(genes_of_interest_list))
genes_of_interest_marker <- list2df_narrow(genes_of_interest_list)
colnames(genes_of_interest_marker) <-c('celltype','markers')
genes_of_interest_marker$weight <- 1;
## cell marker from LKY
tmp1 <- read.delim('D:/写写文章/GD-Driver/BrainCortexDriver_project/data/cell_marker.txt',
                   header = F)
tmp2 <- list()
for(i in tmp1$V1){
  if(grepl('\\$',i)){ct=i;}else{tmp2[[gsub('\\$(.*)','\\1',ct)]]<-unique(unlist(strsplit(gsub(".*\\] (.*)","\\1",i),',')))}
}
genes_of_interest_marker_LKY <- list2df_narrow(tmp2)
colnames(genes_of_interest_marker_LKY) <-c('celltype','markers')
genes_of_interest_marker_LKY$weight <- 1;
genes_of_interest_marker_LKY$celltype<-gsub('\\`','',genes_of_interest_marker_LKY$celltype)
## cell marker from CY
tmp1 <- read.xlsx('D:/写写文章/GD-Driver/BrainCortexDriver_project/data/Single cell gene list_CY.xlsx')
tmp2 <- list()
for(i in 1:nrow(tmp1)){
  if(is.na(tmp1[i,1])==F){ct=gsub('^ ','',tmp1[i,1])}else{
    tmp2[[ct]] <- c(tmp2[[ct]],toupper(tmp1[i,2]))
  }
}
genes_of_interest_marker_CY <- list2df_narrow(tmp2)
colnames(genes_of_interest_marker_CY) <-c('celltype','markers')
genes_of_interest_marker_CY$weight <- 1;
##
genes_of_interest_marker_panglodb$source <- 'panglodb'
genes_of_interest_marker$source <- 'PMID'
genes_of_interest_marker_LKY$source <- 'LKY'
genes_of_interest_marker_CY$source <- 'CY'
genes_of_interest_marker_combine <- unique(rbind(genes_of_interest_marker_panglodb,
                                          genes_of_interest_marker,
                                          genes_of_interest_marker_LKY,
                                          genes_of_interest_marker_CY))
genes_of_interest_marker_combine$celltype <- sprintf('%s-%s',genes_of_interest_marker_combine$source,
                                                  genes_of_interest_marker_combine$celltype)
save(genes_of_interest_marker_combine,file='D:/analysis_eng/SherryPlot/data/genes_of_interest_marker_combine-Brain.RData')

