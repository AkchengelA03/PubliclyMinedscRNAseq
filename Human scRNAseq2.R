#Loading packages
library (DESeq2)
library(tidyverse)
library(Seurat)
library(BiocManager)
library(Matrix)
library(RCurl)
library(readxl)
library(scales)
library(ggpubr)
library(data.table)
library(AnnotationHub)
library(ensembldb)
library(HGNChelper)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(data.table)
library(dplyr)
library(harmony)

#Loading seurat
seurat_obj<- readRDS(
  "melanoma_patients_TIL_CD4_seurat_object.rds"
)

#Filtering Samples
table(TIL_CD4$Patient_ID)
filtered_seurat <- subset(
  TIL_CD4,
  subset = !Patient_ID %in% c("Mel-T-02","Mel-T-06","Mel-T-07","Mel-T-09","Mel-T-10","Mel-T-11","Mel-T-12")
)
table(filtered_seurat$Patient_ID)

library(Seurat)
library(harmony)

#Clustering
DefaultAssay(filtered_seurat) <- "RNA"

filtered_seurat <- FindVariableFeatures(
  filtered_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

filtered_seurat <- ScaleData(
  filtered_seurat,
  features = VariableFeatures(filtered_seurat)
)

filtered_seurat <- RunPCA(
  filtered_seurat,
  features = VariableFeatures(filtered_seurat),
  npcs = 50,
  verbose = FALSE
)

filtered_seurat <- RunHarmony(
  object = filtered_seurat,
  group.by.vars = "Patient_ID",
  reduction.use = "pca"
)

ElbowPlot(filtered_seurat , ndims = 50)

filtered_seurat <- RunUMAP(
  filtered_seurat,
  reduction = "harmony",
  dims = 1:10,
  min.dist = 0.15,
  spread = 1
)

filtered_seurat <- FindNeighbors(
  filtered_seurat,
  reduction = "harmony",
  dims = 1:10,
  graph.name = c("harmony_nn", "harmony_snn")
)

filtered_seurat <- FindClusters(
  filtered_seurat,
  graph.name = "harmony_snn",
  resolution = c(0.1, 0.2, 0.3,0.4,0.5 ,0.6, 0.7, 0.8, 1, 1.2, 1.4),
  verbose = FALSE
)

filtered_seurat@meta.data %>% 
  dplyr::select(contains("harmony_snn_res")) %>% 
  purrr::map_int(~length(unique(.x)))
library(clustree)
clustree(filtered_seurat, prefix = "harmony_snn_res.")

Idents(filtered_seurat) <- "celltypes"


DimPlot(
    filtered_seurat,
    reduction = "umap",
    split.by = "Prior_immune_checkpoint_blockade",
    label = FALSE,
    repel = TRUE
  )


#Featureplot
FeaturePlot(filtered_seurat, 
            reduction = "umap", 
            features = c(""), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = FALSE)


#DEG
all.markers <- FindAllMarkers(filtered_seurat , only.pos = TRUE)
DEG<- all.markers %>%
  filter(
    avg_log2FC > 0.585,
    p_val_adj < 0.05
  )

library(writexl)
write_xlsx(DEG," DEG.xlsx")


#Dotplot
library(scRNAtoolVis)
jjDotPlot(filtered_seurat,
          gene = c("CD4","CD8A","FOXP3","IL2RA","CTLA4","TNFRSF18","TBX21","EOMES","GZMB","TCF7","LEF1","CD44","SELL","CCR7","KLF2","KLRG1","IL7R",
                   "MKI67"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = "celltypes")


#Naming Clusters 0-7
filtered_seurat@meta.data$celltypes[
  filtered_seurat@meta.data$harmony_snn_res.0.4 == 0
] <- "C1"


#Treg seurat
Tregs <- subset(
  filtered_seurat,
  subset = harmony_snn_res.0.4 %in% c("2")
)


#TOX+ Freq
tox_expression <- FetchData(
  Tregs,
  vars = "TOX",
  layer = "data"
)

Tregs$TOX_status <- ifelse(
  tox_expression$TOX > 0,
  "TOX+",
  "TOX−"
)

table(Tregs$TOX_status)
prop.table(table(Tregs$TOX_status)) * 100


tox_frequency <- Tregs@meta.data %>%
  group_by(
    Patient_ID,
    Prior_immune_checkpoint_blockade
  ) %>%
  summarise(
    total_Tregs = n(),
    TOX_positive = sum(TOX_status == "TOX+"),
    TOX_frequency = 100 * TOX_positive / total_Tregs,
    .groups = "drop"
  )

tox_frequency

write_xlsx(tox_frequency,"toxfreq.xlsx")

#TCF7+ TOX+ and TCF7- TOX+ Frequency
gene_expression <- FetchData(
  Tregs,
  vars = c("TCF7", "TOX"),
  layer = "data"
)

Tregs$TCF7_TOX_status <- case_when(
  gene_expression$TCF7 > 0 & gene_expression$TOX > 0 ~ "TCF7+TOX+",
  gene_expression$TCF7 <= 0 & gene_expression$TOX > 0 ~ "TCF7-TOX+",
  gene_expression$TCF7 > 0 & gene_expression$TOX <= 0 ~ "TCF7+TOX-",
  TRUE                                                ~ "TCF7-TOX-"
)

TCF7_TOX_table <- Tregs@meta.data %>%
  group_by(
    Patient_ID,
    Prior_immune_checkpoint_blockade
  ) %>%
  summarise(
    Total_Tregs = n(),
    
    `TCF7+TOX+_number` = sum(
      TCF7_TOX_status == "TCF7+TOX+",
      na.rm = TRUE
    ),
    
    `TCF7-TOX+_number` = sum(
      TCF7_TOX_status == "TCF7-TOX+",
      na.rm = TRUE
    ),
    
    `TCF7+TOX+_%` = round(
      100 * `TCF7+TOX+_number` / Total_Tregs,
      2
    ),
    
    `TCF7-TOX+_%` = round(
      100 * `TCF7-TOX+_number` / Total_Tregs,
      2
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    Prior_immune_checkpoint_blockade,
    Patient_ID
  )

#TCF7- TOX+ Analysis
TOX_seurat <- subset(
  Tregs,
  subset = TCF7_TOX_status %in% c("TCF7-TOX+")
)

#Filtering samples with <3 cells
filtered_TOX_seurat <- subset(
  TOX_seurat,
  subset = !Patient_ID %in% c("Mel-UT-14","Mel-UT-16","Mel-UT-17")
)

#Heatmap for TCF1- TOX+
Idents(filtered_TOX_seurat) <- "Prior_immune_checkpoint_blockade"
jjDotPlot(object = filtered_TOX_seurat,
          gene = c("IL10","HAVCR2","KLRG1","MKI67","BCL6"),
          anno = T,
          plot.margin = c(3,1,1,1),
          point.geom = F,
          tile.geom = T,
          id = "Prior_immune_checkpoint_blockade")


#
Idents(filtered_TOX_seurats) <- "Prior_immune_checkpoint_blockade"

#Statistical analysis using Clusrank
library(clusrank)
gene <- "BCL6"
expression_data <- FetchData(object = filtered_TOX_seurat, vars = c("Prior_immune_checkpoint_blockade", "Patient_ID", gene))
with(expression_data,pairwise.clusWilcox.test(
  expression_data$BCL6,
  group = Prior_immune_checkpoint_blockade,
  cluster = Patient_ID,
  method = "ds",
  p.adjust.method = "fdr"))

#mean mRNA expression
avg_expr <- FetchData(filtered_TOX_seurat, vars = c("IL10","HAVCR2","KLRG1","MKI67","BCL6"), layer = "data")
avg_expr$Prior_immune_checkpoint_blockade<- filtered_TOX_seurat$Prior_immune_checkpoint_blockade
avg_expr$Patient_ID<- filtered_TOX_seurat$Patient_ID

avg_expr_summary <- avg_expr %>%
  group_by(Prior_immune_checkpoint_blockade, Patient_ID) %>%
  summarise(across(c("IL10","HAVCR2","KLRG1","MKI67","BCL6"), mean), .groups = "drop")

#Exporting 
library(writexl)
write_xlsx(avg_expr_summary," mean mRNA.xlsx")

