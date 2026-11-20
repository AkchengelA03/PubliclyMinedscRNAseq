#Loading Packages
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

#Tag Pre and Post Samples
table(Nir_so$PatientID)
Nir_so$Treatment <- ifelse(
  grepl("^Pre_", Nir_so$PatientID),
  "Naive",
  ifelse(
    grepl("^Post_", Nir_so$PatientID),
    "Post",
    NA
  )
)

#Remove Post Samples
seurat_obj <- subset(
  Nir_so,
  subset = !Treatment %in% c("Post")
)
table(seurat_obj$PatientID)

#Clustering
seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize")

seurat_obj <- FindVariableFeatures(seurat_obj )
seurat_obj  <- ScaleData(seurat_obj , verbose = FALSE)
seurat_obj  <- RunPCA(seurat_obj , npcs = 50, verbose = FALSE)
ElbowPlot(seurat_obj , ndims = 50)
seurat_obj  <- RunUMAP(seurat_obj , 
                      reduction = "pca", 
                      dims = 1:15)
seurat_obj  <- FindNeighbors(seurat_obj , 
                            reduction = "pca", 
                            dims = 1:15)
seurat_obj  <- FindClusters(object = seurat_obj ,
                           resolution = c(0.1,0.2,0.3,0.4,0.5, 0.6,  0.8, 1),
                           verbose = F)
seurat_obj @meta.data %>% 
  dplyr::select(contains("snn")) %>% 
  map_int(~ unique(.x) %>% length)

library(clustree)
clustree(seurat_obj , prefix = "RNA_snn_res.")

Idents(seurat_obj) <- "RNA_snn_res.0.2"

#UMAP
DimPlot(seurat_obj ,
        reduction = "umap",
        label = FALSE,
        label.size = 6)


#Feature Plot
FeaturePlot(seurat_obj, 
            reduction = "umap", 
            features = c("CD4","FOXP3","PDCD1","TOX"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = FALSE)


#Dotplpot
library(scRNAtoolVis)
jjDotPlot(seurat_obj ,
          gene = c("CD3G","CD4","CD8A","TCF7","LEF1","SELL","CCR7","FOXP3","CTLA4","PDCD1","TNFRSF18","IL2RA","CD14","FCGR3A","HLA-DRA","ITGAM","CD79A","CD19","MS4A1","IGHG1","IGHG2","TRGV9","TRDC","MKI67"),
          dot.col = c('blue', 'white', 'red'),
          tile.geom = F,
          ytree = F,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          id = "celltypes_PD1")


#DEG
all.markers <- FindAllMarkers(filtered_seurat , only.pos = TRUE)
DEG<- all.markers %>%
  filter(
    avg_log2FC > 0.585,
    p_val_adj < 0.05
  )


#Naming Clusters 1-8
seurat_obj @meta.data$celltypes[
  seurat_obj @meta.data$RNA_snn_res.0.2 == 0
] <- "C1"

#Separating Tregs based on PDCD1 gene expression
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj$celltypes_PD1 <- seurat_obj$celltypes
Idents(seurat_obj) <- "celltypes"         
cells_CD4 <- WhichCells(seurat_obj, idents = "Tregs")
PD1_vec <- LayerData(
  seurat_obj,
  assay = "RNA",
  layer = "counts"
)["PDCD1", cells_CD4]
seurat_obj$celltypes_PD1[cells_CD4[PD1_vec > 0]]  <- "PD1+ Tregs"
seurat_obj$celltypes_PD1[cells_CD4[PD1_vec == 0]] <- "PD1- Tregs"
Idents(seurat_obj) <- "celltypes_PD1"
table(Idents(seurat_obj))
table(seurat_obj$PatientID, seurat_obj$celltypes_PD1)

#Removing Samples with less than <3 PD1+ /PD1- Treg count
filtered_seurat <- subset(
  seurat_obj ,
  subset = !PatientID %in% c("Pre_P24","Pre_P35","Pre_P4","Pre_P8","Pre_P6"))

#Treg Seurat
Tregs <- subset(
  filtered_seurat ,
  subset = celltypes_PD1 %in% c("PD1+ Tregs","PD1- Tregs"))

#Vlnplot PD1- vs PD1+ Tregs
Idents(Tregs) <- "celltypes_PD1"
VlnPlot(
  Tregs,
  features = "",
  idents = c("PD1- Tregs", "PD1+ Tregs"),
  pt.size = 1,
  adjust = 3
)

#Dotplot for genes of interest
jjDotPlot(object = Tregs,
          gene = c("IL10","HAVCR2","KLRG1","MKI67","BCL6"),
          anno = T,
          plot.margin = c(3,1,1,1),
          point.geom = F,
          tile.geom = T,
          id = "celltypes_PD1")

#Statistical analysis using clusrank
library(clusrank)
gene <- "BCL6"
expression_data <- FetchData(object = Tregs, vars = c("celltypes_PD1", "PatientID", gene))
with(expression_data,pairwise.clusWilcox.test(
  expression_data$BCL6,
  group =celltypes_PD1,
  cluster = PatientID,
  method = "ds",
  p.adjust.method = "fdr"))

#Extracting mean mRNA expression per sample
avg_expr <- FetchData(Tregs, vars = c("IL10","KLRG1","HAVCR2","CD44","SELL","MKI67"), layer = "data")
avg_expr$celltypes_PD1<- Tregs$celltypes_PD1
avg_expr$PatientID<- Tregs$PatientID

avg_expr_summary <- avg_expr %>%
  group_by(celltypes_PD1, PatientID) %>%
  summarise(across(c("IL10","KLRG1","HAVCR2","CD44","SELL","MKI67"), mean), .groups = "drop")

#Exporting
library(writexl)
write_xlsx(avg_expr_summary,"mean mRNA.xlsx")


