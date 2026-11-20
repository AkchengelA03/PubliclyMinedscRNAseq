#Load Packages
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


expMat_C1 <- ReadMtx(mtx = "matrix.mtx", cells = "barcodes.tsv", features = "features.tsv",
                     cell.column = 1, feature.column = 2,
                     cell.sep = "\t", feature.sep = "\t",
                     skip.cell = 0, skip.feature = 0,
                     mtx.transpose = FALSE, unique.features = TRUE)

Mouse1 <- Seurat::CreateSeuratObject(counts = expMat_C1,
                                     min.cells = 0, min.features = 0, assay = "RNA")

expMat_C2 <- ReadMtx(mtx = "matrix.mtx", cells = "barcodes.tsv", features = "features.tsv",
                     cell.column = 1, feature.column = 2,
                     cell.sep = "\t", feature.sep = "\t",
                     skip.cell = 0, skip.feature = 0,
                     mtx.transpose = FALSE, unique.features = TRUE)

Mouse2 <- Seurat::CreateSeuratObject(counts = expMat_C2,
                                     min.cells = 0, min.features = 0, assay = "RNA")

expMat_C3 <- ReadMtx(mtx = "matrix.mtx", cells = "barcodes.tsv", features = "features.tsv",
                     cell.column = 1, feature.column = 2,
                     cell.sep = "\t", feature.sep = "\t",
                     skip.cell = 0, skip.feature = 0,
                     mtx.transpose = FALSE, unique.features = TRUE)

Mouse3 <- Seurat::CreateSeuratObject(counts = expMat_C3,
                                     min.cells = 0, min.features = 0, assay = "RNA")

expMat_C4 <- ReadMtx(mtx = "matrix.mtx", cells = "barcodes.tsv", features = "features.tsv",
                     cell.column = 1, feature.column = 2,
                     cell.sep = "\t", feature.sep = "\t",
                     skip.cell = 0, skip.feature = 0,
                     mtx.transpose = FALSE, unique.features = TRUE)

Mouse4 <- Seurat::CreateSeuratObject(counts = expMat_C4,
                                     min.cells = 0, min.features = 0, assay = "RNA")


count_matrices <- list(
  Mouse1 = expMat_C1,
  Mouse2 = expMat_C2,
  Mouse3 = expMat_C3,
  Mouse4 = expMat_C4
)

references <- list(
  Mouse1 = reference_Mouse1,
  Mouse2 = reference_Mouse2,
  Mouse3 = reference_Mouse3,
  Mouse4 = reference_Mouse4
)


#Assign Tissues
prepare_mouse <- function(counts, reference, sample_id) {
  
  antibody_features <- as.character(
    reference$name[
      reference$feature_type == "Antibody Capture"
    ]
  )
  
  tissue_tags <- as.character(
    reference$name[
      grepl("^TotalSeq_Hashtag", reference$id) &
        reference$feature_type == "Antibody Capture"
    ]
  )
  
  protein_features <- setdiff(antibody_features, tissue_tags)
  

  missing_features <- setdiff(
    antibody_features,
    rownames(counts)
  )
  
  if (length(missing_features) > 0) {
    stop(
      sample_id, ": missing antibody features: ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  stopifnot(length(tissue_tags) >= 2)
  
  gene_features <- setdiff(
    rownames(counts),
    antibody_features
  )
  
  obj <- CreateSeuratObject(
    counts = counts[gene_features, , drop = FALSE],
    assay = "RNA",
    min.cells = 0,
    min.features = 0
  )
  
  obj$Sample_ID <- sample_id
  
  obj[["HTO"]] <- CreateAssayObject(
    counts = counts[tissue_tags, , drop = FALSE]
  )
  
  if (length(protein_features) > 0) {
    obj[["ADT"]] <- CreateAssayObject(
      counts = counts[protein_features, , drop = FALSE]
    )
  }
  
  
  hto <- as.matrix(
    counts[tissue_tags, colnames(obj), drop = FALSE]
  )
  
  highest <- apply(hto, 2, max)
  
  second_highest <- apply(
    hto,
    2,
    function(x) sort(x, decreasing = TRUE)[2]
  )
  
  dominant_tag <- rownames(hto)[
    max.col(t(hto), ties.method = "first")
  ]
  
  passes <- highest > 0 & highest >= 30 * second_highest
  
  obj$Hashtag_max_count <- unname(highest)
  obj$Hashtag_second_count <- unname(second_highest)
  obj$Hashtag_30fold_pass <- unname(passes)
  
  obj$Tissue <- ifelse(
    passes,
    dominant_tag,
    NA_character_
  )
  

  if (!any(grepl("^mt-", rownames(obj[["RNA"]])))) {
    stop(sample_id, ": no mitochondrial genes matching '^mt-' found.")
  }
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    assay = "RNA",
    pattern = "^mt-"
  )
  
  obj$RNA_QC_pass <- (
    obj$nFeature_RNA > 200 &
      obj$percent.mt < 5
  )
  
  obj$Keep_cell <- (
    obj$Hashtag_30fold_pass &
      obj$RNA_QC_pass
  )
  
  return(obj)
}

mice_all <- lapply(names(count_matrices), function(id) {
  prepare_mouse(
    counts = count_matrices[[id]],
    reference = references[[id]],
    sample_id = id
  )
})

names(mice_all) <- names(count_matrices)

qc_summary <- do.call(
  rbind,
  lapply(names(mice_all), function(id) {
    
    obj <- mice_all[[id]]
    
    data.frame(
      Mouse = id,
      Starting_cells = ncol(obj),
      Pass_tissue_assignment = sum(obj$Hashtag_30fold_pass),
      Pass_RNA_QC = sum(obj$RNA_QC_pass),
      Pass_both = sum(obj$Keep_cell)
    )
  })
)

print(qc_summary)

Mouse1_all <- mice_all$Mouse1
Mouse2_all <- mice_all$Mouse2
Mouse3_all <- mice_all$Mouse3
Mouse4_all <- mice_all$Mouse4


if (any(qc_summary$Pass_both == 0)) {
  stop("At least one mouse has no passing cells. Inspect qc_summary.")
}

mice_filtered <- lapply(mice_all, function(obj) {
  subset(obj, subset = Keep_cell == TRUE)
})

Mouse1 <- mice_filtered$Mouse1
Mouse2 <- mice_filtered$Mouse2
Mouse3 <- mice_filtered$Mouse3
Mouse4 <- mice_filtered$Mouse4

#Merge Seurats
seurat_Merged <- merge(
  x = Mouse1,
  y = list(Mouse2, Mouse3, Mouse4),
  merge.data = FALSE
)

DefaultAssay(seurat_Merged) <- "RNA"

if (inherits(seurat_Merged[["RNA"]], "Assay5")) {
  seurat_Merged <- JoinLayers(
    seurat_Merged,
    assay = "RNA"
  )
}

stopifnot(
  ncol(seurat_Merged) == sum(qc_summary$Pass_both)
)

#Combine Sample ID and Tissue
seurat_Merged$Mouse_Tissue <- paste(
  seurat_Merged$Sample_ID,
  seurat_Merged$Tissue,
  sep = "_"
)


#Clustering
seurat_Merged <- NormalizeData(seurat_Merged)

seurat_Merged <- FindVariableFeatures(
  seurat_Merged,
  selection.method = "vst",
  nfeatures = 2000
)

seurat_Merged <- ScaleData(
  seurat_Merged,
  features = VariableFeatures(seurat_Merged)
)

seurat_Merged <- RunPCA(
  seurat_Merged,
  features = VariableFeatures(seurat_Merged),
  npcs = 50,
  verbose = FALSE
)

seurat_Merged <- RunHarmony(
  object = seurat_Merged,
  group.by.vars = "Mouse_Tissue",
  reduction.use = "pca"
)

ElbowPlot(seurat_Merged , ndims = 50)

seurat_Merged <- RunUMAP(
  seurat_Merged,
  reduction = "harmony",
  dims = 1:20,
  min.dist = 0.15,
  spread = 1
)

seurat_Merged <- FindNeighbors(
  seurat_Merged,
  reduction = "harmony",
  dims = 1:20,
  graph.name = c("harmony_nn", "harmony_snn")
)

seurat_Merged <- FindClusters(
  seurat_Merged,
  graph.name = "harmony_snn",
  resolution = c(0.1, 0.2, 0.3,0.4,0.5 ,0.6, 0.7, 0.8, 1, 1.2, 1.4),
  verbose = FALSE
)

seurat_Merged@meta.data %>% 
  dplyr::select(contains("harmony_snn_res")) %>% 
  purrr::map_int(~length(unique(.x)))

#Removing Multiplets
seurat_Merged <- subset(
  seurat_Merged,
  subset = !harmony_snn_res.0.4 %in% c("9","10")
)


#Recluster
seurat_Merged <- NormalizeData(seurat_Merged)

seurat_Merged <- FindVariableFeatures(
  seurat_Merged,
  selection.method = "vst",
  nfeatures = 2000
)

seurat_Merged <- ScaleData(
  seurat_Merged,
  features = VariableFeatures(seurat_Merged)
)

seurat_Merged <- RunPCA(
  seurat_Merged,
  features = VariableFeatures(seurat_Merged),
  npcs = 50,
  verbose = FALSE
)

seurat_Merged <- RunHarmony(
  object = seurat_Merged,
  group.by.vars = "Mouse_Tissue",
  reduction.use = "pca"
)

ElbowPlot(seurat_Merged , ndims = 50)

seurat_Merged <- RunUMAP(
  seurat_Merged,
  reduction = "harmony",
  dims = 1:10,
  min.dist = 0.15,
  spread = 1
)

seurat_Merged <- FindNeighbors(
  seurat_Merged,
  reduction = "harmony",
  dims = 1:10,
  graph.name = c("harmony_nn", "harmony_snn")
)

seurat_Merged <- FindClusters(
  seurat_Merged,
  graph.name = "harmony_snn",
  resolution = c(0.1, 0.2, 0.3,0.4,0.5 ,0.6, 0.7, 0.8, 1, 1.2, 1.4),
  verbose = FALSE
)

seurat_Merged@meta.data %>% 
  dplyr::select(contains("harmony_snn_res")) %>% 
  purrr::map_int(~length(unique(.x)))


Idents(seurat_Merged) <- "harmony_snn_res.0.4"

DimPlot(
  seurat_Merged,
  reduction = "umap",
  label = FALSE,
  repel = TRUE
) 


#UMAP
DimPlot(
  seurat_Merged,
  reduction = "umap",
  group.by = "cellidentity",
  label = FALSE
) 

#UMAP Per Sample
DimPlot(
  seurat_Merged,
  reduction = "umap",
  group.by = "cellidentity",
  split.by = "Sample_ID",
  label = FALSE,
  ncol = 3
)

#UMAP Per Tissue
DimPlot(
  seurat_Merged,
  reduction = "umap",
  group.by = "cellidentity",
  split.by = "Tissue",
  label = FALSE,
  ncol = 3
) 
#Feature Plot
FeaturePlot(seurat_Merged, 
            reduction = "umap", 
            features = c("Pdcd1","Tox","Tcf7"), 
            order = TRUE,
            min.cutoff = 'q10',
            ncol= 2,
            label = FALSE)

#Dot Plot
library(scRNAtoolVis)
jjDotPlot(
  seurat_Merged,
  gene = c(
    "Tcf7", "Lef1", "Bach2", "Sell", "Ccr7","Id3","Slamf6","Tox", "Pdcd1","Tnfrsf4","Tnfrsf9","Tnfrsf18",
     "Nt5e", "Entpd1", "Klrg1", "Il7r","Il2ra","Nrp1", "Ikzf2",
    "Batf","Prdm1","Bcl6","Rorc", "Tbx21","Ifng", "Il4",
    "Tnf","Il10","Gzmb", "Areg", "Rorc", "Il17a", "Il22",
    "Isg15", "Isg20", "Irf7", "Mki67", "Stmn1", "Tyms"
  ),
  dot.col = c("blue", "white", "red"),
  tile.geom = FALSE,
  ytree = FALSE,
  rescale.min = -2,
  rescale.max = 2,
  midpoint = 0,
  id = "cellidentity"
) 

#DEG
all.markers <- FindAllMarkers(seurat_Merged, only.pos = TRUE)
DEG.filtered <- all.markers %>%
  filter(
    avg_log2FC > 0.585,
    p_val_adj < 0.05
  )

#Exporting DEG
library(writexl)
write_xlsx(DEG.filtered,"DEG.xlsx")


#Naming clusters applied for clusters 0-8
seurat_Merged@meta.data$cellidentity[
  seurat_Merged@meta.data$harmony_snn_res.0.4 == 0
] <- "Treg C1"


#Frequency Plot
library(dplyr)
library(ggplot2)

tbl <- table(
  seurat_Merged@meta.data$Mouse_Tissue,
  seurat_Merged@meta.data$cellidentity
)

df <- as.data.frame(tbl)
colnames(df) <- c("Mouse_Tissue", "id", "freq")

df <- df %>%
  group_by(Mouse_Tissue) %>%
  mutate(freq_pct = freq / sum(freq)) %>%
  ungroup()


ggplot(df, aes(Mouse_Tissue, freq_pct, fill = id)) +
  geom_col() +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = NULL, y = "Frequency (%)", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Separating Tissues (Tissue= Pancreas/LPL/Kidney)
Tissue <- subset(
  seurat_Merged,
  subset = Tissue %in% c("Tissue")
)

#Pseudotime Trajectory Analysis
# Pseudotime analysis for one tissue
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(patchwork)

set.seed(1234)

seu <- Tissue

dims_use <- 1:10
cluster_column <- "cellidentity"

naive_markers <- c(
  "Tcf7", "Bach2", "Lef1", "Ccr7", "Sell"
)

cells <- colnames(seu)

DefaultAssay(seu) <- "RNA"

sce <- as.SingleCellExperiment(
  seu,
  assay = "RNA"
)

sce <- sce[, cells]


reducedDim(sce, "DR") <- Embeddings(
  seu, "harmony"
)[cells, dims_use, drop = FALSE]

reducedDim(sce, "UMAP") <- Embeddings(
  seu, "umap"
)[cells, 1:2, drop = FALSE]

sce$cluster <- droplevels(
  factor(seu[[]][cells, cluster_column])
)

naive_present <- intersect(
  naive_markers,
  rownames(sce)
)

naive_missing <- setdiff(
  naive_markers,
  naive_present
)
naive_score <- Matrix::colMeans(
  logcounts(sce)[naive_present, , drop = FALSE]
)

cluster_naive_scores <- tapply(
  naive_score,
  sce$cluster,
  mean
)

start_cluster <- names(
  which.max(cluster_naive_scores)
)

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim = "DR",
  start.clus = start_cluster,
  approx_points = 150,
  extend = "n",
  allow.breaks = TRUE
)

sling_tissue <- sce

pt_mat <- slingPseudotime(
  sling_tissue,
  na = TRUE
)

weight_mat <- slingCurveWeights(sling_tissue)

avg_pt <- as.numeric(
  slingAvgPseudotime(sling_tissue)
)

avg_pt[!is.finite(avg_pt)] <- NA_real_

rownames(pt_mat) <- colnames(sling_tissue)
rownames(weight_mat) <- colnames(sling_tissue)

colnames(pt_mat) <- paste0(
  "tissue_Slingshot_PT",
  seq_len(ncol(pt_mat))
)

colnames(weight_mat) <- paste0(
  "tissue_Slingshot_weight",
  seq_len(ncol(weight_mat))
)

new_metadata <- data.frame(
  pt_mat,
  weight_mat,
  tissue_Slingshot_AvgPseudotime = avg_pt,
  row.names = colnames(sling_tissue),
  check.names = FALSE
)

Tissue <- AddMetaData(
  Tissue,
  metadata = new_metadata
)

# UMAPs with Slingshot trajectories
umap <- reducedDim(sling_tissue, "UMAP")

sling_tissue_umap <- embedCurves(
  sling_tissue,
  newDimRed = umap
)

curves <- slingCurves(sling_tissue_umap)

curve_df <- do.call(
  rbind,
  lapply(seq_along(curves), function(i) {
    
    cc <- curves[[i]]
    coords <- cc$s[cc$ord, , drop = FALSE]
    
    data.frame(
      UMAP_1 = coords[, 1],
      UMAP_2 = coords[, 2],
      Lineage = paste0("Lineage ", i)
    )
  })
)

umap_df <- data.frame(
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  Cluster = sling_tissue$cluster,
  Pseudotime = avg_pt,
  row.names = colnames(sling_tissue)
)

p_pt <- ggplot(
  umap_df,
  aes(UMAP_1, UMAP_2)
) +
  geom_point(
    aes(colour = Pseudotime),
    size = 0.6
  ) +
  geom_path(
    data = curve_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      group = Lineage
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.7
  ) +
  scale_colour_gradientn(
    colours = c("#440154", "#21908C", "#FDE725"),
    na.value = "grey85",
    name = "Average\npseudotime"
  ) +
  coord_equal() +
  theme_classic(base_size = 12) +
  labs(
    title = "tissue: Slingshot pseudotime",
    x = "UMAP_1",
    y = "UMAP_2"
  )

p_cl <- ggplot(
  umap_df,
  aes(UMAP_1, UMAP_2)
) +
  geom_point(
    aes(colour = Cluster),
    size = 0.6
  ) +
  geom_path(
    data = curve_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      group = Lineage
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.7
  ) +
  coord_equal() +
  theme_classic(base_size = 12) +
  labs(
    title = "tissue: clusters",
    x = "UMAP_1",
    y = "UMAP_2",
    colour = "Cluster"
  )

print(p_pt | p_cl)

# Pseudotime UMAP for each lineage
pt_plots_tissue <- lapply(
  seq_len(ncol(pt_mat)),
  function(i) {
    
    df <- umap_df
    df$LineagePseudotime <- pt_mat[, i]
    
    assigned_df <- df[
      is.finite(df$LineagePseudotime),
      ,
      drop = FALSE
    ]
    
    this_curve <- curve_df[
      curve_df$Lineage == paste0("Lineage ", i),
      ,
      drop = FALSE
    ]
    
    ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(
        colour = "grey85",
        size = 0.6
      ) +
      geom_point(
        data = assigned_df,
        aes(colour = LineagePseudotime),
        size = 0.6
      ) +
      geom_path(
        data = this_curve,
        aes(
          x = UMAP_1,
          y = UMAP_2,
          group = Lineage
        ),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.7
      ) +
      scale_colour_gradientn(
        colours = c("#440154", "#21908C", "#FDE725"),
        name = "Pseudotime"
      ) +
      coord_equal() +
      theme_classic(base_size = 12) +
      labs(
        title = paste("tissue: Lineage", i),
        x = "UMAP_1",
        y = "UMAP_2"
      )
  }
)

print(wrap_plots(pt_plots_tissue, ncol = 2))
Idents(tissue)= "cellidentity"
seu <- tissue
lineages_use <- 1:3

genes_use <- c(
  "Tcf7", "Lef1", "Bach2", "Sell", "Ccr7",
  "Pdcd1", "Tox", "Havcr2", "Tigit",
  "Tnfrsf9", "Tnfrsf18", "Icos","Il2ra",
  "Foxp3", "Ikzf2"
)

expression_colours <- c(
  "#053061", "#2166AC", "#4393C3",
  "#92C5DE", "#D1E5F0", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D",
  "#B2182B", "#67001F"
)

pseudotime_colours <- c(
  "#440154", "#21908C", "#FDE725"
)

lineages_tissue <- slingLineages(sling_tissue)

meta <- seu[[]]

pt_columns <- paste0(
  "tissue_Slingshot_PT",
  lineages_use
)


DefaultAssay(seu) <- "RNA"


genes_present <- genes_use[
  genes_use %in% rownames(seu[["RNA"]])
]

missing_genes <- setdiff(genes_use, genes_present)


expr_all <- FetchData(
  seu,
  vars = genes_present,
  cells = rownames(meta),
  layer = "data",
  clean = FALSE
)

expr_all <- as.matrix(
  expr_all[rownames(meta), genes_present, drop = FALSE]
)

stopifnot(all(is.finite(expr_all)))

# Gene-expression heatmap for each lineage
make_tissue_heatmap <- function(lineage_number) {
  
  pt_column <- paste0(
    "tissue_Slingshot_PT",
    lineage_number
  )
  
  lineage_clusters <- lineages_tissue[[lineage_number]]
  
  cell_info <- data.frame(
    Cell = rownames(meta),
    Cluster = as.character(meta[[cluster_column]]),
    Pseudotime = meta[[pt_column]],
    stringsAsFactors = FALSE
  )
  
  cell_info <- cell_info[
    is.finite(cell_info$Pseudotime) &
      !is.na(cell_info$Cluster) &
      cell_info$Cluster %in% lineage_clusters,
    ,
    drop = FALSE
  ]
  
  cluster_order <- lineage_clusters[
    lineage_clusters %in% cell_info$Cluster
  ]
  
  n_clusters <- length(cluster_order)
  
  expr <- expr_all[
    cell_info$Cell,
    genes_present,
    drop = FALSE
  ]
  
  cluster_means <- do.call(
    cbind,
    lapply(cluster_order, function(cl) {
      colMeans(
        expr[cell_info$Cluster == cl, , drop = FALSE]
      )
    })
  )
  
  rownames(cluster_means) <- genes_present
  colnames(cluster_means) <- cluster_order
  
  cluster_info <- data.frame(
    Position = seq_len(n_clusters),
    Cluster = cluster_order,
    Pseudotime = vapply(
      cluster_order,
      function(cl) {
        mean(cell_info$Pseudotime[cell_info$Cluster == cl])
      },
      numeric(1)
    )
  )
  
  gene_sd <- apply(cluster_means, 1, sd)
  gene_sd[!is.finite(gene_sd) | gene_sd == 0] <- 1
  
  z_mat <- sweep(
    cluster_means,
    1,
    rowMeans(cluster_means),
    "-"
  )
  
  z_mat <- sweep(
    z_mat,
    1,
    gene_sd,
    "/"
  )
  
  z_mat[z_mat < -2] <- -2
  z_mat[z_mat > 2] <- 2
  
  heatmap_df <- data.frame(
    Gene = rep(
      rownames(z_mat),
      times = ncol(z_mat)
    ),
    Position = rep(
      seq_len(ncol(z_mat)),
      each = nrow(z_mat)
    ),
    Expression = as.vector(z_mat)
  )
  
  heatmap_df$Gene <- factor(
    heatmap_df$Gene,
    levels = rev(genes_present)
  )
  
  strip_edges <- seq(
    0.5,
    n_clusters + 0.5,
    length.out = 1001
  )
  
  strip_midpoints <- (
    head(strip_edges, -1) +
      tail(strip_edges, -1)
  ) / 2
  
  strip_df <- data.frame(
    xmin = head(strip_edges, -1),
    xmax = tail(strip_edges, -1),
    Pseudotime = approx(
      x = cluster_info$Position,
      y = cluster_info$Pseudotime,
      xout = strip_midpoints,
      method = "linear",
      rule = 2
    )$y
  )
  
  p_top <- ggplot(strip_df) +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 0,
        ymax = 1,
        fill = Pseudotime
      ),
      colour = NA
    ) +
    scale_fill_gradientn(
      colours = pseudotime_colours,
      name = "Pseudotime",
      guide = guide_colourbar(
        barheight = grid::unit(22, "mm"),
        barwidth = grid::unit(3, "mm")
      )
    ) +
    scale_x_continuous(
      limits = c(0.5, n_clusters + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    labs(
      title = paste("Lineage", lineage_number)
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(size = 14),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  p_expression <- ggplot(
    heatmap_df,
    aes(Position, Gene, fill = Expression)
  ) +
    geom_tile(
      width = 1,
      height = 1
    ) +
    scale_fill_gradientn(
      colours = expression_colours,
      limits = c(-2, 2),
      breaks = c(-2, -1, 0, 1, 2),
      name = NULL,
      guide = guide_colourbar(
        barheight = grid::unit(25, "mm"),
        barwidth = grid::unit(3, "mm")
      )
    ) +
    scale_x_continuous(
      limits = c(0.5, n_clusters + 0.5),
      breaks = cluster_info$Position,
      labels = cluster_info$Cluster,
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(
        face = "italic",
        colour = "black",
        size = 11
      ),
      axis.text.x = element_text(
        colour = "black",
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 10,
        margin = margin(t = 4)
      ),
      legend.text = element_text(size = 9)
    )
  
  (p_top / p_expression) +
    plot_layout(
      heights = c(0.65, 7)
    )
}

tissue_heatmaps <- setNames(
  lapply(
    lineages_use,
    make_tissue_heatmap
  ),
  paste0("Lineage", lineages_use)
)

print(tissue_heatmaps[["Lineage1"]])
print(tissue_heatmaps[["Lineage2"]])
print(tissue_heatmaps[["Lineage3"]])
