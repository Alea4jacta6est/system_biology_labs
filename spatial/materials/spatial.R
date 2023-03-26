library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


lymph_node <- Load10X_Spatial(data.dir = "~/Code/sys_bio/Day 5 - March 22th/spatial/materials", filename = "spatial/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5")
lymph_node

VlnPlot(lymph_node, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(lymph_node, features = "nCount_Spatial") + theme(legend.position = "right")

lymph_node <- SCTransform(lymph_node, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(lymph_node, features = c("MKI67", "TUBB"), alpha = c(0.1, 1))
SpatialFeaturePlot(lymph_node, features = c("CD3D", "IL7R"), alpha = c(0.1, 1))
SpatialFeaturePlot(lymph_node, features = c("CD19", "MS4A1", "CD21"), alpha = c(0.1, 1))

lymph_node <- RunPCA(lymph_node, assay = "SCT", verbose = FALSE)
lymph_node <- FindNeighbors(lymph_node, reduction = "pca", dims = 1:30)
lymph_node <- FindClusters(lymph_node, verbose = FALSE)
lymph_node <- RunUMAP(lymph_node, reduction = "pca", dims = 1:30)

DimPlot(lymph_node, reduction = "umap")
SpatialDimPlot(lymph_node, label = TRUE, label.size = 3)

SpatialPlot(lymph_node, pt.size.factor = 0.6)
SpatialDimPlot(lymph_node, label = TRUE, label.size = 3)

de_markers <- FindAllMarkers(lymph_node, only.pos = T)
head(de_markers)


goodMarkers <- de_markers %>% filter(cluster == 0) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(0)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 1) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(1)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 2) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(2)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 3) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(3)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 4) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(4)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 5) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(5)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 6) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(6)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 7) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(7)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)

goodMarkers <- de_markers %>% filter(cluster == 8) %>% slice_min(order_by = p_val_adj, n = 4) %>% pull(gene)
SpatialDimPlot(lymph_node, cells.highlight = CellsByIdentities(object = lymph_node, idents = c(8)), facet.highlight = TRUE)
SpatialFeaturePlot(object = lymph_node, features = goodMarkers, alpha = c(0.1, 1), ncol = 2)


# Install Rfast2 first
# lymph_node <- FindSpatiallyVariableFeatures(lymph_node, assay = "SCT", features = VariableFeatures(lymph_node)[1:1000], selection.method = "moransi")
# SpatialFeaturePlot(lymph_node, features = head(SpatiallyVariableFeatures(lymph_node, selection.method = "moransi"), 6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")

