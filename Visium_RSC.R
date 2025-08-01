# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                          RSC Visium data                           -----
# -----                                                                    -----
# -----                           Chatterjee Lab                           -----
# -----                         University of Iowa                         -----
# -----                                                                    -----
# ------------------------------------------------------------------------------

# Summary: RSC Spatial Transcriptomics data (Homecage vs SOR)
# Written by: Budhaditya Basu
# Date: 7/30/2024
install.packages('Seurat')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(openxlsx)
library(scCustomize)
library(ggthemes)
library(ggrepel)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  SOR1 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOR.1.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample01SOR/"
SOR.1 <- Load10X_Spatial(
  SOR.1.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "SOR1")
head(SOR.1@meta.data)
SOR.1$orig.ident <- "SOR1"

SpatialDimPlot(SOR.1, interactive = TRUE)

plot1 <- SpatialFeaturePlot(SOR.1, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("SOR1")
VlnPlot(SOR.1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(SOR.1, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR1.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  SOR2 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR.2.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample02SOR/"
SOR.2 <- Load10X_Spatial(
  SOR.2.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "SOR2")
head(SOR.2@meta.data)
SOR.2$orig.ident <- "SOR2"

plot2 <- SpatialFeaturePlot(SOR.2, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("SOR2")
VlnPlot(SOR.2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(SOR.2, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR2.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  SOR3 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR.3.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample03SOR/"
SOR.3 <- Load10X_Spatial(
  SOR.3.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "SOR3")
head(SOR.3@meta.data)
SOR.3$orig.ident <- "SOR3"

plot3 <- SpatialFeaturePlot(SOR.3, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("SOR3")
VlnPlot(SOR.3, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(SOR.3, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR3.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  SOR4 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR.4.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample04SOR/"
SOR.4 <- Load10X_Spatial(
  SOR.4.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "SOR4")
head(SOR.4@meta.data)
SOR.4$orig.ident <- "SOR4"

plot4 <- SpatialFeaturePlot(SOR.4, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("SOR4")
VlnPlot(SOR.4, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(SOR.4, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR4.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  SOR5 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR.5.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample10SOR/"
SOR.5 <- Load10X_Spatial(
  SOR.5.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "SOR5")
head(SOR.5@meta.data)
SOR.5$orig.ident <- "SOR5"

plot5 <- SpatialFeaturePlot(SOR.5, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("SOR5")
VlnPlot(SOR.5, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(SOR.5, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR5.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  SOR6 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR.6.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample14SOR/"
SOR.6 <- Load10X_Spatial(
  SOR.6.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "SOR6")
head(SOR.6@meta.data)
SOR.6$orig.ident <- "SOR6"

plot6 <- SpatialFeaturePlot(SOR.6, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("SOR6")
VlnPlot(SOR.6, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(SOR.6, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR6.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  SOR7 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR.7.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample16SOR/"
SOR.7 <- Load10X_Spatial(
  SOR.7.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "SOR7")
head(SOR.7@meta.data)
SOR.7$orig.ident <- "SOR7"

plot7 <- SpatialFeaturePlot(SOR.7, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("SOR7")
VlnPlot(SOR.7, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(SOR.7, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR7.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  HC1 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC.1.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample05HC/"
HC.1 <- Load10X_Spatial(
  HC.1.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "HC1")
head(HC.1@meta.data)
HC.1$orig.ident <- "HC1"

plot8 <- SpatialFeaturePlot(HC.1, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("HC1")
VlnPlot(HC.1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(HC.1, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC1.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  HC2 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC.2.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample06HC/"
HC.2 <- Load10X_Spatial(
  HC.2.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "HC2")
head(HC.2@meta.data)
HC.2$orig.ident <- "HC2"

plot9 <- SpatialFeaturePlot(HC.2, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("HC2")
VlnPlot(HC.2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(HC.2, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC2.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  HC3 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC.3.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample07HC/"
HC.3 <- Load10X_Spatial(
  HC.3.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "HC3")
head(HC.3@meta.data)
HC.3$orig.ident <- "HC3"

plot10 <- SpatialFeaturePlot(HC.3, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("HC3")
VlnPlot(HC.3, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(HC.3, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC3.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  HC4 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC.4.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample08HC/"
HC.4 <- Load10X_Spatial(
  HC.4.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "HC4")
head(HC.4@meta.data)
HC.4$orig.ident <- "HC4"

plot11 <- SpatialFeaturePlot(HC.4, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("HC4")
VlnPlot(HC.4, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(HC.4, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC4.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  HC5 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC.5.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample09HC/"
HC.5 <- Load10X_Spatial(
  HC.5.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "HC5")
head(HC.5@meta.data)
HC.5$orig.ident <- "HC5"

plot10 <- SpatialFeaturePlot(HC.5, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("HC5")
VlnPlot(HC.5, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(HC.5, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC5.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  HC6 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC.6.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample13HC/"
HC.6 <- Load10X_Spatial(
  HC.6.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "HC6")
head(HC.6@meta.data)
HC.6$orig.ident <- "HC6"

plot10 <- SpatialFeaturePlot(HC.6, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("HC6")
VlnPlot(HC.6, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(HC.6, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC6.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                  HC7 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HC.7.data.dir <- "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/Sample15HC/"
HC.7 <- Load10X_Spatial(
  HC.7.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "HC7")
head(HC.7@meta.data)
HC.7$orig.ident <- "HC7"

plot10 <- SpatialFeaturePlot(HC.7, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("HC7")
VlnPlot(HC.7, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(HC.7, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC7.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import the biological replicates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR1 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR1.rds")
SOR2 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR2.rds")
SOR3 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR3.rds")
SOR4 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR4.rds")
SOR5 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR5.rds")
SOR6 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR6.rds")
SOR7 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SOR7.rds")
HC1 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC1.rds")
HC2 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC2.rds")
HC3 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC3.rds")
HC4 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC4.rds")
HC5 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC5.rds")
HC6 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC6.rds")
HC7 <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/HC7.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subset RSC
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Subset the RSC region from the biological replicates
SOR1.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample01SOR.csv") %>% pull(Barcode)
SOR2.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample02SOR.csv") %>% pull(Barcode)
SOR3.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample03SOR.csv") %>% pull(Barcode)
SOR4.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample04SOR.csv") %>% pull(Barcode)
SOR5.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample10SOR.csv") %>% pull(Barcode)
SOR6.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample14SOR.csv") %>% pull(Barcode)
SOR7.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample16SOR.csv") %>% pull(Barcode)

HC1.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample05HC.csv") %>% pull(Barcode)
HC2.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample06HC.csv") %>% pull(Barcode)
HC3.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample07HC.csv") %>% pull(Barcode)
HC4.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample08HC.csv") %>% pull(Barcode)
HC5.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample09HC.csv") %>% pull(Barcode)
HC6.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample13HC.csv") %>% pull(Barcode)
HC7.RSC <- read.csv("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/rawdata/RSC.sample15HC.csv") %>% pull(Barcode)

# Subset 
SOR1 <- subset(SOR1, cells = SOR1.RSC)
SOR2 <- subset(SOR2, cells = SOR2.RSC)
SOR3 <- subset(SOR3, cells = SOR3.RSC)
SOR4 <- subset(SOR4, cells = SOR4.RSC)
SOR5 <- subset(SOR5, cells = SOR5.RSC)
SOR6 <- subset(SOR6, cells = SOR6.RSC)
SOR7 <- subset(SOR7, cells = SOR7.RSC)

HC1 <- subset(HC1, cells = HC1.RSC)
HC2 <- subset(HC2, cells = HC2.RSC)
HC3 <- subset(HC3, cells = HC3.RSC)
HC4 <- subset(HC4, cells = HC4.RSC)
HC5 <- subset(HC5, cells = HC5.RSC)
HC6 <- subset(HC6, cells = HC6.RSC)
HC7 <- subset(HC7, cells = HC7.RSC)

# Assign condition before integration
SOR1$condition <- "SOR"
SOR2$condition <- "SOR"
SOR3$condition <- "SOR"
SOR4$condition <- "SOR"
SOR5$condition <- "SOR"
SOR6$condition <- "SOR"
SOR7$condition <- "SOR"
HC1$condition <- "HC"
HC2$condition <- "HC"
HC3$condition <- "HC"
HC4$condition <- "HC"
HC5$condition <- "HC"
HC6$condition <- "HC"
HC7$condition <- "HC"

#===============================================================================
# Normalization and Integration 
#===============================================================================

# Make a list
obj.list <- list(SOR1, SOR2, SOR3, SOR4,
                 SOR5, SOR6, SOR7,
                 HC1, HC2, HC3, HC4, 
                 HC5, HC6, HC7)

# perform SCTransform normalization
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- SCTransform(obj.list[[i]], assay = 'Spatial')
}

# Alternative method
# obj.list <- lapply(X = obj.list, FUN = SCTransform, assay = "Spatial")


# Create RNA assay for all spatial objects
for (i in 1:length(obj.list)) {
  obj.list[[i]][["RNA"]] <- obj.list[[i]][["Spatial"]]
}


anchor_features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)

obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = anchor_features)

rsc.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                      anchor.features = anchor_features) # Takes time!

merged <- IntegrateData(anchorset = rsc.anchors, normalization.method = "SCT", k.weight = 40)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
# merged <- FindClusters(merged, resolution = 0.5)
merged <- FindClusters(merged, resolution = 0.1)

# Save the seurat object
#===============================================================================
saveRDS(merged, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/RSC_integrated_visium.rds")

merged <- readRDS("/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/RSC_integrated_visium.rds")

pdf(file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/SpatialDimplot_RSC.pdf", height = 10, width = 10)
SpatialDimPlot(merged, ncol = 4, pt.size.factor = 10, 
               label = TRUE, label.size = 3, image.alpha = )&NoLegend()
dev.off()

head(merged@meta.data)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform DE analysis (pseudobulking)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# pseudobulk the counts based on condition

# pseudo_rsc <- AggregateExpression(merged, assays = "RNA", return.seurat = T, 
#                                   group.by = c("orig.ident", "condition"))
# 
# # Check the count matrix
# head(pseudo_rsc@assays$RNA$counts)
# 
# pseudo_rsc@meta.data
# pseudo_rsc@assays
# ?FindMarkers
# 
# # each 'cell' is a condition-replicate pseudobulk profile
# Cells(pseudo_rsc)
# 
# # Perform DE testing (DESeq2)
# Idents(pseudo_rsc) <- "condition"
# pseudobulk.rsc <- FindMarkers(object = pseudo_rsc,
#                               ident.1 = "SOR", 
#                               ident.2 = "HC",
#                               test.use = "DESeq2")
# head(pseudobulk.rsc, n = 15)
# 
# signif <- pseudobulk.rsc %>%
#   dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2)%>%
#   tibble::rownames_to_column(var = "gene")
# nrow(signif) # 64 genes!
# 
# write.csv(pseudobulk.rsc, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/Pseudobulk_RSC.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform DE testing (wilcox.test)
head(merged@meta.data)

Idents(merged) <- "condition"
merged <- PrepSCTFindMarkers(merged)
pseudobulk.rsc <- FindMarkers(merged,
                              ident.1 = "SOR",
                              ident.2 = "HC",
                              min.pct=0.2, logfc.threshold=0, 
                              test.use = "wilcox",
                              assay = "SCT")

signif.visium <- pseudobulk.rsc %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2)%>%
  tibble::rownames_to_column(var = "gene")
nrow(signif.visium) # 77 genes!

pseudo_vis_up <- pseudobulk.rsc %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.2)
nrow(pseudo_vis_up)# 64

pseudo_vis_down <- pseudobulk.rsc %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.2)
nrow(pseudo_vis_down)# 13


write.csv(pseudobulk.rsc, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/Pseudobulk_RSC_eachdotasreplicate.csv")


# Make a Volcano Plot
#===============================================================================

data <- data.frame(gene = row.names(pseudobulk.rsc),
                   pval = -log10(pseudobulk.rsc$p_val_adj), 
                   lfc = pseudobulk.rsc$avg_log2FC)

data <- mutate(data, color = case_when(data$lfc > 0.2 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < -0.2 & data$pval > 1.3 ~ "Decreased",
                                       data$lfc >= -0.2 & data$lfc <= 0.2 & data$pval > 1.3 ~ "nonsignificant",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))+
  geom_point(size = 3, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#DE3163", 
                                Decreased = "#40E0D0",
                                nonsignificant = "gray90")) +
  theme_base() + # change overall theme
  theme(legend.position = "none") + # change the legend
  # xlim(-2,2)+
  # scale_y_continuous(breaks = seq(0, 350, by=50), limits=c(0,400))+
  theme(plot.title = element_text(size = 16),
        axis.text = element_text(size = 14))+
  geom_hline(yintercept = 1.3,
             colour = "gray",
             linetype="dashed")+
  geom_vline(xintercept = 0.2,
             colour = "gray",
             linetype = "dashed")+
  geom_vline(xintercept = -0.2,
             colour = "gray",
             linetype = "dashed")+
  ggtitle("Pseudobulk DGE")+
  xlab(expression(log[2]~"(Fold Change)"))+
  ylab(expression(-log[10]~"(Adj P Value)"))+
  geom_text_repel(data=data %>%
                    arrange(-lfc)%>%
                    filter(color != "nonsignificant")%>%
                    head(10), aes(label=gene),
                  size = 7,
                  box.padding = unit(.9, "lines"),hjust= 0.30,
                  segment.color = 'black',max.overlaps = Inf,
                  colour = 'black') +
  geom_text_repel(data=data %>%
                    arrange(lfc)%>%
                    filter(color != "nonsignificant")%>%
                    head(5), aes(label=gene),
                  size = 7,
                  box.padding = unit(.9, "lines"),hjust= 0.30,
                  segment.color = 'black',max.overlaps = Inf,
                  colour = 'black')+
  annotate("text", x = 1, y = 150, fontface = "bold",
           label = as.character(paste0("Up: ", nrow(pseudo_vis_up))),
           size = 6)+
  annotate("text", x = -0.6, y = 150, fontface = "bold",
           label = as.character(paste0("Down: ", nrow(pseudo_vis_down))),
           size = 6)+
  ggpubr::labs_pubr()


vol

# ggsave(filename = "/home/bbasu/hpchome/xenium_SOR/Analysis_new/Volcano_plot_visium_pseudobulk.pdf",
#        vol,
#        height = 8,
#        width = 10,
#        units = "in",
#        dpi = 600)

ggsave(filename = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/Volcano_plot_visium_pseudobulk.pdf",
       vol,
       height = 8,
       width = 10,
       units = "in",
       dpi = 600)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make a heatmap
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ComplexHeatmap)
library(circlize)

# Plot top 50 genes
data <- pseudobulk.rsc %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2)%>%
  dplyr::select(avg_log2FC)%>%
  head(50)

mat <- as.matrix(data)

my_palette <- colorRamp2(
  breaks = c(-1, 0, 2.5),  # Define the range and midpoint
  colors = c("#40E0D0", "white", "#DE3163")  # Colors for low, midpoint, and high
)

ht1 <- Heatmap(mat, cluster_columns = F, cluster_rows = F,
               width = unit(1.5, "cm"), 
               height = unit(35, "cm"),
               show_column_names = F,
               show_row_names = T,
               col = my_palette,
               name = "logFC",
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10,
                                 fontface = "bold")))

pdf(file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/heatmap_visium_pseudobulk.pdf",
    height = 15, width = 6)
ht1
dev.off()



#===============================================================================
# Perform GO enrichment 
#===============================================================================
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(tidyverse)

ego <- enrichGO(gene          = signif.visium$gene,
                OrgDb         = org.Mm.eg.db, # or Org.Hs.eg.db
                ont           = "MF",
                #one of “BP”, “MF”, “CC” or “ALL”
                pAdjustMethod = "BH",
                #one of “bonferroni”, “BH”, “BY”, “fdr”, “none”
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType = "SYMBOL",
                #“ENSEMBL”, “ENTREZID”, “SYMBOL”
                readable      = TRUE)
ego_result <- ego@result

# Make CNET plot
fc <- signif.visium %>%
  dplyr::select(gene, avg_log2FC)%>%
  deframe()

cnet_plot <- cnetplot(ego, circular = TRUE, colorEdge = TRUE, foldChange = fc,
                      showCategory = 10,
                      cex_category = 3.5,
                      cex_gene = 3.5,
                      cex_label_category = 3.5,
                      cex_label_gene = 3.5)+
  scale_color_gradient2(low = "#40E0D0", mid = "white", high = "#DE3163", 
                        midpoint = 0, 
                        breaks = c(-1, 0, 2),
                        limits = c(-1, 2),
                        labels = c(-1, 0, 2))+
  theme(legend.text = element_text(size = 12, face= "bold"),
        legend.title = element_text(size = 14, face = "bold"))+
  labs(color = "log2FC")

cnet_plot
write.xlsx(ego_result, file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/MF_enrichment_visium_pseudobulk.xlsx", 
           rowNames = F)

pdf(file = "/home/bbasu/LSS/lss_schatterj/rawdata_visium_RSC/pseudobulk_MF_CNET_plot.pdf", 
    height = 26, width = 35)
cnet_plot
dev.off()


# paletteer::scale_color_paletteer_c("ggthemes::Red-Blue Diverging", direction = -1,
#                                    limits = c(-2,2))
# 
# scale_color_gradient2(low = "#40E0D0", mid = "white", high = "#DE3163", midpoint = 0, breaks = breaks)+


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sankey Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
install.packages("ggalluvial")
library(ggalluvial)
# Load integrated Xenium object
merged_xenium <- readRDS("/home/bbasu/hpchome/xenium_SOR/RSC_integrated_xenium_new.rds")

xenium_markers <- read.xlsx("/home/bbasu/hpchome/xenium_SOR/Xenium mouse gene panel.xlsx")
xenium_50_gene_panel <- read.xlsx("/home/bbasu/hpchome/xenium_SOR/Updated Xenium 50 custom genes.xlsx")
xenium.gene <- cbind(xenium_markers$Genes, xenium_50_gene_panel$Gene)


# int.genes = intersect(signif.visium$gene, 
#                       xenium.gene) # 38 genes!

pseudo_vis_up <- pseudobulk.rsc %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.2)
nrow(pseudo_vis_up)# 64

pseudo_vis_down <- pseudobulk.rsc %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.2)
nrow(pseudo_vis_down)# 13

# up.int <- intersect(rownames(pseudo_vis_up), xenium.gene)#35 genes!
# down.int <- intersect(rownames(pseudo_vis_down), xenium.gene)#3 genes! # "Ttr"   "Zc3h6" "Cldn5"

head(merged_xenium@meta.data)
Idents(merged_xenium) <- merged_xenium$class
merged_xenium$class.condition <- paste(Idents(merged_xenium), merged_xenium$condition, sep = "_")

DEG <- list()
Idents(merged_xenium) <- merged_xenium$class
merged_xenium <- PrepSCTFindMarkers(merged_xenium)
for(i in levels(merged_xenium)){
  cluster = i
  message("Calculating DEG for ",cluster)
  Idents(merged_xenium) <- merged_xenium$class.condition
  cluster.marker <- FindMarkers(merged_xenium, ident.1 = paste0(cluster,"_SOR"), 
                                ident.2 = paste0(cluster, "_HC"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox",
                                assay = "SCT") 
  DEG[[cluster]] <- cluster.marker
}


int.exn <- cbind(intersect(rownames(pseudo_vis_up), 
                           rownames(DEG$ExN %>%
                                      dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.2))),
                 intersect(rownames(pseudo_vis_down), 
                           rownames(DEG$ExN %>%
                                      dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.2))))%>%
  paste(collapse = ", ")

int.InN <- cbind(intersect(rownames(pseudo_vis_up), 
                           rownames(DEG$InN %>%
                                      dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.2))),
                 intersect(rownames(pseudo_vis_down), 
                           rownames(DEG$InN %>%
                                      dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.2))))%>%
  paste(collapse = ", ")

int.astro <- cbind(intersect(rownames(pseudo_vis_up), 
                             rownames(DEG$Astro %>%
                                        dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.2))),
                   intersect(rownames(pseudo_vis_down), 
                             rownames(DEG$Astro %>%
                                        dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.2))))%>%
  paste(collapse = ", ")

int.oligo <- cbind(intersect(rownames(pseudo_vis_up), 
                             rownames(DEG$Oligo %>%
                                        dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.2))),
                   intersect(rownames(pseudo_vis_down), 
                             rownames(DEG$Oligo %>%
                                        dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.2))))%>%
  paste(collapse = ", ")

int.microglia <- cbind(intersect(rownames(pseudo_vis_up), 
                                 rownames(DEG$Microglia %>%
                                            dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.2))),
                       intersect(rownames(pseudo_vis_down), 
                                 rownames(DEG$Microglia %>%
                                            dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.2))))%>%
  paste(collapse = ", ")


clusters <- c(int.exn, int.InN, int.astro, int.oligo, int.microglia) # int.endo, int.VLMC, int.OPC
class <- c("ExN", "InN", "Astro", "Oligo", "Microglia") # "Endo", "VLMC", "OPC"


data <- data.frame(Visium = rep(int.genes, each = length(clusters))  ,
                   Xenium = rep(clusters, times = length(int.genes)),
                   Class = rep(class, times = length(int.genes)),
                   stringsAsFactors = FALSE)

data$Class <- factor(data$Class,
                     levels = c("ExN", "InN", "Astro", "Oligo", "Microglia"),
                     ordered = T)

# Add a score column: 1 if Gene is in the Cluster, otherwise 0
data$Score <- mapply(function(Visium, Xenium) {
  if (grepl(Visium, Xenium)) 1 else 0
}, data$Visium, data$Xenium)


plot <- ggplot(data = data,
               aes(axis1 = Visium, axis2 = Class,
                   y = Score)) +
  scale_x_discrete(limits = c("Visium", "Xenium")) +
  xlab("") +
  geom_alluvium(aes(fill = Xenium), width = 1/8) +
  geom_stratum(width = 1/8, alpha = 1, aes(color = Xenium))+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size = 3,discern=TRUE)+
  scale_fill_manual(values = c("#E76BF3", "#00B0F6", "#F8766D", "#A3A500","#00BF7D"))+
  scale_color_manual(values = c("#E76BF3", "#00B0F6", "#F8766D", "#A3A500","#00BF7D"))+
  theme_void()+
  theme(legend.position = "none")

plot
ggsave(filename = "/home/bbasu/hpchome/xenium_SOR/Analysis_new/Sankey_plot.svg",
       plot,
       height = 10,
       width = 8,
       units = "in",
       dpi = 600)
