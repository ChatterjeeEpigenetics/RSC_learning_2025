# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                          Xenium SOR Project                        -----
# -----                                                                    -----
# -----                           Chatterjee Lab                           -----
# -----                         University of Iowa                         -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
# Budhaditya Basu
library(Seurat)
library(ggplot2)
library(future)
library(dplyr)
library(scCustomize)
library(patchwork)
library(openxlsx)
library(ggthemes)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(UpSetR)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the xenium data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.dir <- "/home/bbasu/LSS/lss_schatterj/xenium_SOR/10283-SB/10283-SB/output-XETG00077__0040309__10283-SB-1__20240626__212859/"

xenium.obj <- LoadXenium(data.dir = data.dir, fov = "fov", assay = "Xenium")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
head(xenium.obj@meta.data)
# Save the xenium object
saveRDS(xenium.obj, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/xenium_object.rds")
# Load in data 
xenium.obj <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/xenium_object.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ROI based cell selection in xenium explorer v3 and import the cells for each biological replicates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sor1.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_SOR1_cells_stats.csv",
                           header = T, skip = 2)%>%
  pull(Cell.ID)

sor2.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_SOR2_cells_stats.csv",
                           header = T, skip = 2)%>%
  pull(Cell.ID)

sor3.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_SOR3_cells_stats.csv",
                           header = T, skip = 2)%>%
  pull(Cell.ID)

sor4.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_SOR4_cells_stats.csv",
                           header = T, skip = 2)%>%
  pull(Cell.ID)

hc1.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_HC1_cells_stats.csv",
                          header = T, skip = 2)%>%
  pull(Cell.ID)

hc2.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_HC2_cells_stats.csv",
                          header = T, skip = 2)%>%
  pull(Cell.ID)

hc3.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_HC3_cells_stats.csv",
                          header = T, skip = 2)%>%
  pull(Cell.ID)

hc4.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_HC4_cells_stats.csv",
                          header = T, skip = 2)%>%
  pull(Cell.ID)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subset the biological replicates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOR1 <- subset(xenium.obj, cells = sor1.RSC.cells)
SOR2 <- subset(xenium.obj, cells = sor2.RSC.cells)
SOR3 <- subset(xenium.obj, cells = sor3.RSC.cells)
SOR4 <- subset(xenium.obj, cells = sor4.RSC.cells)

HC1 <- subset(xenium.obj, cells = hc1.RSC.cells)
HC2 <- subset(xenium.obj, cells = hc2.RSC.cells)
HC3 <- subset(xenium.obj, cells = hc3.RSC.cells)
HC4 <- subset(xenium.obj, cells = hc4.RSC.cells)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Assign condition before integration
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOR1$orig.ident <- "SOR1"
SOR2$orig.ident <- "SOR2"
SOR3$orig.ident <- "SOR3"
SOR4$orig.ident <- "SOR4"
HC1$orig.ident <- "HC1"
HC2$orig.ident <- "HC2"
HC3$orig.ident <- "HC3"
HC4$orig.ident <- "HC4"

# Condition
SOR1$condition <- "SOR"
SOR2$condition <- "SOR"
SOR3$condition <- "SOR"
SOR4$condition <- "SOR"
HC1$condition <- "HC"
HC2$condition <- "HC"
HC3$condition <- "HC"
HC4$condition <- "HC"


head(SOR1@meta.data)
# Save as seurat object
saveRDS(SOR1, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR1_RSC.rds")
saveRDS(SOR2, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR2_RSC.rds")
saveRDS(SOR3, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR3_RSC.rds")
saveRDS(SOR4, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR4_RSC.rds")

saveRDS(HC1, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC1_RSC.rds")
saveRDS(HC2, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC2_RSC.rds")
saveRDS(HC3, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC3_RSC.rds")
saveRDS(HC4, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC4_RSC.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import the data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOR1 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR1_RSC.rds")
SOR2 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR2_RSC.rds")
SOR3 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR3_RSC.rds")
SOR4 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/SOR4_RSC.rds")

HC1 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC1_RSC.rds")
HC2 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC2_RSC.rds")
HC3 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC3_RSC.rds")
HC4 <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/HC4_RSC.rds")

head(SOR1@meta.data)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SCTransform normalization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_list <- list(SOR1, SOR2, SOR3, SOR4,
                 HC1, HC2, HC3, HC4)

# perform SCTransform normalization
for (i in 1:length(obj_list)) {
  obj_list[[i]] <- SCTransform(obj_list[[i]], assay = "Xenium")
}

# Create RNA assay for all spatial objects
for (i in 1:length(obj_list)) {
  obj_list[[i]][["RNA"]] <- obj_list[[i]][["Xenium"]]
}

anchor_features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)

obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = anchor_features)

rsc.anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                      anchor.features = anchor_features) # Takes time!

merged <- IntegrateData(anchorset = rsc.anchors, normalization.method = "SCT")
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.1)
# Prior to performing differential expression, we first run PrepSCTFindMarkers, which ensures that the fixed value is set properly. 
# Then we use FindMarkers(assay="SCT") to find differentially expressed genes. 

merged <- PrepSCTFindMarkers(merged)

# Save the seurat object
saveRDS(merged, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cluster annotation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Marker genes for all clusters
#===============================================================================
# Import the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium.rds")

cluster.markers <- FindAllMarkers(merged, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.20,
                                  test.use = "wilcox")

# Retrieve the top 10 marker genes per cluster
# Use whichever genes have the highest values under the AVG_LOG column
# top5 <- cluster.markers %>% group_by(cluster) %>%
#   dplyr::slice_max(get(grep("^avg_log", colnames(cluster.markers), value = TRUE)),
#                    n = 5)
# top10 markers
top10 <- cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Create the dot plot
pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Cluster_markers.pdf",
    height = 15, width = 25)
Seurat::DotPlot(merged, features = unique(top10$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                     size = 8, hjust = 1)) +
  Seurat::NoLegend()+
  ggpubr::labs_pubr()
dev.off()

#===============================================================================
# Listing cluster based markers
#===============================================================================
# Cluster 0: Oligo: Sox10, Opalin, Sgk1
# Cluster 1: Astro: Aqp4
# CLuster 2: L2/3: Gsg1l, Nwd2    
# Cluster 3: Endo: Cldn5, Sox17
# Cluster 4: L4: Kcnh5, Gfra2     
# Cluster 5: L6: Rprm, Trbc2
# Cluster 6: L4 RSP: Nell1
# Cluster 7: VLMC: Igf2, Aldh1a2
# Cluster 8: Sst: Gad1, Gad2, Sst ..........to be subclustered
# Cluster 9: Microglia: Cd53, Trem2
# Cluster 10: OPC: Pdgfra
# Cluster 11: L5: Vat1l
# Cluster 12: Vip: Vip .................... to be subclustered
# Cluster 13: NP: Vwc2l
# Cluster 14: L4 RSP: Cbln1
# Cluster 15: L6: Ccn2, Cplx3

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subclustering clusters 8 and 12 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
merged@graphs

# Suclustering 8
subclustered_8 <- FindSubCluster(merged, cluster=c("8"), subcluster.name = "subclusters_8", 
                                 graph.name = "integrated_snn", resolution = 0.1)
unique(subclustered_8$subclusters_8)

Idents(subclustered_8) <- subclustered_8$subclusters_8
my_levels <- c("0","1","2","3","4","5","6","7","8_0","8_1","9","10","11","12","13","14","15")
Idents(subclustered_8) <- factor(Idents(subclustered_8), levels= my_levels)

Stacked_VlnPlot(subclustered_8, 
                features = c("Opalin", "Aqp4","Gsg1l",
                             "Cldn5", "Gfra2",
                             "Rprm", "Nell1",
                             "Igf2", "Pvalb","Sst",
                             "Cd53", "Pdgfra", "Vat1l",
                             "Vip", "Vwc2l", "Cbln1", "Ccn2"),
                x_lab_rotate = 45,
                plot_legend = TRUE)


# Subclustering 12
Idents(subclustered_8) <- subclustered_8$subclusters_8

subclustered_12 <- FindSubCluster(subclustered_8, cluster=c("12"), subcluster.name = "subclusters_12", 
                                  graph.name = "integrated_snn", resolution = 0.1)
unique(subclustered_12$subclusters_12)
head(subclustered_12@meta.data)
Idents(subclustered_12) <- subclustered_12$subclusters_12
cluster.markers <- FindAllMarkers(subclustered_12, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.20,
                                  test.use = "wilcox")

Seurat::DotPlot(subclustered_12, features = unique(top10$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                     size = 8, hjust = 1)) +
  Seurat::NoLegend()+
  ggpubr::labs_pubr()

# Save the seurat object
saveRDS(subclustered_12, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")

# Add cluster information to the metadata
# Annotate cluster
# Cluster 0: Oligo: Sox10, Opalin, Sgk1
# Cluster 1: Astro: Aqp4
# CLuster 2: L2/3: Gsg1l, Nwd2            # Glut
# Cluster 3: Endo: Cldn5, Sox17
# Cluster 4: L4: Kcnh5, Gfra2             # Glut
# Cluster 5: L6: Rprm, Trbc2              # Glut
# Cluster 6: L2/3 RSP: Rnf152,Slc17a6     # Glut
# Cluster 7: VLMC: Igf2, Aldh1a2
# Cluster 8_0: Pvalb: Pvalb 
# Cluster 8_1: Sst: Sst 
# Cluster 9: Microglia: Cd53, Trem2
# Cluster 10: OPC: Pdgfra
# Cluster 11: L5: Vat1l
# Cluster 12_0: Vip: Vip
# Cluster 12_1: Lamp5: Pde11a, Lamp5
# Cluster 12_2: Sncg: Sncg
# Cluster 12_3: Vip: Penk
# Cluster 13: NP SUB: Vwc2l
# Cluster 14: L4 RSP: Cbln1, Nell1
# Cluster 15: L6: Ccn2, Cplx3

# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")
head(merged@meta.data)

# Annotate Clusters
Idents(merged) <- merged$subclusters_12

celltype <- rep(NA, length = ncol(merged))
celltype[which(Idents(merged) %in% c(0))] <- 'Oligo' 
celltype[which(Idents(merged) %in% c(1))] <- 'Astro'
celltype[which(Idents(merged) %in% c(2))] <- 'L2_3'
celltype[which(Idents(merged) %in% c(3))] <- 'Endo'
celltype[which(Idents(merged) %in% c(4))] <- 'L4'
celltype[which(Idents(merged) %in% c(6))] <- 'L2_3 RSP'
celltype[which(Idents(merged) %in% c(14))] <- 'L4 RSP'
celltype[which(Idents(merged) %in% c(5, 15))] <- 'L6'
celltype[which(Idents(merged) %in% c(7))] <- 'VLMC'
celltype[which(Idents(merged) %in% c("8_0"))] <- 'Pvalb'
celltype[which(Idents(merged) %in% c("8_1"))] <- 'Sst'
celltype[which(Idents(merged) %in% c(9))] <- 'Microglia'
celltype[which(Idents(merged) %in% c(10))] <- 'OPC'
celltype[which(Idents(merged) %in% c(11))] <- 'L5'
celltype[which(Idents(merged) %in% c("12_0", "12_3"))] <- 'Vip'
celltype[which(Idents(merged) %in% c("12_1"))] <- 'Lamp5'
celltype[which(Idents(merged) %in% c("12_2"))] <- 'Sncg'
celltype[which(Idents(merged) %in% c(13))] <- 'NP SUB'


celltype <- factor(celltype, 
                   levels = c('L2_3', 'L4', 'L5', 'L6', 'NP SUB', 'L2_3 RSP', 'L4 RSP',
                              'Pvalb', 'Sst', 'Lamp5', 'Sncg', 'Vip',
                              'Oligo', 'Astro', 'Microglia', 'OPC',
                              'Endo', 'VLMC'), 
                   ordered = T)
head(merged@meta.data)
merged$celltype <- celltype

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Save the seurat object
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveRDS(merged, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")

colors_celltype <- c("#1F77B4FF", "#E377C2FF","#D62728FF", "#6677d4", "#720ba9","#FF7F0EFF","#2CA02CFF",
                     "#17BECFFF", "#b7950b", "#98DF8AFF", "#FF9896FF", "#411eee",
                     "#da5e0e", "#9EDAE5FF","#6495ED", "#FFBB78FF", "#8C564BFF", "#C49C94FF")

Idents(merged) <- merged$celltype

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Celltype_Dimplot_new.pdf", height = 6, width = 8)
DimPlot_scCustom(merged, colors_use = colors_celltype)
dev.off()


# Save representative ImageDimplots
pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Celltype_Image_Dimplot_new.pdf", height = 12, width = 15)
ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.3", "fov.8"), 
             dark.background = FALSE, split.by = "orig.ident",
             border.size = NA)
dev.off()


DefaultAssay(merged) <- "Xenium"
head(merged@meta.data)
pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Stacked_VlnPlot_markers.pdf",
    height = 12, width = 12)
Stacked_VlnPlot(merged, 
                features = c("Gsg1l", "Gfra2","Vat1l","Rprm",
                             "Vwc2l", "Rnf152", "Cbln1", "Pvalb","Sst","Pde11a", 
                             "Sncg","Vip", "Opalin", "Aqp4",
                             "Cd53","Pdgfra","Cldn5","Igf2"),
                x_lab_rotate = 45,
                layer = "count",
                plot_legend = TRUE,
                colors_use = colors_celltype)
dev.off()



Idents(merged) <- merged$celltype
colors_celltype <- paletteer::paletteer_d("ggthemes::Classic_20")

p1 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p2 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.2"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p3 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.3"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p4 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.4"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p5 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.5"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p6 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.6"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p7 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.7"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p8 <- ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.8"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)


pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Image_Dimplot_all_replicate.pdf", height = 15, width = 20)
p1+p2+p3+p4+p5+p6+p7+p8
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Image_Dimplot_all_replicate_original_color_scheme.pdf", height = 15, width = 20)
p1+p2+p3+p4+p5+p6+p7+p8
dev.off()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Separate visualization of celltypes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Show only excitatory layers
#===============================================================================

colors_celltype <-  c("#1F77B4FF", "#E377C2FF","#D62728FF", "#6677d4", "#720ba9","#FF7F0EFF","#2CA02CFF",
                      "#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", 
                      "#C7C7C7FF","#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF")

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Excitatory_Image_Dimplot.pdf", height = 10, width = 12)
ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.3", "fov.8"), 
             dark.background = FALSE, split.by = "orig.ident",
             border.size = NA)
dev.off()

# Show only Inhibitory neurons
#===============================================================================
colors_celltype <-  c("#C7C7C7FF","#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", 
                      "#17BECFFF", "#b7950b", "#98DF8AFF", "#FF9896FF", "#411eee", 
                      "#C7C7C7FF","#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF")

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Inhibitory_Image_Dimplot.pdf", height = 10, width = 12)
ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.3", "fov.8"), 
             dark.background = FALSE, split.by = "orig.ident",
             border.size = NA, cells = )
dev.off()

# Show only Glial cells
#===============================================================================
colors_celltype <-  c("#C7C7C7FF","#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", 
                      "#C7C7C7FF","#C7C7C7FF","#C7C7C7FF", "#C7C7C7FF","#C7C7C7FF", 
                      "#da5e0e", "#9EDAE5FF","#6495ED", "#FFBB78FF","#C7C7C7FF", "#C7C7C7FF")

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Glial_Image_Dimplot.pdf", height = 10, width = 12)
ImageDimPlot(merged, cols = colors_celltype, size = 2, fov = c("fov.3", "fov.8"), 
             dark.background = FALSE, split.by = "orig.ident",
             border.size = NA, cells = )
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make Proportion plot for Celltype
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")

head(merged@meta.data)

Idents(merged) <- merged$celltype

table(Idents(merged), merged$condition)

df <- as.data.frame(prop.table(table(Idents(merged), 
                                     merged$condition), margin = 2))


plot <- ggbarplot(df, x = "Var2", y = "Freq", fill = "Var1", 
                  palette = alpha(colors_celltype, 0.9), 
                  xlab = "", 
                  ylab = "Proportion", label = F)+
  labs_pubr()+
  theme(legend.position = "right",
        legend.title = element_blank())+
  rotate_x_text(45)

plot
ggsave(filename = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/RSC_WT_Celltype_Proportion.pdf",
       plot,
       height = 6,
       width = 6,
       units = "in",
       dpi = 600)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform DEG across celltype
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Idents(merged) <- merged$celltype
merged$celltype.condition <- paste(Idents(merged), merged$condition, sep = "_")

head(merged@meta.data)

DEG <- list()
Idents(merged) <- merged$celltype
merged <- PrepSCTFindMarkers(merged)
for(i in levels(merged)){
  cluster = i
  message("Calculating DEG for ",cluster)
  Idents(merged) <- merged$celltype.condition
  cluster.marker <- FindMarkers(merged, ident.1 = paste0(cluster,"_SOR"), 
                                ident.2 = paste0(cluster, "_HC"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox",
                                assay = "SCT") 
  DEG[[cluster]] <- cluster.marker
}

saveRDS(DEG, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_celltype.rds")
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_celltype.rds")
# Store in spreadsheet
write.xlsx(DEG, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_celltype_new.xlsx", 
           rowNames = T)

DEG_upreg <- list()
for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(avg_log2FC > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  DEG_upreg[[cluster]] <- gene.list$gene
}

DEG_upreg_ExN <- DEG_upreg[c("NP SUB","L6","L5","L4 RSP", "L4","L2_3 RSP","L2_3")]
DEG_upreg_InN <- DEG_upreg[c("Pvalb", "Sst", "Lamp5", "Vip")]
#===============================================================================
# Upset Plot (Upregulated)
#===============================================================================
pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/InN_celltypes_Upset.pdf",
    height = 8, width = 14, onefile = FALSE)
upset_plot <- print(upset(fromList(DEG_upreg_InN),
                          sets.bar.color = c("#411eee","#b7950b",  
                                             "#17BECFFF","#98DF8AFF"),
                          order.by = "freq",
                          text.scale = 2,
                          mainbar.y.label = "Unique and overlapping genes",
                          sets.x.label = "# Genes upregulated by learning",
                          point.size = 4,
                          set_size.show = TRUE))
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/ExN_celltypes_Upset.pdf",
    height = 8, width = 14, onefile = FALSE)
upset_plot <- print(upset(fromList(DEG_upreg_ExN),
                          sets.bar.color = c("#720ba9","#6677d4","#D62728FF","#2CA02CFF",
                                             "#E377C2FF","#FF7F0EFF","#1F77B4FF"),
                          nsets = 7,
                          order.by = "freq",
                          text.scale = 2,
                          sets = c("NP SUB","L6","L5","L4 RSP", "L4","L2_3 RSP","L2_3"),
                          mainbar.y.label = "Unique and overlapping genes",
                          sets.x.label = "Upregulated by learning",
                          point.size = 4,
                          set_size.show = TRUE,
                          keep.order = T))
dev.off()

head(merged@meta.data)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Broad cluster annotation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ExN: Slc17a6, Slc17a7
# InN: Gad1, Gad2
# Non-neuro
Idents(merged) <- merged$celltype

Stacked_VlnPlot(merged, 
                features = c("Slc17a6", "Slc17a7", 
                             "Gad1", "Gad2"),
                x_lab_rotate = 45,
                plot_legend = TRUE)


# Add Class information to the metadata
# Annotate Class
Idents(merged) <- merged$celltype
class <- rep(NA, length = ncol(merged))

class[which(Idents(merged) %in% c('L2_3'))] <- 'ExN'
class[which(Idents(merged) %in% c('L4'))] <- 'ExN'
class[which(Idents(merged) %in% c('L5'))] <- 'ExN'
class[which(Idents(merged) %in% c('L6'))] <- 'ExN'
class[which(Idents(merged) %in% c('NP SUB'))] <- 'ExN'
class[which(Idents(merged) %in% c('L2_3 RSP'))] <- 'ExN'
class[which(Idents(merged) %in% c('L4 RSP'))] <- 'ExN'

class[which(Idents(merged) %in% c('Pvalb'))] <- 'InN'
class[which(Idents(merged) %in% c('Sst'))] <- 'InN'
class[which(Idents(merged) %in% c('Vip'))] <- 'InN'
class[which(Idents(merged) %in% c('Lamp5'))] <- 'InN'
class[which(Idents(merged) %in% c('Sncg'))] <- 'InN'

class[which(Idents(merged) %in% c('Oligo'))] <- 'Oligo' 
class[which(Idents(merged) %in% c('Astro'))] <- 'Astro'
class[which(Idents(merged) %in% c('Microglia'))] <- 'Microglia'
class[which(Idents(merged) %in% c('OPC'))] <- 'OPC'
class[which(Idents(merged) %in% c('Endo'))] <- 'Endo'
class[which(Idents(merged) %in% c('VLMC'))] <- 'VLMC'

class <- factor(class, 
                levels = c("ExN","InN","Astro","Oligo","Microglia","Endo","VLMC","OPC"), 
                ordered = T)

merged$class <- class
head(merged@meta.data)

unique(merged$class)
Idents(merged) <- merged$class

# colors_class <- c("#E76BF3", "#00B0F6", "#00BF7D", "#A3A500", "#F8766D", "#CAB2D6", "#B2DF8A", "#A6CEE3")

colors_class <- c("ExN"="#E76BF3",
                  "InN"="#00B0F6",
                  "Astro"="#00BF7D",
                  "Oligo"="#A3A500",
                  "Microglia"="#F8766D",
                  "Endo"="#CAB2D6",
                  "VLMC"="#B2DF8A",
                  "OPC"="#A6CEE3")

DefaultAssay(merged) <- "Xenium"

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Stacked_VlnPlot_markers_Broader_class_WT_RSC.pdf",
    height = 12, width = 12)
Stacked_VlnPlot(merged, 
                features = c("Slc17a7",
                             "Gad2",
                             "Aqp4",
                             "Opalin",
                             "Cd53",
                             "Cldn5",
                             "Igf2",
                             "Pdgfra"),
                x_lab_rotate = 45,
                plot_legend = TRUE,
                layer = "count",
                colors_use = colors_class)
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Broad_class_annotation_Dimplot.pdf", height = 6, width = 8)
DimPlot_scCustom(merged, colors_use = colors_class, figure_plot = TRUE, label = F)
dev.off()


# Save representative ImageDimplots
pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Broad_Class_Image_Dimplot.pdf", height = 10, width = 12)
ImageDimPlot(merged, cols = colors_class, size = 3, fov = c("fov.3", "fov.8"), 
             dark.background = FALSE, split.by = "orig.ident",
             border.size = NA)
dev.off()

# Save the seurat object
saveRDS(merged, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make Proportion plot for broader class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")

Idents(merged) <- merged$class

table(Idents(merged), merged$condition)

df <- as.data.frame(prop.table(table(Idents(merged), 
                                     merged$condition), margin = 2))


plot <- ggbarplot(df, x = "Var2", y = "Freq", fill = "Var1", 
                  palette = alpha(colors_class, 0.9), 
                  xlab = "", 
                  ylab = "Proportion", label = F)+
  labs_pubr()+
  theme(legend.position = "right",
        legend.title = element_blank())+
  rotate_x_text(45)

plot
ggsave(filename = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Broad_Class_Proportion.pdf",
       plot,
       height = 6,
       width = 6,
       units = "in",
       dpi = 600)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform DEG across different class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")

Idents(merged) <- merged$class
merged$class.condition <- paste(Idents(merged), merged$condition, sep = "_")

head(merged@meta.data)

DEG <- list()
Idents(merged) <- merged$class
merged <- PrepSCTFindMarkers(merged)
for(i in levels(merged)){
  cluster = i
  message("Calculating DEG for ",cluster)
  Idents(merged) <- merged$class.condition
  cluster.marker <- FindMarkers(merged, ident.1 = paste0(cluster,"_SOR"), 
                                ident.2 = paste0(cluster, "_HC"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox",
                                assay = "SCT") 
  DEG[[cluster]] <- cluster.marker
}

saveRDS(DEG, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_major_class.rds")
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_major_class.rds")
# Store in spreadsheet
write.xlsx(DEG, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_major_class.xlsx", 
           rowNames = T)

# Store significant DEGs across celltype
DEG_sig <- list()
for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  DEG_sig[[cluster]] <- gene.list
}

# Store in spreadsheet
write.xlsx(DEG_sig, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Significant_DEG_major_class.xlsx", 
           rowNames = F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform GO enrichment 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)


go_results <- list()
dot_plots <- list()
bar_plots <- list()

for (cluster in names(DEG_sig)) {
  genes <- DEG_sig[[cluster]]$gene
  message("Analysing GO:MF enrichment for ", cluster)
  tryCatch(ego <- enrichGO(gene          = genes,
                           OrgDb         = org.Mm.eg.db, # or Org.Mm.eg.db
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           keyType = "SYMBOL",
                           readable      = TRUE), error=function(e) NULL)
  
  tryCatch(go_results[[cluster]] <- ego@result, error=function(e) NULL)
  message("Plotting for ", cluster)
  dot_plots[[cluster]] <- dotplot(ego, showCategory=5, font.size = 14)+
    labs(title = paste0("GO:MF enrichment for ", cluster))
  bar_plots[[cluster]] <- barplot(ego, showCategory=5, font.size = 14)+labs(title = paste0("GO:MF enrichment for ", cluster))
}

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/GO-MF_plot_ExN_InN_Astro.pdf", height = 10, width = 20)
dot_plots$ExN+dot_plots$InN+dot_plots$Astro
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/GO-MF_plot_Oligo_microglia.pdf", height = 10, width = 20)
dot_plots$Oligo+dot_plots$Microglia+dot_plots$Astro
dev.off()

# Store in spreadsheet
write.xlsx(go_results, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/GO_MF_enrichment_major_class.xlsx", 
           rowNames = F)

#===============================================================================
# Make Volcano Plot
#===============================================================================
colors_class <- c("#E76BF3", "#00B0F6", "#00BF7D", "#A3A500", "#F8766D", "#CAB2D6", "#B2DF8A", "#A6CEE3")
plots <- list()
for(i in names(DEG)){
  cluster = i
  message("Plotting for ",cluster)
  cluster_index = which(names(DEG) == cluster)
  #Data wrangling
  cluster.marker <- DEG[[cluster]]
  
  data <- data.frame(gene = row.names(cluster.marker),
                     pval = -log10(cluster.marker$p_val_adj+2.225074e-308), #Add a tiny value to  avoid inf 2.225074e-308
                     lfc = cluster.marker$avg_log2FC)
  
  data <- mutate(data, color = case_when(data$lfc > 0.2 & data$pval > 1.3 ~ "Increased",
                                         data$lfc < -0.2 & data$pval > 1.3 ~ "Decreased",
                                         data$lfc >= -0.2 & data$lfc <= 0.2 & data$pval > 1.3 ~ "nonsignificant",
                                         data$pval < 1.3 ~ "nonsignificant"))
  
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color))+
    geom_point(size = 3, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = colors_class[cluster_index], 
                                  Decreased = colors_class[cluster_index],
                                  nonsignificant = "gray90")) +
    theme_base() + # change overall theme
    theme(legend.position = "none") + # change the legend
    # xlim(-2,2)+
    scale_y_continuous(breaks = seq(0, 350, by=50), limits=c(0,400))+
    theme(plot.title = element_text(size = 16),
          axis.text = element_text(size = 14))+
    geom_hline(yintercept = 1.3,
               colour = "gray",
               linetype="dashed")+
    geom_vline(xintercept = 0,
               colour = "gray",
               linetype = "dashed")+
    ggtitle(cluster)+
    xlab(expression(log[2]~"(Fold Change)"))+
    ylab(expression(-log[10]~"(Adj P Value)"))+
    geom_text_repel(data=data %>%
                      arrange(-pval)%>%
                      filter(color != "nonsignificant")%>%
                      head(10), aes(label=gene),
                    size = 7,
                    box.padding = unit(.9, "lines"),hjust= 0.30,
                    segment.color = 'black',max.overlaps = Inf,
                    colour = 'black') + ggpubr::labs_pubr()
  plots[[cluster]] <- vol
}


pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Volcano_Major_class.pdf", height = 10, width = 12)
plots$ExN+plots$InN+plots$Astro+plots$Oligo
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Volcano_Major_class_new.pdf", height = 10, width = 12)
plots$ExN+plots$InN+plots$Astro+plots$Oligo
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Volcano_Microglia_new.pdf", height = 7, width = 8)
plots$Microglia
dev.off()


pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Volcano_ExN.pdf", height = 6, width = 6)
plots$ExN
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Volcano_Oligo.pdf", height = 6, width = 6)
plots$Oligo
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Consistency across replicates
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/RSC_integrated_xenium_new.rds")

head(merged@meta.data)
pseudobulk_RSC_xenium <- AggregateExpression(merged, assays = "RNA",
                                             group.by = c("orig.ident", "condition"))
head(pseudobulk_RSC_xenium)

avg_expr <- pseudobulk_RSC_xenium$RNA
avg_expr <- as.matrix(avg_expr)
cor_matrix <- cor(avg_expr)

head(cor_matrix)

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Replicate_Xeniumdata_RSC_consistency.pdf",
    height = 8, width = 10)
pheatmap::pheatmap(cor_matrix, main = "RSC Gene Expression Correlation Across Replicates")
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# genes detected in <10% cells per replicate
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# proportion of cell types in each biological replicate
head(merged@meta.data)

Idents(merged) <- merged$class

prop.table(table(Idents(merged), merged$orig.ident), margin = 2)

# Get replicate info
replicates <- unique(merged$orig.ident)

# Initialize a list to store genes detected in <10% cells per replicate
low_detected_genes <- list()

DefaultAssay(merged) <- "RNA"
# rep = "SOR1"
# 
for (rep in replicates) {
  # Subset Seurat object by replicate
  rep_obj <- subset(merged, subset = orig.ident == rep)
  
  # expression matrix
  mat <- GetAssayData(rep_obj, layer = "counts")
  
  # fraction of cells where gene is detected (count > 0)
  detection_frac <- Matrix::rowSums(mat > 0) / ncol(mat)
  
  # Genes detected in less than 10% of cells
  low_detected_genes[[rep]] <- names(detection_frac[detection_frac < 0.1])
}

sor_low_detected_genes <- low_detected_genes[c("SOR1", "SOR2", "SOR3", "SOR4")]
hc_low_detected_genes <- low_detected_genes[c("HC1", "HC2", "HC3", "HC4")]

# Find genes that are lowly detected in ALL replicates
genes_low_all_reps_sor <- Reduce(intersect, sor_low_detected_genes)
genes_low_all_reps_hc <- Reduce(intersect, hc_low_detected_genes)

# check how many genes
length(genes_low_all_reps_sor)
length(genes_low_all_reps_hc)

# Combine the genes
genes_low_all_reps <- union(genes_low_all_reps_hc, genes_low_all_reps_sor)


# check if these genes were present in the significant DEG list
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_major_class.rds")
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_celltype.rds")
# Store significant DEGs across celltype
DEG_sig <- list()
for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  DEG_sig[[cluster]] <- gene.list
}


# lowly-expressed DEGs per cluster
low_expr_deg_by_cluster <- list()

# Loop over clusters
for (clust in names(DEG_sig)) {
  
  low_expr_genes <- genes_low_all_reps
  deg_genes <- DEG_sig[[clust]]
  
  # Find overlap
  matched_genes <- intersect(low_expr_genes, deg_genes$gene)
  
  # Store
  low_expr_deg_by_cluster[[clust]] <- matched_genes
  
}

# No gene was found in the DEG list of major cell type which was detected in less than 10% of cells of each replicate !!
# No gene was found in the DEG list of RSC neuronal sub type which was detected in less than 10% of cells of each replicate!!

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Fos positive versus Fos Negative 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DefaultAssay(merged) <- "Xenium"

DefaultAssay(merged) <- "RNA"
fos_exp <- FetchData(merged, vars = c("Fos"))
fos_exp$sample <- merged@meta.data$orig.ident
fos_exp$condition <- merged@meta.data$condition

fos_SOR <- fos_exp[fos_exp$condition == "SOR", ]

# Get cell barcodes where Fos > 0
cells_SOR_fos_positive <- rownames(fos_SOR)[fos_SOR$Fos > 0]

merged$Fos_exp <- ifelse(colnames(merged) %in% cells_SOR_fos_positive, "FosPos","FosNeg")

merged_SOR <- subset(merged, subset = condition == "SOR")

head(merged_SOR@meta.data)
unique(merged_SOR$condition)

# Calculate the number of fos_pos cells in each class
table(merged_SOR$class, merged_SOR$Fos_exp)
df <- data.frame(table(merged_SOR$class, merged_SOR$Fos_exp, merged_SOR$orig.ident))%>%
  tidyr::spread(key = "Var1", value = "Freq")

write.csv(df, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/fos_table_across_majorClass_WTSOR.csv",
          row.names = FALSE)


Idents(merged_SOR) <- merged_SOR$class
merged_SOR$class.fos <- paste(merged_SOR$class, merged_SOR$Fos_exp, sep = "_")
unique(merged_SOR$class.fos)

#Perform DGE across major class between FosPos vs FosNeg
#===============================================================================

DEG <- list()
Idents(merged_SOR) <- merged_SOR$class
for(i in levels(merged_SOR)){
  cluster = i
  message("Calculating DEG for ",cluster)
  Idents(merged_SOR) <- merged_SOR$class.fos
  cluster.marker <- FindMarkers(merged_SOR, ident.1 = paste0(cluster,"_FosPos"), 
                                ident.2 = paste0(cluster, "_FosNeg"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox",
                                assay = "SCT") 
  DEG[[cluster]] <- cluster.marker
}

saveRDS(DEG, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_FosPos_SOR.rds")
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_FosPos_SOR.rds")
# Store in spreadsheet
write.xlsx(DEG, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/DEG_FosPos_wtSOR.xlsx", 
           rowNames = T)

# Store significant DEGs across celltype
DEG_sig <- list()
for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  DEG_sig[[cluster]] <- gene.list
}

# Store in spreadsheet
write.xlsx(DEG_sig, file = "/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/Significant_DEG_FosPos_wtSOR.xlsx", 
           rowNames = F)


#===============================================================================
# Make Volcano Plot
#===============================================================================
library(ggbreak)
colors_class <- c("#E76BF3", "#00B0F6", "#00BF7D", "#A3A500", "#F8766D", "#CAB2D6", "#B2DF8A", "#A6CEE3")
plots <- list()
for(i in names(DEG)){
  cluster = i
  message("Plotting for ",cluster)
  cluster_index = which(names(DEG) == cluster)
  #Data wrangling
  cluster.marker <- DEG[[cluster]]
  
  data <- data.frame(gene = row.names(cluster.marker),
                     pval = -log10(cluster.marker$p_val_adj+2.225074e-308), #Add a tiny value to  avoid inf 2.225074e-308
                     lfc = cluster.marker$avg_log2FC)
  
  data <- mutate(data, color = case_when(data$lfc > 0.2 & data$pval > 1.3 ~ "Increased",
                                         data$lfc < -0.2 & data$pval > 1.3 ~ "Decreased",
                                         data$lfc >= -0.2 & data$lfc <= 0.2 & data$pval > 1.3 ~ "nonsignificant",
                                         data$pval < 1.3 ~ "nonsignificant"))
  
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color))+
    geom_point(size = 3, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = colors_class[cluster_index], 
                                  Decreased = colors_class[cluster_index],
                                  nonsignificant = "gray90")) +
    theme_classic() + # change overall theme
    theme(legend.position = "none") + # change the legend
    # xlim(-2,2)+
    scale_y_continuous(breaks = seq(0, 350, by=50), limits=c(0,400))+
    scale_x_break(breaks = c(1, 11), scales = 0.1)+ #ExN:c(1,11), InN:c(1.2, 8.9) scales = 0.1
    theme(plot.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    geom_hline(yintercept = 1.3,
               colour = "grey",
               linetype="dashed")+
    geom_vline(xintercept = 0,
               colour = "red2",
               linetype = "dashed")+
    geom_vline(xintercept = 0.2,
               colour = "grey",
               linetype = "dashed")+
    geom_vline(xintercept = -0.2,
               colour = "grey",
               linetype = "dashed")+
    ggtitle(cluster)+
    xlab(expression(log[2]~"(Fold Change)"))+
    ylab(expression(-log[10]~"(Adj P Value)"))+
    geom_text_repel(data=data %>%
                      arrange(-pval)%>%
                      filter(color != "nonsignificant" & lfc > 0)%>%
                      head(7), aes(label=gene),
                    size = 7,
                    box.padding = unit(0.9, "lines"), hjust= 0.30,
                    segment.color = 'black',max.overlaps = Inf,
                    colour = 'black') +
    geom_text_repel(data=data %>%
                      arrange(-pval)%>%
                      filter(color != "nonsignificant" & lfc < 0)%>%
                      head(5), aes(label=gene),
                    size = 7,
                    box.padding = unit(0.9, "lines"), hjust= 0.30,
                    segment.color = 'black',max.overlaps = Inf,
                    colour = 'black') +
    # geom_text_repel(data=data %>%
    #                   arrange(-pval)%>%
    #                   filter(color != "nonsignificant" & lfc < 0)%>%
    #                   head(3), aes(label=gene),
    #                 size = 7,
    #                 box.padding = unit(1.2, "lines"), hjust= 0.30,
    #                 segment.color = 'black',max.overlaps = Inf,
    #                 colour = 'black') + 
    ggpubr::labs_pubr()
  
  plots[[cluster]] <- vol
}

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/WT_SOR_FosPos_vs_FosNeg_volcano_ExN.pdf",
    height = 8, width = 10)
plots$ExN
dev.off()

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/WT_SOR_FosPos_vs_FosNeg_volcano_InN.pdf",
    height = 8, width = 12)
plots$InN
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Show Fos expression
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
head(merged@meta.data)
Assays(merged) 
DefaultAssay(merged) <- "Xenium" 
Idents(merged) <- merged$class

# Rename the FOV images
names(merged@images) <- c("SOR1", "SOR2", "SOR3", "SOR4",
                          "HC1", "HC2", "HC3", "HC4")

pdf("/home/bbasu/LSS/lss_schatterj/xenium_SOR/RSC/Analysis_new/RSC_WT_SOR_Fos_exp.pdf", height = 8, width = 12)
ImageFeaturePlot(merged, fov = c("SOR3", "HC4"),
                 features = c("Fos"),
                 cells = WhichCells(merged, idents = "ExN"),
                 size =3, cols = c("gray95", "firebrick2"), min.cutoff = "q10", max.cutoff = "q90",
                 scale = "all",
                 border.size = NA, dark.background = F)
dev.off()

