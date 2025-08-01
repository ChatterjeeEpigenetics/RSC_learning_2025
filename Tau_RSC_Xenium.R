# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                          Tau Xenium Project                        -----
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the xenium data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.dir <- "/home/bbasu/LSS/lss_schatterj/Tau_xenium/20250127__192401__12664-SB_12712-EP/output-XETG00077__0040312__12664-SB-1_ROI_A__20250127__192422/"

library(arrow)
transcripts <- read_parquet(file.path(data.dir, "transcripts.parquet"))
write.csv(transcripts, gzfile(file.path(data.dir, "transcripts.csv.gz")), row.names = FALSE)

#===============================================================================
xenium.obj <- LoadXenium(data.dir = data.dir, fov = "fov", assay = "Xenium")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
head(xenium.obj@meta.data)
dir.create("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/")

# Save the xenium object
saveRDS(xenium.obj, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/tau_xenium_object.rds")
# Load in data 
xenium.obj <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/tau_xenium_object.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ROI based RSC selection in xenium explorer v3 and import the cells for each biological replicates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau1.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau1_RSC_cells_stats.csv",
                           header = T, skip = 2)%>%
  pull(Cell.ID)

tau2.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau2_RSC_cells_stats.csv",
                           header = T, skip = 2)%>%
  pull(Cell.ID)

tau3.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau3_RSC_cells_stats.csv",
                           header = T, skip = 2)%>%
  pull(Cell.ID)

Cntrl1.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/C1_RSC_cells_stats.csv",
                             header = T, skip = 2)%>%
  pull(Cell.ID)

Cntrl2.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/C2_RSC_cells_stats.csv",
                             header = T, skip = 2)%>%
  pull(Cell.ID)

Cntrl3.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/C3_RSC_cells_stats.csv",
                             header = T, skip = 2)%>%
  pull(Cell.ID)

Cntrl4.RSC.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/C4_RSC_cells_stats.csv",
                             header = T, skip = 2)%>%
  pull(Cell.ID)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subset the biological replicates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau1 <- subset(xenium.obj, cells = tau1.RSC.cells)
Tau2 <- subset(xenium.obj, cells = tau2.RSC.cells)
Tau3 <- subset(xenium.obj, cells = tau3.RSC.cells)

Cntrl1 <- subset(xenium.obj, cells = Cntrl1.RSC.cells)
Cntrl2 <- subset(xenium.obj, cells = Cntrl2.RSC.cells)
Cntrl3 <- subset(xenium.obj, cells = Cntrl3.RSC.cells)
Cntrl4 <- subset(xenium.obj, cells = Cntrl4.RSC.cells)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Assign condition before integration
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau1$orig.ident <- "Tau1"
Tau2$orig.ident <- "Tau2"
Tau3$orig.ident <- "Tau3"

Cntrl1$orig.ident <- "Control1"
Cntrl2$orig.ident <- "Control2"
Cntrl3$orig.ident <- "Control3"
Cntrl4$orig.ident <- "Control4"

#===============================================================================
# Condition
#===============================================================================
Tau1$condition <- "Tau"
Tau2$condition <- "Tau"
Tau3$condition <- "Tau"

Cntrl1$condition <- "Control"
Cntrl2$condition <- "Control"
Cntrl3$condition <- "Control"
Cntrl4$condition <- "Control"

head(Tau1@meta.data)
# Save as seurat object
saveRDS(Tau1, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau1_RSC.rds")
saveRDS(Tau2, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau2_RSC.rds")
saveRDS(Tau3, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau3_RSC.rds")

saveRDS(Cntrl1, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Control1_RSC.rds")
saveRDS(Cntrl2, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Control2_RSC.rds")
saveRDS(Cntrl3, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Control3_RSC.rds")
saveRDS(Cntrl4, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Control4_RSC.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SCTransform normalization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_list <- list(Tau1, Tau2, Tau3, 
                 Cntrl1, Cntrl2, Cntrl3, Cntrl4)

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


# Save the seurat object
saveRDS(merged, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Broader class annotation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium.rds")
head(merged@meta.data)


# Cluster 0: Oligo: Sox10, Opalin, Sgk1
# Cluster 1: Astro: Aqp4
# CLuster 2: Endo: Cldn5, Sox17
# Cluster 3: ExN:Slc17a7
# Cluster 4: ExN:Slc17a7            
# Cluster 5: ExN:Slc17a7              
# Cluster 6: ExN:Slc17a7     
# Cluster 7: InN: Gad1, Gad2     
# Cluster 8: VLMC: Igf2, Aldh1a2 
# Cluster 9: Microglia: Cd53, Trem2
# Cluster 10: OPC: Pdgfra
# Cluster 11: ExN:Slc17a7
# Cluster 12: InN: Gad1, Gad2
# Cluster 13: ExN:Slc17a7
# Cluster 14: ExN:Slc17a7
# Cluster 15: ExN:Slc17a7

Stacked_VlnPlot(merged, 
                features = c("Sox10", "Opalin", 
                             "Aqp4", 
                             "Cldn5", "Sox17",
                             "Slc17a7",
                             "Gad1", "Gad2",
                             "Igf2", "Aldh1a2",
                             "Cd53", "Trem2",
                             "Pdgfra"),
                x_lab_rotate = 45,
                plot_legend = TRUE,
                colors_use = colors_celltype)


# Annotate Clusters
Idents(merged) <- merged$seurat_clusters

class <- rep(NA, length = ncol(merged))
class[which(Idents(merged) %in% c(0))] <- 'Oligo' 
class[which(Idents(merged) %in% c(1))] <- 'Astro'
class[which(Idents(merged) %in% c(2))] <- 'Endo'
class[which(Idents(merged) %in% c(3))] <- 'ExN'
class[which(Idents(merged) %in% c(4))] <- 'ExN'
class[which(Idents(merged) %in% c(5))] <- 'ExN'
class[which(Idents(merged) %in% c(6))] <- 'ExN'
class[which(Idents(merged) %in% c(7))] <- 'InN'
class[which(Idents(merged) %in% c(8))] <- 'VLMC'
class[which(Idents(merged) %in% c(9))] <- 'Microglia'
class[which(Idents(merged) %in% c(10))] <- 'OPC'
class[which(Idents(merged) %in% c(11))] <- 'ExN'
class[which(Idents(merged) %in% c(12))] <- 'InN'
class[which(Idents(merged) %in% c(13))] <- 'ExN'
class[which(Idents(merged) %in% c(14))] <- 'ExN'
class[which(Idents(merged) %in% c(15))] <- 'ExN'


class <- factor(class, 
                levels = c("ExN","InN","Astro","Oligo","Microglia","Endo","VLMC","OPC"), 
                ordered = T)

merged$class <- class

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Save the seurat object
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveRDS(merged, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium.rds")
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium.rds")

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

pdf("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Broad_class_annotation_Dimplot_Tau_RSC.pdf", height = 6, width = 8)
DimPlot_scCustom(merged, colors_use = colors_class, figure_plot = TRUE, label = F)
dev.off()


# Save representative ImageDimplots
pdf("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Broad_Class_Image_Dimplot_Tau_RSC.pdf", 
    height = 10, width = 12)
ImageDimPlot(merged, cols = colors_class, size = 3, fov = c("fov.2", "fov.6"), 
             dark.background = FALSE, split.by = "orig.ident",
             border.size = NA)
dev.off()

DefaultAssay(merged) <- "Xenium"

pdf("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Stacked_VlnPlot_markers_Broader_class_Tau_RSC.pdf",
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make Proportion plot for broader class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium.rds")

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
ggsave(filename = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Broad_Class_Proportion_Tau_RSC.pdf",
       plot,
       height = 6,
       width = 6,
       units = "in",
       dpi = 600)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform DEG across different class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the data
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium.rds")

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
  cluster.marker <- FindMarkers(merged, ident.1 = paste0(cluster,"_Tau"), 
                                ident.2 = paste0(cluster, "_Control"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox",
                                assay = "SCT") 
  DEG[[cluster]] <- cluster.marker
}

saveRDS(DEG, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/DEG_Major_class_Tau_RSC.rds")
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/DEG_Major_class_Tau_RSC.rds")

# Store in spreadsheet
write.xlsx(DEG, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/DEG_major_class_Tau_RSC.xlsx", 
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
write.xlsx(DEG_sig, file = "/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Significant_DEG_major_class_Tau_RSC.xlsx", 
           rowNames = F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Show Fos expression
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

merged <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium_new.rds")

head(merged@meta.data)
unique(merged$condition)
Assays(merged)

DefaultAssay(merged) <- "Xenium"
Idents(merged) <- merged$class
unique(merged$class)

seu_obj <- subset(merged, subset = class %in% c("ExN", "InN", "Astro", "Oligo", "Microglia"))
head(seu_obj@meta.data)
unique(seu_obj$class)

pdf("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Tau_Fos_VlnPlot_split.pdf",
    height = 5, width = 10)
VlnPlot(seu_obj, features = c("Fos"), split.by = "condition",
        assay = "Xenium",
        layer = "count",
        pt.size = 0.01,
        cols = c("Control" = "grey",
                 "Tau" = "red3"), 
        split.plot = FALSE) +
  xlab("")+
  # geom_boxplot(aes(alpha = 0.8))+
  ggpubr::labs_pubr()
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Show Fos expression
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
head(merged@meta.data)
Assays(merged) 
DefaultAssay(merged) <- "Xenium" 
names(merged@images)
Idents(merged) <- merged$class

pdf("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Fos_Exp_spatialPlot_ExN_specific.pdf", height = 8, width = 8)
ImageFeaturePlot(merged, fov = c("fov.2", "fov.6"),
                 features = c("Fos"),
                 cells = WhichCells(merged, idents = "ExN"),
                 size =2, cols = c("gray95", "firebrick2"), min.cutoff = "q10", max.cutoff = "q90",
                 scale = "all",
                 border.size = NA, dark.background = F)
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Show all the replicates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
head(merged@meta.data)
Idents(merged) <- merged$class

colors_class <- c("#E76BF3", "#00B0F6", "#00BF7D", "#A3A500", "#F8766D", "#CAB2D6", "#B2DF8A", "#A6CEE3")

ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov", "fov.2", "fov.3",
                                                            "fov.4", "fov.5", "fov.6", "fov.7"),
             dark.background = FALSE, split.by = "orig.ident",
             border.size = NA)

p1 <- ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p2 <- ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov.2"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p3 <- ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov.3"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p4 <- ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov.4"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p5 <- ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov.5"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p6 <- ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov.6"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()
p7 <- ImageDimPlot(merged, cols = colors_class, size = 2, fov = c("fov.7"), 
                   dark.background = FALSE, split.by = "orig.ident",
                   border.size = NA)+ NoLegend()

pdf("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/RSC_Tau_xeniumData_analysis/Image_Dimplot_all_replicate_Tau_major_class.pdf", height = 5, width = 5)
p1
p2
p3
p4
p5
p6
p7
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# genes detected in <10% cells per replicate
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# proportion of cell types in each biological replicate
merged <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/Tau_RSC_integrated_xenium_new.rds")
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

tau_low_detected_genes <- low_detected_genes[c("Tau1", "Tau2", "Tau3")]
control_low_detected_genes <- low_detected_genes[c("Control1", "Control2", "Control3", "Control4")]

# Find genes that are lowly detected in ALL replicates
genes_low_all_reps_tau <- Reduce(intersect, tau_low_detected_genes)
genes_low_all_reps_control <- Reduce(intersect, control_low_detected_genes)

# check how many genes
length(genes_low_all_reps_tau)
length(genes_low_all_reps_control)

# Combine the genes
genes_low_all_reps <- union(genes_low_all_reps_control, genes_low_all_reps_tau)


# check if these genes were present in the significant DEG list
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/Tau_xenium/RSC/DEG_Major_class_Tau_RSC.rds")
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
