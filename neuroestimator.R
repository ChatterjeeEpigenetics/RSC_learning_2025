################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
## function for correcting file path to adapt to file system: Argon, Topaz, IDAS, or personal
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
## this is Muhammad's wrapper functions, aes for viz, and other customized scripts
## It also loads libraries of interest and makes sure they're installed
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(patchwork)
library(Seurat);library(Signac)
registerDoMC(30)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/rsc-mem")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
### colors from Budhaditya
class.colors <- c("ExN"="#E76BF3", "InN"="#00B0F6", 
                  "Astro"="#00BF7D","Oligo"="#A3A500","Microglia"="#F8766D",
                  "Endo"="#CAB2D6","VLMC"="#B2DF8A","OPC"="#A6CEE3")
subclass.colors <- c("L2_3"="#1F77B4FF","L4"="#E377C2FF","L5"="#D62728FF","L6"="#6677d4",
                     "NP SUB"= "#720ba9","L2_3 RSP"="#FF7F0EFF","L4 RSP"="#2CA02CFF",
                     "Pvalb"="#17BECFFF","Sst"="#b7950b","Lamp5"="#98DF8AFF", "Sncg"="#FF9896FF","Vip"="#411eee",
                     "Oligo"="#da5e0e", "Astro"="#9EDAE5FF","Microglia"="#6495ED","OPC"="#FFBB78FF", 
                     "Endo"="#8C564BFF","VLMC"="#C49C94FF")

## NEUROeSTIMator genes
ne.genes <- c("Arc", "Btg2", "Coq10b", "Crem", "Dusp1", "Dusp5", "Egr1", "Egr3", 
              "Fbxo33", "Fos", "Fosl2", "Gadd45g", "Gmeb2", "Grasp", "Junb", 
              "Nr4a1", "Nr4a2", "Nr4a3", "Per1", "Rgs2", "Sertad1", "Tiparp")
################################################################################
################################################################################
################################################################################
## read data
# exp 1
sor.hc <- read_rds("data/raw/RSC_integrated_xenium_new.rds");gc()
sor.hc.clean <- GetAssayData(JoinLayers(sor.hc[["RNA"]]))
sor.hc.df <- sor.hc.clean %>% as.data.frame()
sor.hc.meta <- sor.hc@meta.data %>% 
  rownames_to_column("cell") %>%
  select(cell, orig.ident, condition, celltype, class)
# exp 2
tau.ctrl <- read_rds("data/raw/Tau_RSC_integrated_xenium_new.rds");gc()
tau.ctrl.clean <- GetAssayData(JoinLayers(tau.ctrl[["RNA"]]))
tau.ctrl.df <- tau.ctrl.clean %>% as.data.frame()
tau.ctrl.meta <- tau.ctrl@meta.data %>% 
  rownames_to_column("cell") %>%
  select(cell, orig.ident, condition, celltype, class)
################################################################################
################################################################################
################################################################################
################################################################################
library(neuroestimator)
# exp 1
sor.hc.act <- neuroestimator(sor.hc.df, species = "mmusculus") 
# exp 2
tau.ctrl.act <- neuroestimator(tau.ctrl.df, species = "mmusculus") 

## combine and save
ne.res <- inner_join(sor.hc.meta, sor.hc.act %>% rownames_to_column("cell")) %>%
  mutate(experiment = "SOR-HC") %>%
  rbind(inner_join(tau.ctrl.meta, tau.ctrl.act %>% rownames_to_column("cell")) %>%
          mutate(experiment = "Tau-Ctrl"))

write_rds(ne.res, "data/derivatives/ne-results.rds", compress = "gz")
ne.res <- read_rds("data/derivatives/ne-results.rds")
################################################################################
################################################################################
################################################################################
################################################################################
## show how many genes from the data are the main list of NE
table(ne.genes %in% rownames(sor.hc.df))
pdf("figs/neuroestimator-yes-no-genes.pdf", width = 3, height = 8)
data.frame(NE_gene = ne.genes) %>%
  mutate(present = NE_gene %in% rownames(sor.hc.df)) %>%
  ggplot(aes(x = "", y = NE_gene, fill = present)) +
  geom_tile() + scale_fill_manual(values = abstract.colors[c(4,2)]) +
  bw.theme + labs(x="", y = "NEUROeSTIMator genes")
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
exp <- unique(ne.res$experiment)
classes <- unique(ne.res$class)
subclasses <- unique(ne.res$celltype)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## for each celltype/class per experiment, build a model to predict condition
##    and get effect size
registerDoMC(4)
ne.LR.res <- foreach(ee = 1:length(exp), .combine = rbind) %dopar% {
  exp.n <- exp[ee]
  class.res <- foreach(cc = 1:length(classes), .combine = rbind) %dopar% {
    class.n <- classes[cc]
    df <- ne.res %>% filter(experiment == exp.n, class == class.n) %>%
      mutate(condition_binary = ifelse(condition %in% c("SOR","Tau"), 1,0),
             predicted_activity = scale(predicted_activity, T,T)[,1])
    coefs_table(glm(condition_binary ~ predicted_activity, data = df,family = binomial)) %>%
      mutate(experiment = exp.n, source = "class", name = class.n)
  }
  subclass.res <- foreach(sc = 1:length(subclasses), .combine = rbind) %dopar% {
    subclass.n <- subclasses[sc]
    df <- ne.res %>% filter(experiment == exp.n, celltype == subclass.n) %>%
      mutate(condition_binary = ifelse(condition %in% c("SOR","Tau"), 1,0),
             predicted_activity = scale(predicted_activity, T,T)[,1])
    coefs_table(glm(condition_binary ~ predicted_activity, data = df,family = binomial)) %>%
      mutate(experiment = exp.n, source = "subclass", name = subclass.n)
  }
  rbind(class.res, subclass.res)
}

## forest plot for class
p50 <- ne.LR.res %>%
  filter(x=="predicted_activity", source =="class",
         !name %in% c("Endo", "VLMC","OPC")) %>%
  mutate(FDR = p.adjust(pval, method="fdr"),
         sig = case_when(FDR < 0.001 ~ "***",FDR < 0.05 ~ "*")) %>%
  ggplot(aes(x = Estimate, y = factor(name,levels = names(class.colors)[8:1]), color = name,
             alpha = !is.na(sig))) +
  geom_vline(xintercept = 0, color = "red",linetype=2) +
  geom_point(show.legend = F) + 
  geom_text(aes(label = sig), vjust = 0, show.legend = F) +
  geom_errorbarh(aes(xmin = confin_min, xmax = confin_max), height = 0.2, show.legend = F) +
  scale_alpha_manual(values = c(0.3,1))+
  scale_color_manual(values = class.colors) +
  facet_wrap(~experiment, scales = "free") +
  bw.theme + labs(x="Coefficient estimate of SOR effect on predicted activity",
                  y = "", caption=paste0("*     FDR < 0.05\n","***   FDR < 0.001\n"))
## forest plot for subclass
p51 <- ne.LR.res %>%
  filter(x=="predicted_activity", source =="subclass", 
         !name %in% c("Astro", "Oligo", "Microglia", "Endo", "VLMC", "OPC")) %>%
  mutate(FDR = p.adjust(pval, method="fdr"),
         sig = case_when(FDR < 0.001 ~ "***",FDR < 0.05 ~ "*")) %>%
  ggplot(aes(x = Estimate, y = factor(name,levels = names(subclass.colors)[18:1]), color = name,
             alpha = !is.na(sig))) +
  geom_vline(xintercept = 0, color = "red",linetype=2) +
  geom_point(show.legend = F) + 
  geom_text(aes(label = sig), vjust = 0, show.legend = F) +
  geom_errorbarh(aes(xmin = confin_min, xmax = confin_max), height = 0.2, show.legend = F) +
  scale_alpha_manual(values = c(0.3,1))+
  scale_color_manual(values = subclass.colors) +
  facet_wrap(~experiment, scales = "free") +
  bw.theme + labs(x="Coefficient estimate of SOR effect on predicted activity",
                  y = "", caption=paste0("*     FDR < 0.05\n","***   FDR < 0.001\n"))

pdf("figs/ne-activity-forest-plots.pdf", width = 8, height =10)
wrap_plots(p50,p51, ncol = 1,heights = c(1,2))
dev.off()

## new fig for Jyoti, as requested
ne.LR.res %>%
  filter(x=="predicted_activity", source =="subclass", 
         !name %in% c("Astro", "Oligo", "Microglia", "Endo", "VLMC", "OPC")) %>%
  mutate(FDR = p.adjust(pval, method="fdr"),
         sigg = case_when(FDR < 0.001 ~ "***",FDR < 0.05 ~ "*")) %>%
  ggplot(aes(x = Estimate, y = factor(name,levels = names(subclass.colors)), color = name,
             alpha = !is.na(sigg))) +
  geom_vline(xintercept = 0, color = "red",linetype=2) +
  geom_point(show.legend = F) + 
  geom_text(aes(label = sigg), angle = 90, vjust = 0, show.legend = F) +
  geom_errorbarh(aes(xmin = confin_min, xmax = confin_max), height = 0.2, show.legend = F) +
  scale_alpha_manual(values = c(0.3,1))+
  scale_color_manual(values = subclass.colors) +
  facet_wrap(~experiment, scales = "free") + coord_flip() +
  labs(x="Coefficient estimate of SOR effect on predicted activity",
       y = "", caption=paste0("*     FDR < 0.05\n","***   FDR < 0.001\n")) +
  bw.theme + theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave2("figs/fig-3f.pdf", width = 12, height= 5)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## spatial viz

# exp 1
sor.hc <- AddMetaData(object = sor.hc, metadata = sor.hc.act, col.name = "neuroestimator_activity")
DefaultAssay(sor.hc) <- "Xenium"
pdf("figs/ne-spatial-viz-sor-hc.pdf", width = 12, height = 13)
ImageFeaturePlot(sor.hc, features = "neuroestimator_activity", 
                 fov = c("fov",paste0("fov.",c(2:8))),cols = c("gray95", "firebrick3"),
                 size =2, scale = "all",border.size = NA, dark.background = F) 
dev.off()

# exp 2
tau.ctrl <- AddMetaData(object = tau.ctrl, metadata = tau.ctrl.act, col.name = "neuroestimator_activity")
DefaultAssay(tau.ctrl) <- "Xenium"
pdf("figs/ne-spatial-viz-tau-hc.pdf", width = 12, height = 13)
ImageFeaturePlot(tau.ctrl, features = "neuroestimator_activity", 
                 fov = c("fov",paste0("fov.",c(2:7))),cols = c("gray95", "firebrick3"),
                 size =1.5, scale = "all",border.size = NA, dark.background = F) 
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
