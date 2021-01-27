# rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape)
library(corrplot)
library(Matrix)
library(Seurat)
library(SingleR)
library(tidyverse)
library(cowplot)
library(SingleCellExperiment)
library(CiteFuse)
library(hrbrthemes)


Sys.setenv(language="en")
getwd()
setwd("C:/Users/jtao/work_dir/CTCL/")
data.all<- readRDS(paste(getwd(), "data/CTCL-tumor-normal-seurat.rds", sep = "/"))
data <- readRDS(paste(getwd(), "data/output/CTCL.cell-typed.rds", sep = "/"))

all.genes <- rownames(data.all)


data.all[["percent.mt"]] <- PercentageFeatureSet(data.all, pattern = '^MT-')
data.all <- subset(data.all, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 10000)
data <- ScaleData(data)
data.all <- LogNormalize(data.all)
data.all <- ScaleData(data.all)
df.meta.all <- data.all@meta.data
cluster <- data$RNA_snn_res.0.5
df.meta.all[["cluster"]] <- cluster
data.all@meta.data <- df.meta.all


MyColorsCellTypes <- c("B cell" = "#1F77B4", "CD8 T cell"="#AEC7E8",  "fibroblast"="#FF7F0E", "myeloid dendritic cell"="#FFBB78", 
                       "plasma cell"="#2CA02C", "epithelial cell"="#98DF8A", "natural killer cell"="#D62728", "macrophage"="#FF9896", "CD4 T cell"="#9467BD", "endothelial cell"="#C5B0D5", "regulatory T cell"="#8C564B", "monocyte"="#E377C2", "unclassified"="#BCBD22")

fix.sc = scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 4.), breaks = c(0, 1.0, 2.0), guide = guide_colourbar(direction = "horizontal"))
DimPlot(data, reduction = "umap",
        label = TRUE, 
        # cols = MyColorsCellTypes,
        pt.size = 1, 
        repel = TRUE,
        label.size = 4)+theme_void()+theme(legend.position="top")



FeaturePlot(data, c("MS4A1"),reduction = "umap",  pt.size = 1)+theme_void()
Idents(data) <- data@meta.data[["RNA_snn_res.0.5"]]
Idents(data) <- data@meta.data["cluster.type.singler"]

data.sce <-as.SingleCellExperiment(data)
visualiseExprs(data.sce, 
               plot = "pairwise",
               feature_subset = c("TRAC", "TRBC1"), group_by ="cluster.type.singler" )
data.sce[["RNA_snn_res.0.5"]]

sce.cluster8 <- data.sce[, data.sce$RNA_snn_res.0.5==8]
sce.cluster12<-data.sce[, data.sce$RNA_snn_res.0.5==12]
sce.patient1 <- data.sce[, data.sce$orig.ident=="SC205mf12"]

visualiseExprs(sce.patient1,plot = "pairwise", feature_subset = c("TRAC", "TRBC2"), group_by ="cluster.type.singler")

Idents(data) <- data$RNA_snn_res.0.5
data_tumor_ctcl <- subset(data, ident = c("8","12"))
data.sce
data_tumor_ctcl



Idents(data) <- data$RNA_snn_res.0.5
data_subset <- subset(data, ident=c("8", "12"))
Idents(data) <- data$cell.ident
data_tumor <- subset(data, ident = "tumor")
Idents(data_subset) <- data_subset$orig.ident

sample_1 <- subset(data_subset, idents = "SC205mf12")

# data.sce.subset <-as.SingleCellExperiment(data_subset)
sample.sce <- list()
sampletype <- unique(data_tumor$orig.ident)
sampletype

for (i in 1:length(sampletype)){
  print(i)
  sample.sce[[i]] <- as.SingleCellExperiment(subset(data_subset,
                                                  ident = sampletype[[i]]))
}
Idents(data_tumor) <- data_tumor$orig.ident
sample_1.sce <- as.SingleCellExperiment(subset(data_tumor, idents = sampletype[1]))
sample_2.sce <- as.SingleCellExperiment(subset(data_tumor, idents = sampletype[2]))
sample_3.sce <- as.SingleCellExperiment(subset(data_tumor, idents = sampletype[3]))
sample_4.sce <- as.SingleCellExperiment(subset(data_tumor, idents = sampletype[4]))
sample_5.sce <- as.SingleCellExperiment(subset(data_tumor, idents = sampletype[5]))


sampletype
p1 <- visualiseExprs(sample_5.sce,plot = "pairwise", feature_subset = c("TRBC1", "TRBC2"),
               group_by = "seurat_clusters") 
p1

p1.2 <- visualiseExprs(sample_4.sce,plot = "pairwise", feature_subset = c("TRBC1", "TRBC2"),
               group_by = "seurat_clusters")
plot_row <- plot_grid(p1,p1.2, align = "h", rel_widths = c(1,1)) 
sampletype[[5]]
title <- ggdraw() + 
  draw_label(
    paste(sampletype[3], sampletype[4], sep = "                                                          "),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
plot_grid(title, plot_row, ncol=1, rel_heights = c(0.05, 1))

Idents(data) <- data$RNA_snn_res.0.5
p3 = RidgePlot(data, "TRGC2")+theme_bw()+theme(legend.text = element_text(6))
p3
p4 = RidgePlot(SeuratObj1, "CD40")+theme_bw()#+theme(legend.text=element_text(size=20))+scale_fill_manual(values=MyColorsCellTypes)

#heatmap#
#on T cell markers 
t_markers <- c("CD3D", "CD3E", "CD3F", "CD3G", "CD4", "CD8A", "CD8B")
t_markers<- c("CD3", "CD4", "CD8a", "CD8b")
tcr_markers <-c("TRBC1", "TRBC2")
pro_markers <- c("ACTG1", "ANP32B", "ATP5C1", "DUT", "HMGN1", "HN1",
                "NPM1", "NUSAP1", "PCNA", "PPA1", "PPIA", "PSMB2", "RAN",
                "RANBP1", "SET", "SMC4", "STMN1") #proliferating cell markers
pro_list <- list(pro_markers)

# add modular score
data_tcell <- AddModuleScore(object = data_tcell, features = pro_list, ctrl = 5,
                        name = "proliferation_features")
umap <- as.data.frame(data_tcell@reductions[["umap"]]@cell.embeddings)
pro_scatter <- ggplot (umap, 
                       aes(x= UMAP_1, y=UMAP_2, 
                                  )) 
pro_scatter + labs(colour = "proliferation markers avg. expression") + theme_void()+ theme(legend.position = "bottom") +geom_point(size = 0.5, 
                                       aes(color = data_tcell$proliferation_features1)) + scale_color_gradientn(colours = rainbow(5))
Idents(data_tcell)<- data_tcell$seurat_clusters

Idents(data_tumor) <- data$seurat_clusters
DoHeatmap(data_tcell, 
          features = t_markers, 
          size = 4)  + labs(colours = "Expression") + theme(legend.position= "bottom")  
?DoHeatmap
data@assays
saveRDS(data, file = paste(working_dir, "data/output/CTCL.cell-typed.rds" ,sep = ""))
save.image(file = paste(working_dir, "15012021.RData"))

working_dir<- ("C:/Users/jtao/work_dir/CTCL/")
setwd(working_dir)

data_all <- readRDS(gzfile(paste(working_dir, "data/output/CTCL.cell-typed.rds", sep = "")))
data_all <- NormalizeData(data_all)
all.genes <- rownames(data_all)
# length(all.genes)
data_all <-FindVariableFeatures(data_all, selection.method = "vst", nfeatures = 10000)
data_all<- ScaleData(data_all)


