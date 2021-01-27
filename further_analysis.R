library(data.table)
library(ggplot2)
library(reshape)
library(corrplot)
library(Matrix)
library(Seurat)
library(tidyverse)
library(cowplot)
library(biomaRt)
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(comprehenr)


# global_variables --------------------------------------------------------
project<-"CTCL"
Sys.setenv(language="en")
if (project == "CTCL"){working_dir<-"C:/Users/jtao/work_dir/CTCL/"}
print(working_dir)

getwd()
setwd(working_dir)
data_dir<-paste(working_dir,"data/output/CTCL-tumor-normal-seurat-SCT.rds", sep = "")

# sample_annotation -------------------------------------------------------
# ----------------------------------------------------------------------------
sample_no <- c("CTCL-2", "CTCL-5", "CTCL-6", "CTCL-8", "CTCL-12", 
               "HC-1", "HC-2", "HC-3", "HC-4")
sample_ID <- c("SC67mf2", "SC82mf5", "SC157mf62", "SC158mf8", "SC205mf12", 
               "SC50nor", "SC68nor", "SC124nor", "SC125nor")
df_id = data.frame(row.names=sample_ID, val = sample_no)
# function to assign sample no to according to sample ID
assign_sample_no <- function(sample_ID){
  sample_no <- df_id[sample_ID, ]
  return(sample_no)
}


# load the data
data <- readRDS(paste(working_dir, "data/output/CTCL-seurat-celltyped.rds", sep = "/"))
data_ctcl <- readRDS(paste(working_dir, "data/output/CTCL_tcell.rds", sep = "/"))

# select T cells only subset for analysis
Idents(data) <- data$cell.ident
data_tumor <- subset(data, ident = "tumor")

sample_IDs <- data_tcell$orig.ident
#with list comprehension 
data_tcell@meta.data[["sample_no"]] <- to_vec(for(i in sample_IDs) assign_sample_no((i)) )

tcell.cluster <- c("0", "3", "4", "11", "14")
Idents(data_tumor) <- data_tumor$seurat_clusters
data_tcell <- subset(data_tumor, ident = tcell.cluster)
data_tcell<- FindVariableFeatures(data_tcell, nfeatures = 10000)

data_tcell <- FindNeighbors(data_tcell)
data_tcell <- FindClusters(data_tcell, genes.use = data_tcell@assays[["SCT"]], resolution= 0.4)
names(data_tcell@meta.data)[7] <- "original_cluster(res0.5)"

data_tcell <- RunUMAP(data_tcell, features = VariableFeatures(object = data_tcell))

# head(cluster.markers.list[[4]], 20)
Idents(data_tcell) <- data_tcell$orig.ident

Idents(data_tcell) <-data_tcell$seurat_clusters
cluster.markers.list <- list()
for (i in 0:11){
  j = i+1
  cluster.markers.list[[j]] <- FindMarkers(data_tcell, ident.1 = i, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
}



Idents(data_tcell) <- data_tcell$orig.ident
Idents(data_tcell) <- data_tcell$seurat_clusters
DimPlot(data_tcell, reduction = "umap",
        label = TRUE, 
        # cols = MyColorsCellTypes,
        pt.size = 1, 
        repel = TRUE,
        label.size = 4)+theme_void()+theme(legend.position="top")

# data_tcell <- RunPCA(data_tcell, features = VariableFeatures(object = data))
#another differential expression test#

head(cluster.markers.list[[12]], 30)
Idents(data_tcell) <- data_tcell$orig.ident      

FeaturePlot(data_tcell, c("JUNB"),reduction = "umap",  pt.size = 1)+theme_void()


t_markers <- c("CD3D", "CD3E", "CD3F", "CD3G", "CD4", "CD8A", "CD8B")
# t_markers<- c("CD3", "CD4", "CD8a", "CD8b")
tcr_markers <-c("TRBC1", "TRBC2")
pro_markers <- c("ACTG1", "ANP32B", "ATP5C1", "DUT", "HMGN1", "HN1",
                 "NPM1", "NUSAP1", "PCNA", "PPA1", "PPIA", "PSMB2", "RAN",
                 "RANBP1", "SET", "SMC4", "STMN1") #proliferating cell markers
tumor_markers <- c
pro_list <- list(pro_markers)

data_tcell <- AddModuleScore(object = data_tcell, features = pro_list, ctrl = 5,
                             name = "proliferation_features")
umap <- as.data.frame(data_tcell@reductions[["umap"]]@cell.embeddings)
pro_scatter <- ggplot (umap, 
                       aes(x= UMAP_1, y=UMAP_2, 
                       )) 
pro_scatter + labs(colour = "proliferation markers avg. expression") + theme_void()+ theme(legend.position = "bottom") +geom_point(size = 0.5, 
                                                                                                                                   aes(color = data_tcell$proliferation_features1)) + scale_color_gradientn(colours = rainbow(5))
Idents(data_tcell)<- data_tcell$seurat_clusters

DoHeatmap(data_tcell, 
          features = tcr_markers, 
          size = 4)  + labs(colours = "Expression") + theme(legend.position= "bottom") 
#keep subsetting
ctcl_cluster <- c("0", "1", "2", "3", "4", "6", "7",
                  "8", "11")
Idents(data_tcell) <- data_tcell$seurat_clusters
data_ctcl <- subset(data_tcell, ident = ctcl_cluster)

# transform to SCE object

saveRDS(data_tcell, file = paste(working_dir, "data/output/CTCL_tcell.rds", sep = "/"))
sampletype <- unique(data_ctcl$sample_no)
sampletype

Idents(data_tcell) <- data_tcell$sample_no
sample_1.sce <- as.SingleCellExperiment(subset(data_tcell, idents = sampletype[1]))
sample_2.sce <- as.SingleCellExperiment(subset(data_tcell, idents = sampletype[2]))
sample_3.sce <- as.SingleCellExperiment(subset(data_tcell, idents = sampletype[3]))
sample_4.sce <- as.SingleCellExperiment(subset(data_tcell, idents = sampletype[4]))
sample_5.sce <- as.SingleCellExperiment(subset(data_tcell, idents = sampletype[5]))


sampletype
p1 <- visualiseExprs(sample_3.sce,plot = "pairwise", feature_subset = c("TRBC1", "TRBC2"),
                     group_by = "seurat_clusters") 
p1

p1.2 <- visualiseExprs(sample_4.sce,plot = "pairwise", feature_subset = c("TRBC1", "TRBC2"),
                       group_by = "seurat_clusters")
plot_row <- plot_grid(p1,p1.2, align = "h", rel_widths = c(1,1)) 

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


# to quantify
Idents(data_ctcl) <- data_ctcl$sample_no
sample.seurat <- list()
for (i in 1:length(sampletype)){
  print(i)
  sample.seurat[[i]] <- subset(data_ctcl, ident = sampletype[[i]])
}

aa <- sample.seurat[[1]]

VlnPlot(data_ctcl, c("TRBC1", "TRBC2"))

df<- aa@assays[["SCT"]]@data[tcr_markers,]
df_1 <- df %>% filter("TRBC1" > 1)
df[, df["TRBC1",] >1 & df["TRBC2", ] >1]
