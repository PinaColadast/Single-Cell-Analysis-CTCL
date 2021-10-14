
rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape)
library(corrplot)
library(Matrix)
library(singlecell)
library(Seurat)
library(SingleR)
library(tidyverse)
library(cowplot)
library(biomaRt)
library(SingleCellExperiment)

setwd("C:/data/10x datasets/To-process/QianLambrechts-Pan-cancer/")


#===
#Seurat includes a list of cellcycle markers 
#==================  
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#==============================
#plot from Seurat data sets
#================================

#DataMatrix <- Read10X(data.dir="C:/data/10x datasets/To-process/Lymphoma-Bcells/DLBCL3")
#MetaData <- read.table("C:/data/10x datasets/To-process/QianLambrechts-Pan-cancer/BC_metadata.csv", sep=",", header=TRUE)
#rownames(MetaData)<- MetaData$CellID

#SeuratObj <- CreateSeuratObject(counts=DataMatrix, project="DLBC", assay="RNA")
SeuratObj <- readRDS("C:/data/10x datasets/To-process/GSE128531_CTCL/CTCL-tumor-normal-seurat.rds")
SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-")
#SeuratObj <- AddMetaData(SeuratObj, MetaData)
SeuratObj <- subset(SeuratObj, subset = nFeature_RNA > 300 & percent.mt < 10) 

SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)

#add cell cycle scores and regress them out 
SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


SeuratObj <- ScaleData(object = SeuratObj, vars.to.regress = c("S.Score", "G2M.Score"), use.umi = TRUE)
SeuratObj <- RunPCA(object = SeuratObj, pc.genes = SeuratObj@var.genes, do.print = TRUE, pcs.print = 1:5, 
                    genes.print = 5)
SeuratObj<- RunTSNE(object = SeuratObj, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)
SeuratObj <- FindNeighbors(SeuratObj, dims = 1:10) 
SeuratObj <- FindClusters(object = SeuratObj, resolution=0.5)
#SeuratObj1 <- AddMetaData(SeuratObj1, MDsubset)
SeuratObj <- RunUMAP(SeuratObj, dim = 1:10)
DimPlot(SeuratObj, reduction = "umap")


fix.sc = scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 4.), breaks = c(0, 1.0, 2.0), guide = guide_colourbar(direction = "horizontal"))

p2 = FeaturePlot(SeuratObj, "FAP", pt.size = 1, reduction = "umap")+theme_void()+theme(legend.text=element_text(size=20))+fix.sc
p2

VlnPlot(SeuratObj,"TRBC2")

#-----------------------
#run SingleR with personalized Encode - requires having an sce object
#--------------------------------------------------


SceObj <- readRDS( "C:/data/10x datasets/Sce objects/Lambrechts_BC-cell-types-CL-sce.rds")
SeuratObj <- readRDS("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/ToTransfer/Citeseq-lung-nonimmune-panel-with-doublet-calls-seurat.rds")
#SeuratObj <- subset(SeuratObj, SAMPLE_TISSUE_OF_ORIGIN != "mouse")
#saveRDS(SeuratObj, "C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/ToTransfer/Citeseq-lung-nonimmune-panel-with-doublet-calls-seurat.rds")

df.encode <- BlueprintEncodeData()
annotation <- read.table("C:/data/Broad_CCLE_TPM/Annotation-Blueprint-encode.txt", header=TRUE, sep = "\t")
annotation$cell.type <- recode(annotation$cell.type, "natural killer cells"="natural killer cell")
df.encode@colData@listData<-annotation

#df.encode <- df.encode[, df.encode$cell.type %in% c("muscle cell", "hematopoietic stem cell", "B cell", "CD4 T cell","CD8 T cell", "endothelial cell", "eosinophil", "epithelial cell", "monocyte", "myeloid dendritic cell","natural killer cell", "neutrophil","plasma cell", "regulatory T cell" , "fibroblast", "macrophage")]
#saveRDS(df.encode,"C:/data/SingleR-reference-sets/encode-personalized-celltypes.rds")

df.encode <- readRDS("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/ToTransfer/encode-personalized-celltypes.rds")
tt <-as.matrix(SeuratObj@assays[["SCT"]]@data)
singler1<-SingleR(tt, ref=df.encode, labels = df.encode$cell.type,  method = c("cluster"),
                  clusters=SeuratObj@meta.data$seurat_cluster,genes = "de", quantile = 0.8, fine.tune = TRUE,
                  tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                  check.missing = TRUE)


ClusterCellTypes <- data.frame(singler1@listData[["pruned.labels"]])
Clusters <- data.frame(seurat_cluster=rownames(ClusterCellTypes), cell.type=ClusterCellTypes[1])
colnames(Clusters)[2]<-"cell.type"
rm(tt)



MetaDataM <- data.frame(CellID=rownames(SeuratObj@meta.data),SeuratObj@meta.data)
#colnames(MetaDataM)[29]<-"cluster_subtype"
MetaDataM$cluster.type.singler <- "unclassified"
for (i in 1: nrow(MetaDataM)){
  for (k in 1: nrow(Clusters)){
    if (MetaDataM$seurat_cluster[i]==Clusters$seurat_cluster[k]) MetaDataM$cluster.type.singler[i]<-Clusters$cell.type[k]
    
    
  } 
}


SeuratObj <- AddMetaData(SeuratObj, MetaDataM)
SeuratObj <- SetIdent(SeuratObj, value="cluster.type.singler")
MyColorsCellTypes <- c("B cell" = "#1F77B4", "CD8 T cell"="#AEC7E8",  "fibroblast"="#FF7F0E", "myeloid dendritic cell"="#FFBB78", 
                       "plasma cell"="#2CA02C", "epithelial cell"="#98DF8A", "natural killer cell"="#D62728", "macrophage"="#FF9896", "CD4 T cell"="#9467BD", "endothelial cell"="#C5B0D5", "regulatory T cell"="#8C564B", "monocyte"="#E377C2", "unclassified"="#BCBD22")
p1 <-DimPlot(SeuratObj, reduction = "wnn.umap", cols = MyColorsCellTypes)
SeuratObj <- SetIdent(SeuratObj, value="orig.ident")
p2<-DimPlot(SeuratObj, reduction = "wnn.umap")
plot_grid(p1,p2)

FeatureScatter(SeuratObj, "FAP", "FN1", cells = WhichCells(SeuratObj, idents  = "tumor"))+theme_bw()+theme(legend.text=element_text(size=20))

#saveRDS(SeuratObj, "C:/data/10x datasets/Seurat objects/Lambrechts_BC-cell-types-singler.rds")


SCEObj <- as.SingleCellExperiment(SeuratObj)
#saveRDS(SCEObj, "C:/data/10x datasets/Sce objects/Lambrechts_BC-cell-types-singler.rds")


#====================
#biomarts for gene names conversions

#genesMat <- read.table("C:/Users/i0301956/OneDrive - Sanofi/QIAGEN_ProcessedData/GSE117156_GPL18573/features.tsv", sep="\t")
#colnames(genesMat)[1]<-"ensembl_gene_id"

#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#genes <- genesMat$ensembl_gene_id
#df$gene_name <- NA
#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
#                                                          "external_gene_name", "description"),values=genes,mart= mart)
#genesMat2<-merge(genesMat, G_list, by = "ensembl_gene_id")
#write.table(dgBioturing2, "de_genes_2182019_neutrophils_blood_vs_tumor_bioturing_venice_gene_id.txt", sep="\t", row.names = FALSE)



