rm(list=ls())
library(ggplot2)
library(data.table)
library(reshape)
library(corrplot)
library(Matrix)
library(Seurat)
library(SingleR)
library(tidyverse)
library(cowplot)
library(sctransform)
library(SeuratDisk)
library(patchwork)


reference <- LoadH5Seurat("/home/tjinyan/work_dir/CTCL/data/pbmc_multimodal.h5seurat")
data_dir <- "/home/tjinyan/work_dir/CTCL/data/Human_PBMC_10k/filtered_feature_bc_matrix"
#==========================================================================
#read in mtx 
#==========================================================================
# Read in `matrix.mtx`
counts <- readMM(paste(data_dir,
                       "matrix.mtx.gz",
                       sep= "/")
)

# Read in `genes.tsv`
genes <- read_tsv(paste(data_dir,
                        "features.tsv.gz",
                        sep= "/"), col_names = FALSE)
gene_ids <- genes$X2

# Read in `barcodes.tsv`
cell_ids <- read_tsv(paste(data_dir, "barcodes.tsv.gz", sep = "/"),
                     col_names = FALSE)$X1

rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
SeuratObj <- CreateSeuratObject(counts, project = "SeuratProject", assay = "RNA",
                           min.cells = 0, min.features = 0, names.field = 1,
                           names.delim = "_", meta.data = NULL)

DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()




#recode Azimuth names to cell ontology

# p1 = DimPlot(SeuratObj, reduction = "ref.umap", group.by = "cell.type", label = TRUE, label.size = 3, repel = TRUE)# + NoLegend()
# p2 = DimPlot(SeuratObj, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3 ,repel = TRUE) #+ NoLegend()
# p1 + p2

#saveRDS()




Seurat.STnorm.pca <- function(SeuratObj){
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-")
  SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SeuratObj.ST <- SCTransform(SeuratObj, assay = "RNA",vars.to.regress =c( "percent.mt", "S.Score", "G2M.Score" ), return.only.var.genes = FALSE)
  SeuratObj.ST <- FindVariableFeatures(SeuratObj.ST, nfeatures = 8000) 
  SeuratObj.ST <- ScaleData(SeuratObj.ST)
  SeuratObj.ST <- RunPCA(SeuratObj.ST)
  SeuratObj.ST <- RunUMAP(SeuratObj.ST, 
                          features = VariableFeatures(object = SeuratObj.ST))
  
  DefaultAssay(SeuratObj.ST) <- "SCT"
  
  return(SeuratObj.ST)
}

SeuratObj<- Seurat.STnorm.pca(SeuratObj)

anchors <- FindTransferAnchors(
  reference = reference,
  query = SeuratObj,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)
#for multi-modal data
SeuratObj <- MapQuery(
  anchorset = anchors,
  query = SeuratObj,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)
