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


# global_variables --------------------------------------------------------
project<-"CTCL"
Sys.setenv(language="en")
if (project == "CTCL"){working_dir<-"C:/Users/jtao/work_dir/CTCL/"}
print(working_dir)

getwd()
setwd(working_dir)
data_dir<-paste(working_dir,"data/output/CTCL-tumor-normal-seurat-SCT.rds", sep = "")


# data_loading ------------------------------------------------------------
data <- readRDS(gzfile(data_dir))


#============================================================================
  #ST normalize data: Seurat workflow - ST transform for RNA and CLR protein and do PCA for both
  #this is what is recommend in the vignette . Regress ccell cycle scores too
#============================================================================
Seurat.STnorm.pca <- function(SeuratObj){
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  
  SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SeuratObj.ST <- SCTransform(SeuratObj, assay = "RNA",vars.to.regress =c( "percent.mt", "S.Score", "G2M.Score" ), return.only.var.genes = FALSE)
  SeuratObj.ST <- FindVariableFeatures(SeuratObj.ST, nfeatures = 15000) 
  SeuratObj.ST <- ScaleData(SeuratObj.ST)
  SeuratObj.ST <- RunPCA(SeuratObj.ST)
  
  DefaultAssay(SeuratObj.ST) <- 'ADT'
  # we will use all ADT features for dimensional reduction
  # we set a dimensional reduction name to avoid overwriting the 
  # VariableFeatures(SeuratObj.ST) <- rownames(SeuratObj.ST[["ADT"]])
  SeuratObj.ST <- NormalizeData(SeuratObj.ST, assay = "ADT",normalization.method = 'CLR', margin = 2) %>% 
    ScaleData() %>% RunPCA(reduction.name = 'apca')
  return(SeuratObj.ST)
}

#============================================================================
# cell identity assignment: normal or tumor
#=============================================================================
tumor_function <- function (orig_ident) {
  if (grepl('mf', orig_ident)) {
    new_iden = 'tumor'
    
  } else if (grepl('nor', orig_ident)){new_iden = 'norm'
  
  }
  
  return(new_iden)
}


#==============================================
#call SingleR on clusters to assign cell types (data is SCTransformed)
#===================================================

Seurat.Singler.SCT <- function(SeuratObj, ref) {
  
  #requires reference file input for SingleR so the same function can be used for human and mice 
  #also requires seurat clusters to be present
  
  
  tt <-SeuratObj@assays[["SCT"]]@data
  singler1<-SingleR(tt, ref=ref, labels = ref$cell.type,  method = c("cluster"),
                    clusters=SeuratObj@meta.data$seurat_cluster,genes = "de", quantile = 0.8, fine.tune = TRUE,
                    tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                    check.missing = TRUE)
  
  
  ClusterCellTypes <- data.frame(singler1@listData[["pruned.labels"]])
  Clusters <- data.frame(seurat_cluster=as.numeric(rownames(ClusterCellTypes))-1, cell.type=ClusterCellTypes[1])
  
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
  
  
  return(SeuratObj)
  
}


#----------------------------------------------------------------------------
# clustering 
#----------------------------------------------------------------------------
data<-FindNeighbors(data,
                    # dim = 1:2000,
                    features = VariableFeatures(object = data),
                    k.param = 30, annoy.metric = "euclidean" 
)
# ?FindNeighbors
data<-FindClusters(data, resolution = 0.5) #default 0.8
DimPlot(data, reduction = 'umap')

# assign cell identity 

df_meta <- data@meta.data

new_iden <- c()
orig_iden <- df_meta$orig.ident
for (i in seq_along(orig_iden)) {
  new = tumor_function(orig_iden[[i]])
  new_iden <- append(new_iden, new)
}


df_meta['cell.ident'] <- new_iden
data@meta.data <- df_meta
Idents(data) <- data$seurat_clusters
data_tcell <- subset(data, ident = c("0", "3", "4", "6", "11", "13", "14"))

#---------------------------------------------------------------------------
#data normalization SCT transformation and PCA
#---------------------------------------------------------------------------
df_ref <- readRDS(gzfile("C:/Users/jtao/work_dir/CTCL/data/encode-personalized-celltypes.rds"))
data <- Seurat.Singler.SCT(data, df_ref)

Idents(data) <- data$cluster.type.singler


#find all cluster expression difference by building a large list

Idents(data_tcell) <-data_tcell$seurat_clusters
cluster.markers.list <- list()
for (i in 0:17){
  j = i+1
  cluster.markers.list[[j]] <- FindMarkers(data_tcell, ident.1 = i, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
}

head(cluster.markers.list[[18]], 30)


#===============================================================================
# further clustering on T-cell clusters
#===============================================================================

# saveRDS(data, paste(working_dir, "data/output/CTCL_clustered_SCT.rds", sep = "/"))
data <- readRDS(paste(working_dir, "data/output/CTCL-seurat-celltyped.rds", sep = "/"))
Idents(data) <- data$cell.ident
Idents(data) <- data$cluster.type.singler
Idents(data) <- data$cell.type.by.markers
DimPlot(data_tcell, reduction = "umap",
        label = TRUE, 
        # cols = MyColorsCellTypes,
        pt.size = 1, 
        repel = TRUE,
        label.size = 4)+theme_void()+theme(legend.position="top")

# head(cluster.markers.list[[3]], 30)
data_tumor <- subset(data, ident = "tumor")

# cell_types <- c("T cell lymphoma", "keratinocytes", "skin cells",
#                           "damaged/tumor NK/T cells", "NK/T cells", "keratinocytes", "monocytes", 
#                           "keratinocytes", "fibroblasts", "damaged skin cells",
#                           "stromal cells", "NK/T cells", "stromal cells",
#                           "microphages/monocytes", "NK/T cells", "damaged skin cells",
#                           "fibroblasts", "smooth muscle cells ")
# length(cell_types)

# cell_type_by_markers <- c()
# df_meta <- as.data.frame(data@meta.data)
# clusters <- as.vector(df_meta$seurat_clusters)

cluster_to_celltypes <- function(cluster_number){
  idx <- as.integer(cluster_number) +1 
  cell_type <- cell_types[[idx]]
  
  return(cell_type) 
}

for (i in seq_along(clusters)){
  cell_type_by_markers <- append(cell_type_by_markers, cluster_to_celltypes(clusters[[i]]))
  
} 

df_meta$cell.type.by.markers <- cell_type_by_markers
# df_meta <- subset(df_meta, select = -c(cell_type_by_markers))
data@meta.data <- df_meta

write.table(as.data.frame(data@meta.data),  file = paste(working_dir, 
                                                         "data/output/CTCL-meta-data.txt", sep = "/"),
sep = "\t", col.names = TRUE, row.names = TRUE)

tcell.cluster <- c("0", "3", "4", "11", "14")
data_tumor <- subset(data, subset = cell.ident == "tumor")
Idents(data) <- data$seurat_clusters
data_tcell <- subset(data_tumor, ident = tcell.cluster)

data_tcell <- FindNeighbors(data_tcell)
data_tcell <- FindClusters(data_tcell, genes.use = data_tcell@assays[["SCT"]], resolution= 0.7)
# names(data_tcell@meta.data)[7] <- "original_cluster(res0.5)

data_tcell <- RunUMAP(data_tcell, features = VariableFeatures(object = data))

head(cluster.markers.list[[4]], 20)
Idents(data_tcell) <- data_tcell$orig.ident
