#install dependencies 
#install.packages("rjags", repos ="https://CRAN.R-project.org/package=rjags") #necessary for installing infercnv 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repo="http://cran.us.r-project.org")
#BiocManager::install("infercnv")
#install.packages("dplyr")
#install.packages("biomaRt", "http://cran.us.r-project.org")

library(infercnv)
library(biomaRt)
library(dplyr)
library(Seurat)

#set working directory
setwd("/home/tjinyan/work_dir/CTCL")
working_dir <- getwd()
Sys.setenv(language="en")
#----------------------------------------------------------------------------
# #Create Infercnv Object # -----------------------------------------------
#----------------------------------------------------------------------------

# 1. raw_counts matrix 
# extract from Seurat object 
# Read rds dataset
data <- readRDS(paste(working_dir, "/data/cite-seq-ctcl-ctrl.rds",sep =""))
# run SCT to take out doplets#
Seurat.STnorm.pca <- function(SeuratObj){
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  
  SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SeuratObj.ST <- SCTransform(SeuratObj, assay = "RNA",vars.to.regress =c( "percent.mt", "S.Score", "G2M.Score" ), return.only.var.genes = FALSE)
  SeuratObj.ST <- FindVariableFeatures(SeuratObj.ST, nfeatures = 15000) 
  SeuratObj.ST <- ScaleData(SeuratObj.ST)
  SeuratObj.ST <- RunPCA(SeuratObj.ST)
  return(SeuratObj.ST)
}

#data <- Seurat.STnorm.pca(data)
#data <- readRDS(paste(getwd(), "/data/cite-seq_ctcl_ctrl.SCT.rds", sep = ""))
#saveRDS(data, file = paste(getwd(), "/data/cite-seq_ctcl_ctrl.SCT.rds", sep = ""))
raw_counts_matrix <- data@assays[["RNA"]]

# 2. cell annotation files

tumor_function <- function (orig_ident) {
  if (grepl('tumor', orig_ident)) {
    new_iden = 'tumor'
    
  } else if (grepl('norm', orig_ident)){new_iden = 'norm'
  
  }
  
  return(new_iden)
}

new_iden <- c()
orig_iden <- as.vector(rownames(data@meta.data))
for (i in seq_along(orig_iden)) {
  new = tumor_function(orig_iden[[i]])
  new_iden <- append(new_iden, new)
}

df_meta <- data@meta.data
df_meta$cell.ident <- new_iden
#df_meta["cell.ident"]
cell.annotation <- df_meta["cell.ident"]
write.table(cell.annotation, paste(getwd(), "data/output/cell.annotation.txt", sep = "/"),
                col.names = FALSE, sep = "\t" )

#------------------------------------------------------
# 3. gene order file (order on Chromosomes)
#------------------------------------------------------
all.genes <- c(rownames(data))

df_gene <- as.data.frame(table(all.genes), row.names = NULL)
colnames(df_gene) <- c("hgnc_symbol", "freq")
head(df_gene,10)
# retrieve chromosomes positions given a list of genes 
query_genes = df_gene
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results_id <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'chromosome_name',
                     'start_position', 'end_position'),
      filters = "hgnc_symbol", 
      values = df_gene$hgnc_symbol, 
      mart = ensembl)

chromo_list <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15", "16", "17", "18", "19",
                 "20", "21", "22", "X", "Y")

results <- results_id %>% filter(chromosome_name %in% chromo_list)
results <- results %>% select(hgnc_symbol, chromosome_name, start_position, end_position)

#check if any duplicates in gene position
rep_gene <- data.frame(table(results$hgnc_symbol))
#results[results$hgnc_symbol %in% rep_gene[rep_gene$Freq>1, ]$Var1, ]

#clear replicates
results_unique <- results[!duplicated(results$hgnc_symbol), ]
results_unique
# write table of gene notations
write.table(results_unique, paste(getwd(), "data/output/gene_chromopos.txt", sep = "/"),
            col.names = FALSE, row.names = FALSE,  sep = "\t", quote = FALSE)

# filter the counts matrix according to results of chromosome positions
counts_matrix <- raw_counts_matrix[c(results$hgnc_symbol), ]
# write.table(counts_matrix, file = paste(getwd(), "data/output/cnt_matrix", sep = "/")
#              , sep = "\t", col.names= FALSE, row.names = FALSE)


#-------------------------------------------------------------------------------
# Create InferCNV object and run ------------------------------------------
#-------------------------------------------------------------------------------
out_dir <- paste(getwd(), "/data/output/InferCNV/", sep ="")
if (dir.exists(out_dir)){
	    out_dir <- out_dir
  } else {dir.create(out_dir)}

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=paste(getwd(), "data/output/cell.annotation.txt", sep = "/"),
                                    delim="\t",
                                    gene_order_file= paste(getwd(), "data/output/gene_chromopos.txt", sep = "/"),
                                    ref_group_names=c("norm"))
infercnv_obj = infercnv::run(infercnv_obj, cutoff=0.1,
			     # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)
