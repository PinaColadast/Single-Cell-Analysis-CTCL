#install dependencies 
install.packages("rjags") #necessary for installing infercnv 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")
 
library(infercnv)
library(biomaRt)
library(dplyr)
library(Seurat)

#set working directory
setwd("C:/Users/jtao/work_dir/CTCL")
working_dir <- getwd()
Sys.setenv(language="en")
#----------------------------------------------------------------------------
# #Create Infercnv Object # -----------------------------------------------
#----------------------------------------------------------------------------

# 1. raw_counts matrix 
# extract from Seurat object 
# Read rds dataset
data <- readRDS(paste(working_dir, "/data/CTCL-tumor-normal-seurat.rds",sep =""))

data = data
raw_counts_matrix <- data@assays[["RNA"]]

# 2. cell annotation files

tumor_function <- function (orig_ident) {
  if (grepl('mf', orig_ident)) {
    new_iden = 'tumor'
    
  } else if (grepl('nor', orig_ident)){new_iden = 'norm'
  
  }
  
  return(new_iden)
}

new_iden <- c()
orig_iden <- data$orig.ident
for (i in seq_along(orig_iden)) {
  new = tumor_function(orig_iden[[i]])
  new_iden <- append(new_iden, new)
}

df_meta <- data@meta.data
df_meta['cell.ident'] <- new_iden
df_meta["cell.ident"]
cell.annotation <- data@meta.data["orig.ident"]
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
results[results$hgnc_symbol %in% rep_gene[rep_gene$Freq>1, ]$Var1, ]

#clear replicates
results_unique <- results[!duplicated(results$hgnc_symbol), ]

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
out_dir <- paste(getwd(), "data/output/InferCNV/", sep = "/")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=paste(getwd(), "data/output/cell.annotation.txt", sep = "/"),
                                    delim="\t",
                                    gene_order_file= paste(getwd(), "data/output/gene_chromopos.txt", sep = "/"),
                                    ref_group_names=c("SC124nor","SC50nor", "SC68nor")) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)
