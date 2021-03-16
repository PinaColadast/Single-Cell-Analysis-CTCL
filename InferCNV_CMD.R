install.packages("rjags") #necessary for installing infercnv 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")

install.packages("optparse")


library(infercnv)
library(biomaRt)
library(dplyr)
library(Seurat)
library(optparse)


#=====================================================================#
# arguments parser set-up
#=====================================================================#

option_list <- list(
  make_option(c("-d", "--data_dir"),
              help="data directory for input seurat object file"),
  make_option(c("-f", "--file_name"), 
              help = "file name for RDS(Seurat object) file,
              please label cell with tumor/norm in metadata column '[cell.ident]'")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)

args <- parse_args(parser, positional_arguments = 1)
working_dir <- args$data_dir
file <- args$fine_name
file_dir <- paste(working_dir, file, sep = "/")

if (file.exists(file_dir)==FALSE){
  stop("input file doesn't exist, please show correct file")
}


#set working directory
setwd(working_dir)
working_dir <- getwd()
Sys.setenv(language="en")
#----------------------------------------------------------------------------
# #Create Infercnv Object # -----------------------------------------------
#----------------------------------------------------------------------------

# 1. raw_counts matrix 
# extract from Seurat object 
# Read rds dataset
data <- readRDS(paste(working_dir, "/CTCL-tumor-normal-seurat.rds", sep =""))

raw_counts_matrix <- data@assays[["RNA"]]

# 2. cell annotation files

df_meta <- data@meta.data

cell.annotation <- data@meta.data["cell.ident"]
write.table(cell.annotation, paste(getwd(), "output/cell.annotation.txt", sep = "/"),
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
write.table(results_unique, paste(getwd(), "/output/gene_chromopos.txt", sep = ""),
            col.names = FALSE, row.names = FALSE,  sep = "\t", quote = FALSE)

# filter the counts matrix according to results of chromosome positions
counts_matrix <- raw_counts_matrix[c(results$hgnc_symbol), ]
# write.table(counts_matrix, file = paste(getwd(), "data/output/cnt_matrix", sep = "/")
#              , sep = "\t", col.names= FALSE, row.names = FALSE)


#-------------------------------------------------------------------------------
# Create InferCNV object and run ------------------------------------------
#-------------------------------------------------------------------------------
out_dir <- paste(getwd(), "/output/InferCNV/", sep = "")

if (dir.exists(out_dir)){
  out_dir <- out_dir
} else {dir.create(out_dir)}

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=paste(getwd(), "output/cell.annotation.txt", sep = "/"),
                                    delim="\t",
                                    gene_order_file= paste(getwd(), "output/gene_chromopos.txt", sep = "/"),
                                    ref_group_names=c("normal")) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

print(paste("files are written in:", out_dir, sep = " "))
      