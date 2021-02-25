library(Seurat)
library(devtools)
install_github("navinlabcode/copykat")
library(copykat)

# test copykat if it works 
#copykat.test <- copykat(rawmat=exp.rawdata, sam.name="test")

Sys.setenv(language="en")
getwd()
setwd("/home/tjinyan/work_dir/CTCL/")
data.all<- readRDS(paste(getwd(), "data/CTCL-tumor-normal-seurat.rds", sep = "/"))
exp.rawdata <- as.matrix(data.all@assays$RNA@counts)
#write.table(as.matrix(data.all@assays$RNA@counts), paste(getwd(), "data/exp.rawdata.txt",
 #                                                        sep = ''), sep = "\t", quote = FALSE,
  #          row.names = TRUE)
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=8)
