library(Seurat)
library("tidyverse")
library("argparse")
library(future)
options(future.globals.maxSize = 500 * 1024^3)
p <- ArgumentParser(description="heatmap")
p$add_argument("-r", "--rds", type="character",required=T)
p$add_argument("-o", "--outdir", type="character",required=T)
p$add_argument("-p", "--prefix", type="character", default="scaled")
opt <- p$parse_args()


if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
obj <- readRDS(opt$rds)
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 1e6)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj), verbose = TRUE, vars.to.regress = NULL)
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".rds", sep=""))

