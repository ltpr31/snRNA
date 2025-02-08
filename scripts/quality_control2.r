library("argparse")
library(Seurat)
library("tidyverse")
library(future)
options(future.globals.maxSize = 200 * 1024^3)
parser <- ArgumentParser(description='quality control')
parser$add_argument("-i", "--rds", type="character",required=F)
parser$add_argument("--nGene.min", type="integer", default=NULL)
parser$add_argument("--nUMI.min", type="integer", default=NULL)
parser$add_argument("-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")
opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
#####################################################
obj <- readRDS(opt$rds)
nUMI.min=opt$nUMI.min
nGene.min=opt$nGene.min
obj <- subset(obj, subset = nGene > nGene.min & nUMI > nUMI.min )
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".afterQC.rds",sep=""))

