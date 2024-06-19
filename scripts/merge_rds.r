library("argparse")
library(Seurat)
library("tidyverse")
library("gridExtra")
parser <- ArgumentParser(description='merge Seurat objs')
parser$add_argument( "-i", "--rds", type="character",required=T,nargs="+")
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")
opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
########### function definition ##############
MergeObject <- function(rds.list, obj.list = NULL, ...) {
	if (is.null(obj.list)) {
		rds.list <- lapply(X = rds.list, FUN = function(x) {
					readRDS(x)
				})
	} else {
		rds.list <- obj.list
	}
	merged_obj <- merge(rds.list[[1]], rds.list[-1])
	
	return(merged_obj)
}

obj=MergeObject(opt$rds)
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".rds",sep=""))










