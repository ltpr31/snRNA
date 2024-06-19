library("argparse")
library(Seurat)
library(future)
options(future.globals.maxSize = 100 * 1024^3)
parser <- ArgumentParser(description='subset data')
parser$add_argument( "-i", "--rds", type="character",required=T, default=NULL)
parser$add_argument( "-s","--subset", type="character",required=F, default=NULL)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument( "-p", "--prefix", type="character", default="subset")
opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
obj <- readRDS(opt$rds)
levels(Idents(obj))
t=paste("obj=subset(x = obj, subset=",opt$subset, ",cells=NULL,idents=NULL,downsample=NULL)")
eval(parse(text = t))
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".rds",sep=""))




