library("argparse")
library(Seurat)
library("ggplot2")
library(future)
options(future.globals.maxSize = 500 * 1024^3)

parser <- ArgumentParser(description='clustering and dimensionality reduction')
parser$add_argument( "-i", "--rds", type="character",required=T, default=NULL)
parser$add_argument( "-d", "--dim", type="integer", default=20)
parser$add_argument(  "--resolution", type="double", default=0.5,nargs="+")
parser$add_argument(  "--high.variable.genes", type="integer",required=F, default=2000)
parser$add_argument( "-H", "--height", type="double", default=7)
parser$add_argument("-W", "--width", type="double", default=6)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")

opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
###############################################
obj <- readRDS(opt$rds)
reduction <- "pca"
DefaultAssay(obj) <- "integrated"
obj <- FindVariableFeatures(obj, selection.method = "vst",
	nfeatures = opt$high.variable.genes, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
############# clustering #######################
obj <- FindNeighbors(obj, dims = 1:opt$dim,reduction = reduction)
obj <- FindClusters(obj, resolution = opt$resolution,reduction = reduction)
############# UMAP/tSNE ########################
obj <- RunUMAP(obj, dims = 1:opt$dim,reduction = reduction)
obj <- RunTSNE(obj, dims = 1:opt$dim,reduction = reduction)
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".rds",sep=""))
tsne <- TSNEPlot(obj, label = TRUE,label.size=8)
umap <- UMAPPlot(obj, label = TRUE,label.size=8)

png(filename=paste(opt$outdir,"/","umap_",opt$prefix,".png",sep=""), height=opt$height*300, width=opt$width*300, res=300, units="px")
print(umap)
dev.off()
png(filename=paste(opt$outdir,"/","tsne_",opt$prefix,".png",sep=""), height=opt$height*300, width=opt$width*300, res=300, units="px")
print(tsne)
dev.off()



