library("argparse")
library(Seurat)
library("glmGamPoi")
library(future)
options(future.globals.maxSize = 100 * 1024^3)
parser <- ArgumentParser(description='single sample clustering')
parser$add_argument( "-i", "--rds", type="character",required=T, default=NULL)
parser$add_argument(  "--resolution", type="double", default=0.5)
parser$add_argument(  "--high.variable.genes", type="integer",required=F, default=2000)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")

opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
#######################################################
obj <- readRDS(opt$rds)
DefaultAssay(obj)<-"RNA"

obj <- SCTransform(obj, 
			vars.to.regress = NULL,
			method = "glmGamPoi",
			variable.features.n = opt$high.variable.genes,
			return.only.var.genes=FALSE,
			verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(object = obj),npcs = 100)
reduction = "pca"
#################clustering###########
obj <- FindNeighbors(obj, dims = 1:100,reduction = reduction)
obj <- FindClusters(obj, resolution = opt$resolution,reduction = reduction, n.start = 10)
obj <- RunUMAP(obj, dims = 1:10,reduction = reduction,min.dist = 0.05, n.neighbors = 5, seed.use = 100)
obj <- RunTSNE(obj, dims = 1:20, check_duplicates = FALSE,reduction = reduction)

saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".rds",sep=""))
