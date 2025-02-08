library("argparse")

parser <- ArgumentParser(description='Integration')
parser$add_argument( "-i", "--rds", type="character",required=T, default=NULL)
parser$add_argument( "--integration.reduction", type="character",required=F, default="cca")
parser$add_argument( "-b", "--batch.id", type="character",required=F)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument( "-p", "--prefix", type="character", default="demo")
opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
###############################################
library(Seurat)
library("glmGamPoi")
library(future)
options(future.globals.maxSize = 1000 * 1024^3)

###############################################
obj <- readRDS(opt$rds)
DefaultAssay(obj)<-"RNA"
obj <- SCTransform(obj, 
			vars.to.regress = NULL,
			method = "glmGamPoi",
			variable.features.n = 2000,
			return.only.var.genes=FALSE,
			verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(object = obj),npcs = 100)

##############################integration#########################################################
reduction = "pca"
obj.list <- SplitObject(obj, split.by =opt$batch.id)
	obj.list <- lapply(X = obj.list, FUN = function(x) {
					x <- SCTransform(x, 
							vars.to.regress = NULL,
							method = "glmGamPoi",
							variable.features.n = 2000,
							return.only.var.genes=FALSE,
							verbose = FALSE)
				})		
		features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
		obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
		obj.list <- lapply(X = obj.list, FUN = function(x) {
					x <- RunPCA(x, features = rownames(x), verbose = FALSE)
				})
		anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
				anchor.features = features, dims = 1:30, reduction =  opt$integration.reduction, k.anchor = 5)
		obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)

	reduction <- "pca"
	DefaultAssay(obj) <- "integrated"
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".integrated.rds",sep=""))




