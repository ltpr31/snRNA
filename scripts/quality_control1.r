library("argparse")
library(Seurat)
library("tidyverse")
library(future)
options(future.globals.maxSize = 200 * 1024^3)
parser <- ArgumentParser(description='quality control')
parser$add_argument( "-d", "--data.dir", type="character",required=F)
parser$add_argument(  "--metadata.col.name ", type="character",required=F, default=NULL,nargs="+")
parser$add_argument(  "--metadata.value ", type="character",required=F, default=NULL,nargs="+")
parser$add_argument( "--nGene.min", type="integer", default=NULL)
parser$add_argument(  "--nUMI.min", type="integer", default=NULL)
parser$add_argument(  "--log10GenesPerUMI", type="double",required=F, default=NULL)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")
opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
################### load count ############################
assay = "RNA"
	obj.data <- Read10X(data.dir = opt$data.dir,gene.column=1)
	obj <- CreateSeuratObject(counts = obj.data, assay = assay,project="snRNA",min.cells=3)
	rm(list=c("obj.data"))

############ add meta data to seurat obj ##################
m=cbind(opt$metadata.value,opt$metadata.col.name)
	for(i in 1:nrow(m)){
		x=m[i,]
		obj <- AddMetaData(
				object = obj,
				metadata = x[1],
				col.name = x[2]
		)
	}
#######################################
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
metadata <- obj@meta.data
metadata$cells <- rownames(metadata)

metadata <- metadata %>%
		dplyr::rename(
				nUMI = nCount_RNA,
				nGene = nFeature_RNA)
obj@meta.data <- metadata

nUMI.min=opt$nUMI.min
nGene.min=opt$nGene.min

obj <- subset(obj, subset = nGene > nGene.min & nUMI > nUMI.min )
obj <- subset(obj, subset = log10GenesPerUMI > opt$log10GenesPerUMI )

saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".afterQC.rds",sep=""))

