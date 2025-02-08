library("argparse")
library(Seurat)
library("tidyverse")
library(DoubletFinder)
library(future)
parser <- ArgumentParser(description='doublet removing')
parser$add_argument( "-i", "--rds", type="character",required=T, default=NULL)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")
opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
options(future.globals.maxSize = 500 * 1024^3)
obj <- readRDS(opt$rds)

############### optimal pK #############
sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:10, sct = TRUE)
sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
bcmvn_obj <- find.pK(sweep.stats_obj)
pK_bcmvn <- bcmvn_obj$pK[which.max(bcmvn_obj$BCmetric)] %>% as.character() %>% as.numeric()
############### doublet detection ######
annotations <- obj@meta.data[,"seurat_clusters"]
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.0372*nrow(obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, 
		nExp = nExp_poi.adj, reuse.pANN = F, sct = TRUE)
############### doublet removing #######
meta<-obj@meta.data
group<-colnames(meta)[length(colnames(meta))]
cs=meta[,group]=="Singlet"
cs=meta$cells[cs]
obj=subset(obj,cells =cs )
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".doubletFinder.rds",sep=""))

