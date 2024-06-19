library("argparse")
library(Seurat)
library(monocle)
library(dplyr)
parser <- ArgumentParser(description='single cell trajectory construction')
parser$add_argument("-i", "--rds", type="character",required=T)
parser$add_argument("--gene.select.method", type="character",required=F, default="monocle")
parser$add_argument("--gene.list", type="character",required=F, default=NULL)
parser$add_argument("-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="tissue")
opt <- parser$parse_args()
#####################################################
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create directory:",opt$outdir,sep=""))
	}
}
###############################################
options(future.globals.maxSize = 500 * 1024^3)
############## newCellDataSet ##############
scRNA <- readRDS(opt$rds)
DefaultAssay(scRNA) <- "RNA"
assay="RNA"
expr_matrix <- scRNA@assays[[assay]]@counts
p_data <- scRNA@meta.data
f_data <- data.frame(gene_short_name = row.names(scRNA),row.names = row.names(scRNA))
#####CDS construction
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
		phenoData = pd,
		featureData = fd,
		lowerDetectionLimit = 0.5,
		expressionFamily = negbinomial.size())
cds<-estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

############ setOrderingFilter ##############

if(opt$gene.select.method=="Customization"){
	gene.list<-read.table(opt$gene.list,header = T,sep="\t",row.names=1)
	disp.genes=rownames(gene.list)	
	cds <- setOrderingFilter(cds, disp.genes)
}else if(opt$gene.select.method=="monocle"){
	disp_table <- dispersionTable(cds)
	disp.genes <- disp_table$gene_id
	cds <- setOrderingFilter(cds, disp.genes)
	write.table(disp_table, file=paste(opt$outdir,"/",opt$prefix,".dispersionTable.tsv",sep=""), row.names = F, col.names = T, quote = F)
}
############ reduceDimension ###############
cds <- reduceDimension(cds, max_components = 2,
		method = 'DDRTree')
############ orderCells ####################
cds <- orderCells(cds)
saveRDS(cds, file = paste(opt$outdir,"/",opt$prefix,".CelldataSet.rds",sep=""))
Time_diff <- differentialGeneTest(cds[disp.genes,], cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff=arrange(Time_diff,qval)
write.table(data.frame(ID=rownames(Time_diff),Time_diff),file=paste(opt$outdir,"/",opt$prefix,".Time_diff_all.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)

