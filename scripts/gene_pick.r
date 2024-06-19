library("argparse")
p <- ArgumentParser(description="pick genes for trajectory construction")
p$add_argument("-i", "--input", type="character",required=T)
p$add_argument("-m", "--inputck", type="character",required=T)
p$add_argument("-o", "--output", type="character",required=T)
p$add_argument("--ntr", type="double", default=2000)
p$add_argument("--nck", type="double", default=500)
opt <- p$parse_args()
tr <- read.table(opt$input, sep = "\t", header = T)
tr <- tr[!grepl("^gene", tr$ID), ]
tr <- tr[tr$qval < 0.01 , ]
tr <- tr[order(tr$qval),]
tr <- tr[1:opt$ntr,]

ck <- read.table(opt$inputck, sep = "\t", header = T)
ck <- ck[!grepl("^gene", ck$ID), ]
ck <- ck[ck$qval < 0.01 , ]
ck <- ck[order(ck$qval),]
ck <- ck[1:opt$nck,]

id <- ck$ID

######### gene pick
result <- tr[!tr$ID %in% id, ]
res <- result$ID
res <- as.matrix(res)
colnames(res)[1] <- "ID"
res <- na.omit(res) 
write.table(res, opt$output, col.names = T, row.names = F, quote = F)












