DEGList <-
function (counts = matrix(0, 0, 0), lib.size = colSums(counts,na.rm=TRUE), norm.factors = rep(1, ncol(counts)), group = rep(1, ncol(counts)), length=NULL, length.gene="length.gene.hg19",length.isoform="length.isoform.hg19", rpkms=NULL, genes = NULL, remove.zeros = FALSE, level="gene"){

counts <- as.matrix(counts)

if(is.null(genes)){
genes<-rownames(counts)
}

if (is.null(length)){
if(level=="gene"){
if(is.character(length.gene)){
data(list=length.gene)
length.vec<-get(length.gene)
}else{
length.vec <- length.gene
}
}
if(level=="isoform"){
if(is.character(length.isoform)){
data(list=length.isoform)
length.vec<-get(length.isoform)
}else{
length.vec <- length.isoform
}
}
length<-matrix(rep(length.vec,each=ncol(counts)),ncol=ncol(counts),byrow=TRUE)
rownames(length)<-rownames(length.vec)
colnames(length)<-colnames(counts)
}else{
length<-length[,match(colnames(counts),colnames(length))]
}

genes<-genes[which(genes %in% rownames(length))]
counts<-counts[match(genes,rownames(counts)),,drop=FALSE]
length<-length[match(genes,rownames(length)),,drop=FALSE]

nlib <- ncol(counts)
ntags <- nrow(counts)
if (nlib > 0L && is.null(colnames(counts))){
colnames(counts) <- paste0("Sample", 1L:nlib)
}
if (ntags > 0L && is.null(rownames(counts))){
rownames(counts) <- 1L:ntags
}
if (!is.null(rpkms)){
if (nlib != ncol(rpkms)){
stop("'rpkms' and 'counts' must have equal number of columns")
}
if (ntags != nrow(rpkms)){
stop("'rpkms' and 'counts' must have equal number of rows")
}
}
if (is.null(lib.size)){
lib.size <- colSums(counts,na.rm=TRUE)
}
if (is.null(norm.factors)){
norm.factors <- rep(1, ncol(counts))
}
if (nlib != length(norm.factors)){
stop("Length of 'norm.factors' must equal number of columns in 'counts'")
}
if (is.null(group)){
group <- rep(1, ncol(counts))
}
group <- dropEmptyLevels(group)
if (nlib != length(group)){
stop("Length of 'group' must equal number of columns in 'counts'")
}
samples <- data.frame(group = group, lib.size = lib.size,norm.factors = norm.factors)
if (anyDuplicated(colnames(counts))){
message("Repeated column names found in count matrix")
row.names(samples) <- 1L:nlib
}else{
row.names(samples) <- colnames(counts)
}
x <- new("DGEList", list(counts = counts, length=length, samples = samples))
if (!is.null(genes)) {
genes <- as.data.frame(genes, stringsAsFactors = FALSE)
if (nrow(genes) != ntags){
stop("Counts and genes have different numbers of rows")
}
x$genes <- genes
}
if (remove.zeros) {
all.zeros <- rowSums(counts > 0, na.rm = TRUE) == 0
if (any(all.zeros)) {
x <- x[!all.zeros, ]
message("Removing ", sum(all.zeros), " rows with all zero counts")
}
}
x
}
