roma <-
function(counts, rpkms=NULL, design = NULL, lib.size = NULL, normalize.method = "none", span = 0.5, plot = FALSE, save.plot = FALSE, ...){

if (is.null(rpkms)){
if (is(counts, "DGEList")) {
rpkms<-counts$rpkms
}
}
if (is.null(lib.size)){
if (is(counts, "DGEList")) {
lib.size<-counts$lib.size
}
}

v<-voom(counts,design=design,lib.size = lib.size,normalize.method=normalize.method,span=span,plot=plot,save.plot=save.plot, ...)

v$E<-log2(rpkms)

fit.lm <- lmFit(v,design=design)

fit <- eBayes(fit.lm)

fit
}
