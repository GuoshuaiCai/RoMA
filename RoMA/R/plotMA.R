plotMA <-
function(fit,coef=2,h=1,n=3000,sort = "p", p = 0.05,colors=NULL,legend=NULL){
if(sort=="p"){
sort<-"P.Value"
}
if(sort=="q"){
sort<-"adj.P.Val"
}
if(sort=="M"){
sort<-"logFC"
}
if(sort=="A"){
sort<-"AveExpr"
}

res<-topTable(fit,coef=coef,number=Inf)


if(is.null(colors)){
colors<-rep("black",nrow(res))
}
if(sort %in% c("P.Value","adj.P.Val")){
res<-res[order(res[,sort]),]
sig<-which(res[,sort]<=p)
sig<-sig[1:(if (length(sig)<n) length(sig) else n)]
colors[sig]<-"red"
}
if(sort %in% c("logFC","t")){
res<-res[order((-1)*abs(res[,sort])),]
sig<-which(abs(res[,sort])>=p)
sig<-sig[1:(if (length(sig)<n) length(sig) else n)]
colors[sig]<-"red"
}

plot(res$AveExpr,res$logFC,pch=20,cex=0.5,col=colors,xlab="A value", ylab="M value")
abline(h=0,v=0,lwd=1.5,col="black")
abline(h=c((-1)*h,h),lwd=rep(1.5,2), lty=2)

if(!is.null(legend)){
legend("topright",legend,fill=colors,bty="n")
}
}
