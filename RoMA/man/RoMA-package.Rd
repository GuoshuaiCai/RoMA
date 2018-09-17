\name{RoMA-package}
\alias{RoMA-package}
\alias{RoMA}
\docType{package}
\title{
\packageTitle{RoMA}
}
\description{
\packageDescription{RoMA}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{RoMA}
\packageIndices{RoMA}

}
\author{
\packageAuthor{RoMA}

Maintainer: \packageMaintainer{RoMA}
}

\keyword{ package }

\examples{

data(x)

e<-DEGList(counts=x)

e <- calcNormFactors(e,method="TMM")

e<-calcLibSizes(e)
e<-calcNormRPKMs(e)

des<-factor(c(rep("BR",7),rep("UHR",7)))

des<-as.matrix(cbind(1,as.numeric(des)))

fit<-roma(e, design = des, lib.size = NULL, normalize.method = "none", 
span = 0.5, plot = TRUE)

plotMA(fit,coef=2,h=1,n=3000,sort="p",p=0.05)

res<-topTable(fit,coef=2,n=Inf,sort = "p", p = 0.05)

res<-topTable(fit,coef=2,n=Inf)
}
