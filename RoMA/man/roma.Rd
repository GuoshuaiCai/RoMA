\name{roma}
\alias{roma}

\title{
The Main Function of RoMA
}
\description{
RoMA incorporates information from both mRNA abundance and raw counts by modeling RPKM (reads per kilobase per million), which represents the relative abundance of mRNA transcripts, and borrowing mean-variance dependency from CPM (counts per million) as a precision weight accounting for the variability in sequencing depth. LOWESS curve is fitted to capture the CPM mean-variance dependency and the inverse variance is taken into the modeling as the precision weight. Moderate-t statistics is performed by the eBayes function in limma.
}
\usage{
roma(counts, rpkms = NULL, design = NULL, lib.size = NULL, 
normalize.method = "none", span = 0.5, plot = FALSE, save.plot = FALSE, ...)
}

\arguments{
  \item{counts}{
numeric matrix containing raw counts or a DGEList object.
}
  \item{rpkms}{
numeric vector containing RPKMs for each sample. If NULL and counts is a DGEList then it can be obtained from counts$rpkms.
}
  \item{design}{
design matrix with rows corresponding to samples and columns to coefficients of study condition and covariates to be estimated.
}
  \item{lib.size}{
numeric vector containing total library sizes for each sample. If NULL and counts is a DGEList then library sizes are calculated from the columnwise counts totals.
}
  \item{normalize.method}{
normalization method to be applied to the logCPM values. Choices are as for the method argument of normalizeBetweenArrays of limma. As normalization is performed ahead in RoMA standard analysis, this argument is suggested to set as "none".
}
  \item{span}{
width of the lowess smoothing window as a proportion.
}
  \item{plot}{
logical, should a plot of the mean-variance trend be displayed?
}
  \item{save.plot}{
logical, should the coordinates and line of the plot be saved in the output?
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
