\name{calcNormRPKMs}
\alias{calcNormRPKMs}

\title{
Calculation of Normalized RPKMs
}
\description{
Although RPKM/FPKM normalizes the RNA-seq data, it does not account for composition bias. Composition bias can be accounted for using TMM normalization by equaling the median of log-ratios of read counts between samples. RoMA use this function to normalize RPKMs by the TMM method.
}
\usage{
calcNormRPKMs(e)
}

\arguments{
  \item{e}{
The DGEList object with estimated library sizes and transcript lengths.
}
}

\value{
The DGEList object with calculated RPKMs in e$rpkms.
}
