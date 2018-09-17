\name{DEGList}
\alias{DEGList}

\title{
DEGList Maker
}
\description{
Creates a DGEList object with various information including library sizes, gene length and others. 
}
\usage{
DEGList(counts = matrix(0, 0, 0), lib.size = colSums(counts, na.rm = TRUE), 
norm.factors = rep(1, ncol(counts)), group = rep(1, ncol(counts)), 
length = NULL, length.gene="length.gene.hg19",
length.isoform="length.isoform.hg19", rpkms = NULL, genes = NULL, 
remove.zeros = FALSE, level = "gene")
}

\arguments{
  \item{counts}{
numeric matrix of read counts.
}
  \item{lib.size}{
numeric vector giving sequence depth for each library, which is the total count by default.
}
  \item{norm.factors}{
numeric vector of normalization factors that modify the library sizes, which is 1 by default.
}
  \item{group}{
vector or factor giving the experimental group/condition for each sample/library.
}
  \item{length}{
numeric matrix of gene length, which could be gene expression effective length estimated by RSEM. 
}
  \item{length.gene}{
data frame containing annotation information for each gene, which is from the hg19 human assembly UCSC annotation by default.
}
  \item{length.isoform}{
 frame containing annotation of gene length for each gene isoform, which is from the hg19 human assembly UCSC annotation by default.
}
  \item{rpkms}{
numeric matrix of rpkms, which is NULL by default.
}
  \item{genes}{
data frame containing annotation information for each gene.
}
  \item{remove.zeros}{
logical, whether to remove rows that have 0 total count.
}
  \item{level}{
indicate the data is on gene level or isoform level.
}
}

\value{
DGEList object
}
