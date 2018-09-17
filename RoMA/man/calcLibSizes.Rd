\name{calcLibSizes}
\alias{calcLibSizes}

\title{
Calculation of Library Sizes
}
\description{
The effective library size is estimated by norm factor * total counts, where norm factor can be estimated by TMM.
}
\usage{
calcLibSizes(e)
}

\arguments{
  \item{e}{
The DGEList object with estimated normalization factors.
}
}

\value{
The DGEList object with estimated effective library sizes in e$lib.size.
}
