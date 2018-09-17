\name{plotMA}
\alias{plotMA}

\title{
MA-plot Maker
}
\description{
Visualize the difference between two measurements by log ratio and mean average of data.
}
\usage{
plotMA(fit, coef = 2, h = 1, n = 3000, sort = "p", p = 0.05, colors = NULL, 
legend = NULL)
}

\arguments{
  \item{fit}{
the fit object from RoMA. 
}
  \item{coef}{
column number or column name of coefficient matrix which is of interest.
}
  \item{h}{
threshold line on log ratio
}
  \item{n}{
maximum number of significant genes will be shown
}
  \item{sort}{
values which will be ranked to select signficiant genes.
}
  \item{p}{
shreshold value
}
  \item{colors}{
colors will be display
}
  \item{legend}{
legend will be display
}
}
