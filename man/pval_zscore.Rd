% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pval_zscore.R
\name{pval_zscore}
\alias{pval_zscore}
\title{Compute p-values for a Z-score}
\usage{
pval_zscore(beta, sigma)
}
\arguments{
\item{beta}{the estimate}

\item{sigma}{estimate's estimated variance}
}
\value{
the p-value
}
\description{
Compute p-values for a Z-score assuming normal distribution of the z-score 
under the null Hypothesis H0
}
