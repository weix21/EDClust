% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cluster.R
\name{NormalizeSC}
\alias{NormalizeSC}
\title{Normalization}
\usage{
NormalizeSC(X)
}
\arguments{
\item{X}{a count expression matrix.}
}
\value{
A normalized single-cell expression matrix.
}
\description{
This function serves as the data preprocessing step for SHARP
}
\examples{

XNorm <- NormalizeSC(X)


}
