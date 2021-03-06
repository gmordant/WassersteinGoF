% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{generateGausWassersteinQuantile}
\alias{generateGausWassersteinQuantile}
\title{Generate the quantiles for an isotropic Gaussian.}
\usage{
generateGausWassersteinQuantile(quantiles = 1.2, nRep, samplesize,
  p = 2, B, dim, nCores = 0)
}
\arguments{
\item{quantiles}{a list of quantiles you want to return. By default, the function returns only a sample vector from the distribution under the null.}

\item{nRep}{takes the number of replications}

\item{samplesize}{sample size to be used}

\item{p}{1=W1, 2=W2}

\item{B}{the number of points to approximate the continuous distribution}

\item{dim}{the dimension}

\item{nCores}{number of cores for the parellelisation. By default, the number is the number of cores of your machine minus 1.}
}
\value{
list with sample vector from the distribution under the null and specified quantiles (if any)
}
\description{
Generate the quantiles for an isotropic Gaussian.
}
