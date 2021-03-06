% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getWH_ord}
\alias{getWH_ord}
\title{Get the matrices W and H for the ordinal case}
\usage{
getWH_ord(data, beta)
}
\arguments{
\item{data}{A matrix containing the training samples which we will use in
the ordinal case.}

\item{beta}{A matrix contains the estimated coefficient. Note that the
beta_mat is a n * k matrix which n is the number of the explanatory
variables and k+1 is the number of categories}
}
\value{
a list contains several components including the variance-covariance
  matrix, minimum eigenvalue, W and H.
}
\description{
Get the matrices W and H using the Rcpp package for the ordinal case
}
\details{
getWH_ord uses the current training data and the estimated coefficient under
the ordinal case to obtain the matrices W and H to further get the
variance-covariance matrix and minimum eigenvalue. Note that using the Rcpp
package can significantly reduce the time of operation and get conclusions
faster.
}
