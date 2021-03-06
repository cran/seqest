% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_gee.R
\name{print.seqGEE}
\alias{print.seqGEE}
\title{Print the results by the generalized estimating equations.}
\usage{
\method{print}{seqGEE}(x, ...)
}
\arguments{
\item{x}{A variable of type \code{seqGEE}}

\item{...}{Additional variables to be transferred to functions}
}
\value{
print.seqGEE returns several messages to show the conclusion.
}
\description{
\code{print.seqGEE} print the result of the logistic regression model used by
the method of adaptive shrinkage estimate.
}
\details{
This function is used to present results in a concise way. If we select
enough samples that satisfy the stopping criterion, then we show several
messages to report the conclusion including the length of fixed size
confidence set, the number of samples we choose, the value of coefficient and
the index of the non zero coefficient
}
