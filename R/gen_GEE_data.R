#' @title Generate the datasets with clusters
#'
#' @description
#' \code{gen_GEE_data} generates the clustered data used for the generalized
#' estimating equations with sequential method.
#'
#' @details
#' The gen_GEE_data function is used to generate data. We can get data from two
#' different distributions, corresponding to continuous and discrete cases. In
#' the continuous case, the covariates vector x is created from a multivariate
#' normal distribution with mean 0 and an AR(1) correlation matrix with
#' autocorrelation coefficient and marginal variance. The value of
#' autocorrelation coefficient and marginal variance are two arguments which we
#' need specified. Then, the response y is generated by the equation: y = wx + e
#' where the random error vector e follows a normal distribution with mean 0 and
#' three different covariance structures with corresponding dimensional numbers.
#' These three covariance matrices are the identity matrix, the exchangeable,
#' and the AR(1) autoregressive correlation structure. In the discrete case, we
#' use a logistic model. The covariates vectors x is the same as the continuous
#' case. The binary response vector for each cluster has  an AR(1) correlation
#' structure with correlation coefficient alpha, and the marginal expectation u
#' satisfies the following equation: logit(u) = wx
#' @param numClusters A numeric number represents the number of clusters we will
#'   generated. Note that each cluster has several similar subjects. It should
#'   be a integer.
#' @param clusterSize  A numeric number specifying the number of subjects in
#'   each cluster. The subject in the same cluster is highly correlated to each
#'   other which can be regarded as the longitudinal data.
#' @param clusterRho A numeric parameter in correlation structure for the
#'   clusters. It will be ignored when responseCorstr is independence.
#' @param clusterCorstr A character string specifying the correlation structure
#'   for the clusters. Allowed structures are: "independence", "exchangeable"
#'   and "ar1".
#' @param beta A nummeric vector denotes the true parameter in GEE model.
#' @param family The type of response data, matching one of 'gaussian()' or
#'   'binomial()'. The 'gaussian()' corresponds to the continuous case and
#'   'binomial' corresponds to the discrete case.
#' @param intercept A logical value indicating whether to add intercept term.
#'   The default value is TRUE.
#' @param xCorstr A character string specifying the correlation structure for
#'   the covariate. The default value is 'ar1'.
#' @param xCorRho A numeric parameter indicating the correlation coefficient in
#'   covariables. It does something similar to what the argument clusterRho
#'   does. The default value is 0.5.
#' @param xVariance A numeric number specifying the marginal variance in the
#'   correlation matrix in one clusters. The default value is 0.2.
#' @export
#' @return a list containing the following components
#' \item{x}{the covariate matrices. Note that the number of rows is numClusters
#' *  clusterSize and the number of columns is the length of beta + 1 if
#' intercept is TRUE.}
#' \item{y}{the response data which has the same number of rows to x}
#' \item{clusterID}{the id for each sample. Note that the subjects in the same
#' cluster will have identical id. }
#'
#' @references {
#' Chen, Z., Wang, Z., & Chang, Y. I. (2019). Sequential adaptive variables and
#' subject selection for GEE methods. \emph{Biometrics}. doi:10.1111/biom.13160
#' }
#'
#' @seealso{
#'    \code{\link{gen_multi_data}} for categorical and ordinal case
#'
#'    \code{\link{gen_bin_data}} for binary classification case.
#'}
#'
#' @examples
#' initialSampleSize <-  75
#' clusterSize <-  5
#' responseCorstr <-  "ar1"
#' responseCorRho <-  0.3
#' response <-  gaussian()
#' beta0 <-  c(1, -1.1, 1.5, -2, rep(0, 50))
#' xVariance <-  0.2
#' xCorRho <-  0.5
#' xCorstr <-  "ar1"
#' data <- gen_GEE_data(numClusters = initialSampleSize,
#'                      clusterSize = clusterSize,
#'                      clusterCorstr = responseCorstr,
#'                      clusterRho = responseCorRho,
#'                      beta = beta0,
#'                      family = response,
#'                      intercept = TRUE,
#'                      xVariance = xVariance,
#'                      xCorstr = xCorstr,
#'                      xCorRho = xCorRho)

gen_GEE_data <- function(numClusters,
                    clusterSize,
                    clusterRho,
                    clusterCorstr,
                    beta,
                    family,
                    intercept = TRUE,
                    xCorstr = "ar1",
                    xCorRho = 0.5,
                    xVariance = 0.2
) {
  p <- length(beta)
  n <- numClusters * clusterSize

  ## gen X
  if (intercept) p <- p - 1
  zcor <- genCorMat(xCorstr, xCorRho, p)
  diag(zcor) <- rep(1, p)
  z <- mvtnorm::rmvnorm(n, rep(0, p), zcor * xVariance)
  if (intercept) x <- cbind(1, z) else x <- z

  eta <- x %*% beta

  ## gen y
  if (family$family == "binomial") {
    # family <- get(family)
    mu <- family$linkinv(eta)
    genBinf <- genBin(clusterCorstr)
    y <- c(apply(mu, 2, genBinf, clusterRho))
  } else { ## Gaussian
    cormat <- genCorMat(clusterCorstr, clusterRho, clusterSize)
    e <- c(t(mvtnorm::rmvnorm(numClusters, rep(0, clusterSize), cormat)))
    y <- eta + e
  }
  return(list(x=x, y=y, clusterID=rep(1:numClusters, each=clusterSize)))
}
