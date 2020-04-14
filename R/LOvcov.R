#' Variance-covariance estimator that is robust to heteroskedasticity and many regressors
#'
#' Computes a variance-covariance estimate for a set of linear contrasts \eqn{V\beta} in the regression model \eqn{y = X\beta + \epsilon} using
#' the leave-out estimator from Kline, Saggio, and Sølvsten (2019).  The leave-out estimator
#' provides a variance estimator which is asymptotically valid under many regressors and
#' unrestricted heteroskedasticity.
## #' For further details, see Kline, Saggio, and Sølvsten (2019), https://arxiv.org/abs/1806.01494.
#'
#' @param linmod an object of class "lm" that stems from fitting the regression model \eqn{y = X\beta + \epsilon}.
#' @param V a matrix specifying the linearly independent restrictions on the parameters.
#' The number of rows in \code{V} is equal to the number of contrasts and the number of columns in \code{V} must be equal to the dimension of \eqn{\beta}.
#'
#' @return \code{LOvcov} returns a variance-covariance matrix for the linear contrasts \eqn{V\beta}, using
#' the variance estimator proposed in Kline, Saggio, and Sølvsten (2019).
#' @references Kline, Saggio, and Sølvsten (2019). \emph{Leave-out estimation of variance components}. \url{https://arxiv.org/abs/1806.01494}
#' @examples
#' ## An example of a regression with 640 observations and 512 regressors
#' ## where the last three coefficients are of interest.
## #' \donttest{
#' set.seed(1)
#' X <- cbind(1, (0.5+runif(640))*matrix(exp(rnorm(640*511)), 640,511))
#' y <- X %*% c( -.5, rep(0.002,511) ) + rnorm(640)*.3*rowMeans(X)^2
#' V <- cbind( matrix(0, 3, 509), diag(3))
#' linmod <- lm(y ~ X-1)
#'
#' vcov <- LOvcov(linmod, V)
#'
#' ## Coefficients and standard errors of interest
#' coeff <- V %*% linmod$coefficients
#' cbind(coeff,sqrt(diag(vcov)))
#'
#' ## P-value for joint significance test
#' 1-pchisq( t(coeff) %*% solve(vcov) %*% coeff, df=3)
## #' }

#' @export

LOvcov <- function(linmod, V) {

  if (class(linmod) != "lm")
    stop("please supply an object of class 'lm'")

  if (ncol(V) != length(linmod$coefficients)) {
    stop("incompletely specified linear contrast matrix V; mismatch between the number of
         columns of V and the number of estimated parameters")
  }

  y <- linmod$fitted.values + linmod$residuals; y_<- y-mean(y)
  Pii <- stats::hatvalues(linmod)
  lo_sigma2hats <- y_ * linmod$residuals / (1.0 - Pii)
  lo_sigma2hats[Pii>=.99] <- (y_^2)[Pii>=.99]
  bread <- qr.Q(linmod$qr) %*% backsolve(qr.R(linmod$qr),t(V), transpose = TRUE)

  return( crossprod(bread,bread * lo_sigma2hats) )
}

