% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LOFtest.R
\name{LOFtest}
\alias{LOFtest}
\title{Test a hypothesis that imposes many restrictions in a linear regression with many regressors and heteroskedasticity}
\usage{
LOFtest(linmod, R, q, nCores = 1, size = 0.05)
}
\arguments{
\item{linmod}{an object of class "lm" that stems from fitting the regression model \eqn{y = X\beta + \epsilon}.}

\item{R}{a matrix specifying the linearly independent restrictions on the parameters.
The number of rows in \code{R} is equal to the number of restrictions and the number of columns in \code{R} must be equal to the dimension of \eqn{\beta}.}

\item{q}{a vector of the hypothesized values for the linear combinations of the parameters \eqn{R\beta}.}

\item{nCores}{an optional argument for users who want to use parallel processing to speed up computations.}

\item{size}{the desired nominal size of the test. The standard is 5\%.}
}
\value{
\code{LOFtest} returns the realized value of the usual F-statistic,
the \code{p}-value corresponding to this realization when using the critical value function proposed in Anatolyev and Sølvsten (2020),
and the critical value for the specificed \code{size}.
}
\description{
Computes the leave-out F-test from Anatolyev and Sølvsten (2020) which performs a test of the null hypothesis \eqn{R\beta = q} in the linear regression model \eqn{y = X\beta + \epsilon}.
The test provides an adjustment to the usual critical value used with the F-test, and this adjustment makes the test robust to many restrictions and heteroskedasticity in the error terms.
The test retains asymptotic validity even when the number of regressors and restrictions are proportional to sample size.
}
\examples{
## An example of a regression with 640 observations, 512 regressors,
## and a null hypothesis that imposes 384 restrictions on the model.
\donttest{
set.seed(1)
X <- cbind(1, (0.5+runif(640))*matrix(exp(rnorm(640*511)), 640,511))
y <- X \%*\% c( -.5, rep(0.002,511) ) + rnorm(640)*.3*rowMeans(X)^2
R <- cbind( matrix(0, 384, 128), diag(384))
q <- rep(0.002, 384)
linmod <- lm(y~X-1)

LOFtest(linmod, R, q)
}
## The null is not rejected at the 5\\% level since the value of the
## F statistic is below the critical value returned by LOFtest.
## Relying on the usual critical value qf(.95,384,640-512)=1.279
## would incorrectly lead to a rejection of the null.
}
\references{
Anatolyev and Sølvsten (2020). \emph{Testing Many Restrictions Under Heteroskedasticity}. \url{https://arxiv.org/abs/2003.07320}
}
