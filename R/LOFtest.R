#' Test a hypothesis that imposes many restrictions in a linear regression with many regressors and heteroskedasticity
#'
#' Computes the leave-out F-test from Anatolyev and Sølvsten (2020) which performs a test of the null hypothesis \eqn{R\beta = q} in the linear regression model \eqn{y = X\beta + \epsilon}.
#' The test provides an adjustment to the usual critical value used with the F-test, and this adjustment makes the test robust to many restrictions and heteroskedasticity in the error terms.
#' The test retains asymptotic validity even when the number of regressors and restrictions are proportional to sample size.
## #' For further details, see Anatolyev and Sølvsten (2020), https://arxiv.org/abs/2003.07320.
#'
#' @param linmod an object of class "lm" that stems from fitting the regression model \eqn{y = X\beta + \epsilon}.
#' @param R a matrix specifying the linearly independent restrictions on the parameters.
#' The number of rows in \code{R} is equal to the number of restrictions and the number of columns in \code{R} must be equal to the dimension of \eqn{\beta}.
#' @param q a vector of the hypothesized values for the linear combinations of the parameters \eqn{R\beta}.
#' @param nCores an optional argument for users who want to use parallel processing to speed up computations.
## #' @param frac the fraction of observations to compute the leave-three-out estimates for. A smaller value decreases computation time but decreases the accuary of the critical value.
#' @param size the desired nominal size of the test. The standard is 5\%.
#'
#' @return \code{LOFtest} returns the realized value of the usual F-statistic,
#' the \code{p}-value corresponding to this realization when using the critical value function proposed in Anatolyev and Sølvsten (2020),
#' and the critical value for the specificed \code{size}.
#'
#' @references Anatolyev and Sølvsten (2020). \emph{Testing Many Restrictions Under Heteroskedasticity}. \url{https://arxiv.org/abs/2003.07320}
#'
#' @examples
#' ## An example of a regression with 640 observations, 512 regressors,
#' ## and a null hypothesis that imposes 384 restrictions on the model.
#' \donttest{
#' set.seed(1)
#' X <- cbind(1, (0.5+runif(640))*matrix(exp(rnorm(640*511)), 640,511))
#' y <- X %*% c( -.5, rep(0.002,511) ) + rnorm(640)*.3*rowMeans(X)^2
#' R <- cbind( matrix(0, 384, 128), diag(384))
#' q <- rep(0.002, 384)
#' linmod <- lm(y~X-1)
#'
#' LOFtest(linmod, R, q)
#' }
#' ## The null is not rejected at the 5\% level since the value of the
#' ## F statistic is below the critical value returned by LOFtest.
#' ## Relying on the usual critical value qf(.95,384,640-512)=1.279
#' ## would incorrectly lead to a rejection of the null.

#' @export

LOFtest <- function( linmod, R, q, nCores = 1, size = 0.05 ){

  if (class(linmod) != "lm")
    stop("please supply an object of class 'lm'")

  #Dimensions
  n <- length(linmod$residuals); m <- length(linmod$coefficients); r <- nrow(R)

  if (ncol(R) != m)
    stop("incompletely specified restriction matrix R; mismatch between the number of
         columns of R and the number of regression coefficient")

  if (length(q) != r)
    stop("mismatch between number of restrictions and the values for the restrictions
         to take on; check that the length of q matches the number of rows of R")

  #Projection matrices
  Q <- qr.Q(linmod$qr)
  M <- diag(n)-tcrossprod(Q); dM <- diag(M)
  if ( min(dM) < 0.001 )
    stop("The supplied design matrix loses full rank when dropping a single observation,
          this means that one (or more) coefficients are estimated based on only one observation with leverage of one.
          Consider removing observations with leverage of one and dropping from the model the corresponding parameters
          that are estimated solely from these observations")
  inv <- backsolve(qr.R(linmod$qr),t(R), transpose = TRUE)
  qr.inv <- qr(inv)
  bread <- Q %*% qr.Q(qr.inv)
  B <- tcrossprod( bread ); dB <- diag(B)

  #Numerator of F-statistic
  Fnum <- sum( backsolve( qr.R(qr.inv), R %*% linmod$coefficients - q, transpose = TRUE )^2 )
  #Residuals and leave-one-out variance estimator
  y <- linmod$fitted.values + linmod$residuals; y_ <- y-mean(y)
  e <- linmod$residuals;  hs1 <- e*y_/dM
  #Estimated mean for numerator of F-statistic
  hE <- sum(dB*hs1)

  #Unbiased variance-covariance matrix for rotated linear combinations
  Omega = crossprod(bread, bread * hs1)

  #Variance estimator matrices
  MdBdM <- M * dB/dM
  V <- -2*Matrix::skewpart(MdBdM)
  UV <- 2*( B - Matrix::symmpart(MdBdM) )^2 - V^2
  Vy <- V %*% y_
  #Biased variance estimator
  Var <- as.numeric( crossprod(hs1,UV %*% hs1) + sum( Vy^2 * hs1 ) )

  #clean-up before computing unbiased variance estimator for leave-out F-test
  rm(qr.inv,bread,B,dB,MdBdM)

  #frac < 1 provides a randomized approximation of the leave-three-out variance estimator
  frac <- 1
  #Randomly chosen elements to compute:
  Sel <- Matrix::Matrix(FALSE, nrow = n, ncol = n)
  Sel[upper.tri(Sel, diag = FALSE)] <- (stats::rbinom(n*(n-1)/2, 1, frac) >= 1)
  #Compute only important elements
  # sUV = rowSums(UV^2)
  # sV = rowSums(V^2)
  # Q = Q & ( (UV^2 > outer(sUV,sUV,FUN=function(x,y) x+y)/2/sqrt(n)) | (V^2 > outer(sV,sV,FUN=function(x,y) x+y)/2/sqrt(n) ) )

  #Threshold for treating D2 and D3 as zero
  eps <- 0.01

  `%mode%` <- foreach::`%do%`; i<-1
  if(nCores>1){cl <- parallel::makeCluster(nCores); doParallel::registerDoParallel(cl); `%mode%` <- foreach::`%dopar%`}
  Sp <- foreach::foreach(i=1:n, .combine=rbind, .packages = c("Matrix")) %mode% {
    #Loop over the choosen pairs. Calculate relevant variance terms for both i and j.
    ji <- (1:n)[Sel[i,]]
    if( length(ji)>0 ){
      #Vectors to collect results in
      sp <- 0; seciBias <- rep(0,length(ji)); secjBias <- rep(0,length(ji)); loc <- 1
      #Extract relevant numbers and vectors for i - may be approximated using JLA
      Mi <- M[,i]; Mii <- Mi[i]; MiSq <- Mi^2
      D2i <- Mii * dM - MiSq
      e2i <- (dM * e[i] - Mi * e)/D2i; e2i[D2i < eps^2] <- 0
      e2ki <- (Mii * e - Mi * e[i])/D2i; e2ki[D2i < eps^2] <- 0
      Vyi <- V[,i] * y_; yi <- y_[i]

      for( j in ji){
        #Extract relevant numbers and vectors for j - may be approximated using JLA
        Mj <- M[,j]; Mjj <- Mj[j]; MjSq <- Mj^2;
        D2j <- Mjj * dM - MjSq #D2[,j];
        e2j <- (dM * e[j] - Mj * e)/D2j; e2j[D2j < eps^2] <- 0
        e2kj <- (Mjj * e - Mj * e[j])/D2j; e2kj[D2j < eps^2] <- 0
        Vyj <- V[,j] * y_; yj <- y_[j]
        #Extract relevant joint numbers - may be approximated using JLA
        UVij <- UV[i,j]
        #Computing Dijk
        D3 <- Mii*D2j - Mjj * MiSq - dM * Mi[j]^2 + 2 * Mi[j] * Mi * Mj
        #Triples of observations where leave-three-out fails due to i or j
        G3i <- (D3 + (D2j<eps^3) * D2j[i] * D2i / eps < eps^3)
        G3j <- (D3 + (D2i<eps^3) * D2i[j] * D2j / eps < eps^3)
        #Leave-three-out residuals
        e3i <- (e[i] - Mi[j]*e2j - Mi*e2kj)/ D3 * D2j
        e3i[ D3 < eps^3 ] <- e2i[D3 < eps^3]
        e3i[ G3i == 1 ] <- 0
        e3j <- (e[j] - Mj[i]*e2i - Mj*e2ki)/ D3 * D2i
        e3j[ D3 < eps^3 ] <- e2j[D3 < eps^3]
        e3j[ G3j == 1 ] <- 0
        #Second component of variance estimator minus the biased part in Var
        seciBias[loc] <- yi^2 * Vyi[j] * sum( Vyi[G3i] )
        sp <- sp + yi * Vyi[j] * sum( Vyi * e3i ) - Vyi[j] * Vy[i] * hs1[i]
        secjBias[loc] <- yj^2 * Vyj[i] * sum( Vyj[G3j] )
        sp <- sp + yj * Vyj[i] * sum( Vyj * e3j ) - Vyj[i] * Vy[j] * hs1[j]
        #Variance product estimator weighted by UVij minus the biased part in Var
        if( sum(G3i*G3j) == 0 ){
          sp <- sp + 2 * UVij * ( yi * yj / D2i[j] * sum( (Mjj*Mi - Mi[j]*Mj) * y_ * e3j ) - hs1[i] * hs1[j])
        } else if ( D2i[j] > eps^2 ){
          #Upward biased estimator when the same leave-three-out failure is due to i and j
          sp <- sp + 2 * UVij * ( (UVij >0)* ( .5 * yi^2 * yj * e2j[i] + .5 * yj^2 * yi * e2i[j] ) - hs1[i] * hs1[j])
        } else {
          sp <- sp + 2 * UVij * ( (UVij >0)* yi^2 * yj^2 - hs1[i] * hs1[j])
        }
        loc <- loc +1
      }
      rbind(c(0,sp),c(i,sum(seciBias)),cbind(ji[secjBias!=0],secjBias[secjBias!=0]))
    } else
      c(0,0)
  }
  if(nCores>1){parallel::stopCluster(cl); rm(cl)}
  #Final variance estimate. Unbiased if leave-three-out holds everywhere.
  #Scaling by inverse of probability that Qij is non-zero.
  pos <- tapply(Sp[Sp[,1]>0,][,2], Sp[Sp[,1]>0,][,1], FUN=sum)
  LOFVar <- Var + (sum(Sp[Sp[,1]==0,2]) + sum( pos*(pos >0))) / frac
  #Test statistic uses leave-three-out variance estimate if positive
  #and a guaranteed positive and upward biased alternative otherwise
  if(LOFVar <= 0){    LOFVar <- sum( UV * (UV >= 0) * tcrossprod(y_^2) ) + sum( Vy^2 * y_^2  )  }


  #Dimension robust inference
  hatw <- pmax( eigen( Omega, only.values = TRUE)$values, 0 )
  hatw <- hatw / sum( hatw )
  barF <- replicate( 49999, crossprod( stats::rchisq(r,1), hatw) / stats::rchisq(1,n-m)*(n-m) )


  #Where to store results
  res <- list(Fstat=numeric(), pval=numeric(), critval=numeric() )
  res$Fstat <- Fnum/r/sum(e^2)*(n-m)
  res$pval <- mean( (hE + sqrt(LOFVar)*( barF - 1 )/sqrt( 2*sum(hatw^2) + 2/(n-m) ) )/r/sum(e^2)*(n-m) >= res$Fstat )
  res$critval <- (hE + sqrt(LOFVar)*( stats::quantile(barF,1-size) - 1)/sqrt( 2*sum(hatw^2) + 2/(n-m) ) )/r/sum(e^2)*(n-m)

  res
}
