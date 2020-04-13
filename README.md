## manyRegressors

R package for inference in models with heteroskedasticity and many regressors.

### To install the package, run the following code in R
```
install.packages("devtools") # Only needed if you do not have devtools already
devtools::install_github("mikkelsoelvsten/manyRegressors")
```
### Included functions

```LOFtest``` Test a hypothesis that imposes many restrictions in a linear regression with many regressors and heteroskedasticity.

```LOvcov``` Variance-covariance estimator that is robust to heteroskedasticity and many regressors.

### References

Kline, Saggio, and Sølvsten (2019). *Leave-out estimation of variance components*. [https://arxiv.org/abs/1806.01494](https://arxiv.org/abs/1806.01494)
