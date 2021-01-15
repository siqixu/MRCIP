\name{MRCIP}
\alias{MRCIP}
\title{MRCIP: A robust Mendelian randomization method accounting for correlated and idiosyncratic pleiotropy}
\usage{
MRCIP(data,beta=0,mu_gamma=NA,s_gamma=NA,s_alpha=NA,rho=0,h=NA,h_start=0.01,h_step=0.01,tol_dw=1e-4,tol=1e-8,n_iter=1000,MVNtest="FALSE")
}
\description{
MRCIP is a Mendelian randomization approach for the inference of the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data, which accounts for the correlated pleiotropy and idiosyncratic pleiotropy.
}
\arguments{
  \item{data}{A matrix or data frame consists of four columns: the 1st (2nd) column contains the estimated genetic effects on the outcome (exposure); the 3rd (4th) column contains the estimated standard errors of the estimated genetic effects on the outcome (exposure).}
  \item{beta}{The initial value of the causal effect in the PRW-EM algorithm. The default value is zero.}
  \item{mu_gamma}{The initial value for the mean of the genetic effects on the exposure in the PRW-EM algorithm. If it is missing, the sample mean of the estimated genetic effects on the exposure will be used.}
  \item{s_gamma}{The initial value for the standard error of the genetic effects on the exposure in the PRW-EM algorithm. If it is missing, it will be calculated based on the sample variance of the estimated genetic effects on the exposure.}
  \item{s_alpha}{The initial value for the standard error of the pleiotropy effects on the outcome in the PRW-EM algorithm. If it is missing, it will be calculated based on the sample variance of the  estimated genetic effects on the outcome.}
  \item{rho}{The initial value of the correlated pleiotropy index (CPI) in the PRW-EM algorithm. The default value is zero.}
  \item{h}{The bandwidth used in the kernel estimation for the Pearson residuals. If it is missing, the value "sel_h" provided in the plot_h function will be used, which detects the possible abrupt change in the empirical mean down-weighting level over a grid of h values.}
   \item{h_start}{The start value in the grid search for the bandwidth used in the kernel estimation. The default value is 0.01. It will be ignored if h is provided.}
   \item{h_step}{The step size in the grid search for the bandwidth used in the kernel estimation. The default value is 0.01. It will be ignored if h is provided.}
   \item{tol_dw}{The grid search for the bandwidth stops when the difference in the mean downweighting level between two consecutive steps is less than tol_dw. It will be ignored if h is provided.}
  \item{tol}{The tolerance used in the PRW-EM algorithm. The PRW-EM algorithm stops when the difference in the weighted log-likelihood function between two consecutive iterations is less than tol. The default value is 1e-8.}
  \item{n_iter}{The maximum number of iterations in the PRW-EM algorithm. The default value is 1000.}
  \item{MVNtest}{If MVNtest="TRUE", then the chi-square QQ plot and the Mardia's test are provided for the multivariate normality test of the standardized residuals.}
}

\value{
\item{beta}{The estimated causal effect of the exposure on the outcome, the corresponding p value, the lower limit and upper limit of the 95 percent confidence interval}
\item{rho}{The estimated correlated pleiotropy index (CPI) as well as the corresponding p value}
\item{IVs with weights less than 0.5}{The indices and the weights of the genetic variants that have the weights less than 0.5}
\item{weights}{The weights of the genetic variants calculated by the Pearson residuals}
\item{Normality test for the standardized residuals}{The Mardia’s test for the multivariate skewness and kurtosis of the standardized residuals}
\item{standardized residuals}{The standardized residuals of the MRCIP model}
}

\examples{
  out = MRCIP(data=example)
  out$beta # the estimation and inference for the causal effect
  out$rho # the estimation and inference for correlated pleiotropy index (CPI)
}
