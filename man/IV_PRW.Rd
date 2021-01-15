\name{IV_PRW}
\alias{IV_PRW}
\title{IV-PRW: A weighting scheme to handle the idiosyncratic pleiotropy}
\usage{
IV_PRW(data,beta=0,mu_gamma,s_gamma,s_alpha,rho=0,PRW_EM ="TRUE", h=NA, w=NA, rho_null="FALSE", tol=1e-8, n_iter=1000)
}
\description{
IV-PRW is a weighting scheme to handle the presence of idiosyncratic pleiotropy, where the model parameters are estimated by maximizing a weighted likelihood function with a PRW-EM algorithm
}
\arguments{
  \item{data}{A matrix or data frame consists of four columns: the 1st(2nd) column contains the estimated genetic effects on the outcome (exposure); the 3rd (4th) column contains the estimated standard errors of the estimated genetic effects on the outcome (exposure).}
  \item{beta}{The initial value of the causal effect in the PRW-EM algorithm. The default value is zero}
  \item{mu_gamma}{The initial value for the mean of the genetic effects on the exposure in the PRW-EM algorithm. If it is missing, the sample mean of the estimated genetic effects on the exposure will be used.}
  \item{s_gamma}{The initial value for the standard error of the genetic effects on the exposure in the PRW-EM algorithm. If it is missing, it will be calculated based on the sample variance of the estimated genetic effects on the exposure.}
  \item{s_alpha}{The initial value for the standard error of the pleiotropy effects on the outcome in the PRW-EM algorithm. If it is missing, it will be calculated based on the sample variance of the  estimated genetic effects on the outcome.}
  \item{rho}{The initial value of the correlated pleiotropy index (CPI) in the PRW-EM algorithm. The default value is zero.}
  \item{PRW_EM}{If it is "TRUE", then the PRW-EM algorithm will be used which contains the updating steps for the weights. Otherwise, the weights should be provided in the argument w.}
  \item{h}{The bandwidth used in the kernel estimation for the Pearson residuals. The value of h should be provided when the PRW_EM is set to be "TRUE".}
  \item{w}{The weights for the genetic variants. This argument will be ignored when PRW_EM="TRUE".}
  \item{rho_null}{If it is "TRUE", then the parameters in MRCIP will be estimated when rho is fixed at zero.}
  \item{tol}{The tolerance used in the PRW-EM algorithm. The default value is 1e-8.}
  \item{n_iter}{The maximum number of iterations in the PRW-EM algorithm. The default value is 1000.}
}


