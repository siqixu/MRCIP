\name{plot_h}
\alias{plot_h}
\title{A scatter plot of the mean downweighting level against the bandwidth in the kernel estimation}
\usage{
plot_h(data,beta=0,mu_gamma=NA,s_gamma=NA,s_alpha=NA,rho=0,h_start=0.01,h_step=0.01,tol_dw=1e-4,tol=1e-8,n_iter=1000,plot_h="TRUE")
}
\description{
This function provides a scatter plot of the mean downweighting level against the bandwidth in the kernel estimation}
\arguments{
 \item{data}{A matrix or data frame consists of four columns: the 1st(2nd) column contains the estimated genetic effects on the outcome (exposure); the 3rd (4th) column contains the estimated standard errors of the estimated genetic effects on the outcome (exposure).}
  \item{beta}{The initial value of the causal effect in the PRW-EM algorithm. The default value is zero}
  \item{mu_gamma}{The initial value for the mean of the genetic effects on the exposure in the PRW-EM algorithm. If it is missing, the sample mean of the estimated genetic effects on the exposure will be used.}
  \item{s_gamma}{The initial value for the standard error of the genetic effects on the exposure in the PRW-EM algorithm. If it is missing, it will be calculated based on the sample variance of the estimated genetic effects on the exposure.}
  \item{s_alpha}{The initial value for the standard error of the pleiotropy effects on the outcome in the PRW-EM algorithm. If it is missing, it will be calculated based on the sample variance of the  estimated genetic effects on the outcome.}
  \item{rho}{The initial value of the correlated pleiotropy index (CPI) in the PRW-EM algorithm. The default value is zero.}
  \item{h_start}{The start value in the grid search for the bandwidth used in the kernel estimation. The default value is 0.01.}
\item{h_step}{The step size in the grid search for the bandwidth used in the kernel estimation. The default value is 0.01.}
\item{tol_dw}{The grid search for the bandwidth stops when the difference in the mean downweighting between two consecutive steps is less than  tol_dw.}
  \item{tol}{The tolerance used in the PRW-EM algorithm. The PRW-EM algorithm stops when the difference in the weighted log-likelihood function between two consecutive iterations is less than tol. The default value is 1e-8.}
  \item{n_iter}{The maximum number of iterations in the PRW-EM algorithm. The default value is 1000.}
  \item{plot_h}{If it is "TRUE", then a plot of the mean downweighting level against the bandwidth will be produced.}
}

\value{
\item{dw}{A matrix containing two columns. The 1st column is the bandwidth used in the kernel estimation. The 2nd column is the corresponding empirical mean downweighting level given the bandwidth in the 1st column.}
\item{sel_h}{The possible abrupt change in the empirical mean down-weighting level over a grid of bandwidth values. If no abrupt change is detected, it will be set to the maximum value among the grid of bandwidth values}
}

\examples{
  plot_h(data=example)
}
