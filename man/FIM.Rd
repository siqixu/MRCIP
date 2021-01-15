\name{FIM}
\alias{FIM}
\title{The estimated Fisher information matrix}
\usage{
FIM(sj_Gamma, sj_gamma, theta, w=NA)
}
\description{
This function provides an estimated Fisher information matrix for the parameters in the MRCIP model.
}
\arguments{
  \item{sj_Gamma}{A vector containing the estimated standard errors of the estimated genetic effects on the outcome}
  \item{sj_gamma}{A vector containing the estimated standard errors of the estimated genetic effects on the exposure}
  \item{theta}{A vector containing the estimated parameters in MRCIP}
  \item{w}{The weights of the genetic variants}
}
