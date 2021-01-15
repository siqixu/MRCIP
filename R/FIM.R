FIM <- function(sj_Gamma, sj_gamma, theta, w=NA){

  beta = theta[1]
  mu_gamma = theta[2]
  s_gamma = theta[3]
  s_alpha = theta[4]
  rho = theta[5]

  cov = s_gamma * s_alpha * rho
  du1 = c(mu_gamma,0)
  du2 = c(beta,1)
  du3 = du4 = du5 = c(0,0)
  ds1 = matrix(c(2*beta*s_gamma^2 + 2*cov, s_gamma^2, s_gamma^2, 0),2,2)
  ds2 = matrix(c(0,0,0,0),2,2)
  ds3 = matrix(2*s_gamma*c(beta^2,beta,beta,1),2,2)
  ds4 = matrix(c(2*s_alpha,0,0,0),2,2)
  ds5 = matrix(c(2*beta,1,1,0),2,2)
  D = matrix(c(beta^2*s_gamma^2 + s_alpha^2 + 2*beta*cov ,
               beta*s_gamma^2 + cov,
               beta*s_gamma^2 + cov,
               s_gamma^2),2,2)
  FIM = matrix(0,5,5)

  for(j in 1:length(sj_Gamma)){
    Dj = D +  matrix(c(sj_Gamma[j]^2, 0, 0, sj_gamma[j]^2),2,2)
    iv.Dj = solve(Dj)
    FIM_j = matrix(0,5,5)
    FIM_j[1,1] = du1%*%iv.Dj%*%du1 + 0.5*sum(diag(iv.Dj%*%ds1%*%iv.Dj%*%ds1))
    FIM_j[1,2] = du1%*%iv.Dj%*%du2 + 0.5*sum(diag(iv.Dj%*%ds1%*%iv.Dj%*%ds2))
    FIM_j[1,3] = du1%*%iv.Dj%*%du3 + 0.5*sum(diag(iv.Dj%*%ds1%*%iv.Dj%*%ds3))
    FIM_j[1,4] = du1%*%iv.Dj%*%du4 + 0.5*sum(diag(iv.Dj%*%ds1%*%iv.Dj%*%ds4))
    FIM_j[1,5] = du1%*%iv.Dj%*%du5 + 0.5*sum(diag(iv.Dj%*%ds1%*%iv.Dj%*%ds5))
    FIM_j[2,2] = du2%*%iv.Dj%*%du2 + 0.5*sum(diag(iv.Dj%*%ds2%*%iv.Dj%*%ds2))
    FIM_j[2,3] = du2%*%iv.Dj%*%du3 + 0.5*sum(diag(iv.Dj%*%ds2%*%iv.Dj%*%ds3))
    FIM_j[2,4] = du2%*%iv.Dj%*%du4 + 0.5*sum(diag(iv.Dj%*%ds2%*%iv.Dj%*%ds4))
    FIM_j[2,5] = du2%*%iv.Dj%*%du5 + 0.5*sum(diag(iv.Dj%*%ds2%*%iv.Dj%*%ds5))
    FIM_j[3,3] = du3%*%iv.Dj%*%du3 + 0.5*sum(diag(iv.Dj%*%ds3%*%iv.Dj%*%ds3))
    FIM_j[3,4] = du3%*%iv.Dj%*%du4 + 0.5*sum(diag(iv.Dj%*%ds3%*%iv.Dj%*%ds4))
    FIM_j[3,5] = du3%*%iv.Dj%*%du5 + 0.5*sum(diag(iv.Dj%*%ds3%*%iv.Dj%*%ds5))
    FIM_j[4,4] = du4%*%iv.Dj%*%du4 + 0.5*sum(diag(iv.Dj%*%ds4%*%iv.Dj%*%ds4))
    FIM_j[4,5] = du4%*%iv.Dj%*%du5 + 0.5*sum(diag(iv.Dj%*%ds4%*%iv.Dj%*%ds5))
    FIM_j[5,5] = du5%*%iv.Dj%*%du5 + 0.5*sum(diag(iv.Dj%*%ds5%*%iv.Dj%*%ds5))
    FIM_j = FIM_j + t(FIM_j) - diag(diag(FIM_j))
    FIM = FIM + w[j]*FIM_j
  }
  return(FIM)
}




