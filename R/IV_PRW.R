IV_PRW <- function(data,
                   beta=0, mu_gamma=NA, s_gamma=NA, s_alpha=NA, rho=0,
                   PRW_EM ="TRUE", h=NA, w=NA, rho_null="FALSE", tol=1e-8, n_iter=1000)  {
  Gamma_hat = data[,1]
  gamma_hat = data[,2]
  sj_Gamma = data[,3]
  sj_gamma = data[,4]
  Gamma_hat = sign(gamma_hat)*Gamma_hat
  gamma_hat = abs(gamma_hat)
  s2j_Gamma = sj_Gamma^2
  s2j_gamma = sj_gamma^2
  J=length(Gamma_hat)

  if(is.na(mu_gamma)){mu_gamma = mean(gamma_hat)}
  if(is.na(s_gamma)){s_gamma = sqrt(var(gamma_hat)-mean(s2j_gamma))}
  if(is.na(s_alpha)){s_alpha = sqrt(var(Gamma_hat)-beta^2*s_gamma^2-mean(s2j_Gamma))
  }else if (s_alpha==0){s_alpha = sqrt(var(Gamma_hat)-beta^2*s_gamma^2-mean(s2j_Gamma))}
  if(is.na(s_alpha)){s_alpha = 0.001} # in case that var(Gamma_hat)-beta^2*s_gamma^2-mean(s2j_Gamma)<0

  B = rbind(beta,1)
  C = rbind(1,0)
  Mu_eta = rbind(mu_gamma,0)
  S_eta = matrix(c(s_gamma^2,
                   rho*s_gamma*s_alpha,
                   rho*s_gamma*s_alpha,
                   s_alpha^2),nrow = 2)
  iv.S_eta = solve(S_eta)
  S = matrix(NA,2,2)
  S[1,1] = t(B)%*%S_eta%*%B + mean(s2j_Gamma)
  S[1,2]=S[2,1] = t(B)%*%S_eta%*%C
  S[2,2] = t(C)%*%S_eta%*%C + mean(s2j_gamma)
  Si = matrix(1,J,1) %*% as.vector(S) + cbind(s2j_Gamma,0,0,s2j_gamma)
  det_Si = Si[,1]*Si[,4]-Si[,2]^2
  iv.Si = (1/det_Si)* cbind(Si[,4],-Si[,2],-Si[,3],Si[,1])
  M = c(t(B)%*%Mu_eta,t(C)%*%Mu_eta)
  Mi = cbind(Gamma_hat-M[1], gamma_hat-M[2])

  maxlnL = diff_lnL = 1000
  k=0
  while (diff_lnL>tol)
  {
    if (PRW_EM == "TRUE"){
      # pearson residuals
      res = Mi[,1]/ sqrt(Si[,1])
      uniden1 = (matrix(res)%*%matrix(1,ncol = J) - matrix(1,nrow = J)%*%matrix(res, ncol = J))^2
      f.star = rowSums((1/sqrt(2*pi*h^2))*exp(-uniden1/(2*h^2)))/J
      m.star = dnorm(x = res)
      p.res = f.star/m.star - 1
      w = 1 - p.res^2/(p.res+2)^2
      w[is.nan(w)] = 0
    }

    # posterior mean and covariance
    Sj_eta = Eta_j = numeric()
    iv_Sj_eta = ( matrix(1,J,1)%*%as.vector(iv.S_eta) +
                    as.matrix(1/s2j_Gamma)%*%as.vector(B%*%t(B)) + as.matrix(1/s2j_gamma)%*%as.vector(C%*%t(C)) )
    det_Sj_eta = iv_Sj_eta[,1]*iv_Sj_eta[,4]-(iv_Sj_eta[,2])^2
    Sj_eta = (1/det_Sj_eta) * cbind(iv_Sj_eta[,4],-iv_Sj_eta[,2],-iv_Sj_eta[,3],iv_Sj_eta[,1])
    Eta_j_r = t( rep(1,J)%*%(t(Mu_eta)%*%iv.S_eta) +
                   (Gamma_hat/s2j_Gamma)%*%t(B) + (gamma_hat/s2j_gamma)%*%t(C) )
    Eta_j = rbind(colSums(t(Sj_eta[,1:2])*Eta_j_r),
                  colSums(t(Sj_eta[,3:4])*Eta_j_r))

    # beta
    nub = sum( w*(-Gamma_hat*Eta_j[1,] + Eta_j[2,]*Eta_j[1,] + Sj_eta[,2]) /s2j_Gamma )
    deb = sum( w*( Eta_j[1,]^2 + Sj_eta[,1] )/(-s2j_Gamma) )
    beta = nub/deb
    B = rbind(beta,1)

    #  Mu_eta
    mu_gamma = sum( Eta_j[1,]*w)/sum(w)
    mu_gamma = as.vector(mu_gamma)
    Mu_eta = rbind(mu_gamma,0)

    # S_eta
    if(rho_null=="TRUE"){
      s2_gamma = sum(w*(Eta_j[1,]-Mu_eta[1])^2 + Sj_eta[,1])/sum(w)
      s2_alpha = sum(w*(Eta_j[2,]-Mu_eta[2])^2 + Sj_eta[,4])/sum(w)
      S_eta = matrix(c(s2_gamma,0,0,s2_alpha), 2,2)
    }else{
      S_eta = ( rbind((Eta_j[1,]-Mu_eta[1])*w,(Eta_j[2,]-Mu_eta[2])*w)
                %*%cbind(Eta_j[1,]-Mu_eta[1],Eta_j[2,]-Mu_eta[2])
                + matrix(w%*%Sj_eta, 2, 2) )/ sum(w)
    }
    iv.S_eta = solve(S_eta)

    # lnL
    S = matrix(NA,2,2)
    S[1,1] = t(B)%*%S_eta%*%B
    S[1,2]=S[2,1] = t(B)%*%S_eta%*%C
    S[2,2] = t(C)%*%S_eta%*%C
    Si = matrix(1,J,1) %*% as.vector(S) + cbind(s2j_Gamma,0,0,s2j_gamma)
    det_Si = Si[,1]*Si[,4]-Si[,2]^2
    iv.Si = (1/det_Si)* cbind(Si[,4],-Si[,2],-Si[,3],Si[,1])
    M = c(t(B)%*%Mu_eta,t(C)%*%Mu_eta)
    Mi = cbind(Gamma_hat-M[1], gamma_hat-M[2])
    lnLi = log( (det_Si)^(-0.5) * exp(-0.5 *(Mi[,1]^2*iv.Si[,1]
                                             + 2*Mi[,1]*Mi[,2]*iv.Si[,2]
                                             + Mi[,2]^2*iv.Si[,4]) ))
    lnLi_w = lnLi*w
    lnLi_w[is.infinite(lnLi) & w<0.05] = 0
    diff_lnL = abs(sum(lnLi_w)-maxlnL)
    maxlnL = sum(lnLi_w)
    k=k+1
    if(k>n_iter){break}
  }

  rho = S_eta[1,2]/sqrt(S_eta[1,1])/sqrt(S_eta[2,2])
  s_gamma = sqrt(S_eta[1,1])
  s_alpha = sqrt(S_eta[2,2])
  dw = 1-sum(w)/length(w)
  outlier_id = data.frame("IVs_ID"=which(w<0.5), "weight"=round(w[w<0.5],3))

  std.resid = numeric()
  for(i in 1:J){
    tmp = eigen(matrix(iv.Si[i,],2))
    sqr_ivS = tmp$vectors%*%diag(sqrt(tmp$values))%*%t(tmp$vectors)
    std.resid = rbind(std.resid,Mi[i,]%*%sqr_ivS)
  }

  out = list("theta" = cbind(beta, mu_gamma, s_gamma, s_alpha, rho, maxlnL, k),
             "mean downweight" = cbind(h,dw),
             "weight" = w,
             "outlier_id"= outlier_id,
             "std.resid"= std.resid
  )
  return(out)
}



