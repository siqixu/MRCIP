plot_h <- function(data,
                   beta=0, mu_gamma=NA, s_gamma=NA, s_alpha=NA, rho=0,
                   h_start=0.01, h_step=0.01, tol_dw=1e-4, tol=1e-8, n_iter=1000,
                   plot_h = "TRUE"){
  if(h_start<=0){h_start=0.01}
  if(h_step<=0){h_step=0.01}
  h = h_start
  theta = dw = w = numeric()
  d_dw = 1

  while (d_dw > tol_dw) {
    MRCIP_i = IV_PRW(data,
                     beta, mu_gamma, s_gamma, s_alpha, rho,
                     h=h, tol=tol, n_iter=n_iter)
    theta = rbind(theta,MRCIP_i[[1]])
    dw = rbind(dw,MRCIP_i[[2]])
    if(h!=h_start){
      d_dw = dw[nrow(dw)-1,"dw"] - dw[nrow(dw),"dw"]}
    h = h + h_step
  }
  if(plot_h=="TRUE"){
    plot(dw[,"h"],dw[,"dw"],
         xlab = "h", ylab = "mean downweighting level")
  }

  # detect the possible abrupt change in the mean downweighting level
  h = dw[,"h"]
  dw = dw[,"dw"]
  d_dw = round(dw - c(dw[-1],NA),2)
  d_dw_l = round((c(NA,NA,NA,dw)[1:length(dw)] - dw)/3, 2)
  d_dw_l[2] = round(dw[1]-dw[2],2)
  d_dw_l[3] = round((dw[1]-dw[3])/2,2)
  d_dw_r = round((c(dw[-1],NA) - c(dw[-(1:3)],NA,NA,NA))/2,2)
  tmp = cbind(h,dw,d_dw_l,d_dw,d_dw_r)
  sel = tmp[d_dw_l<d_dw,]
  sel = sel[sel[,"d_dw_r"] < sel[,"d_dw"],]
  sel = sel[is.na(sel[,1])==FALSE,]
  if(is.vector(sel)){sel = t(sel)}
  if(nrow(sel)==0){sel_h = max(h)
  }else if(nrow(sel)==1){sel_h = sel[,"h"]
  }else{
    sel = sel[which.max(round(2*sel[,"d_dw"]-sel[,"d_dw_l"]-sel[,"d_dw_r"],2)) ,]
    sel_h = sel["h"]
  }

  out = list("dw" = cbind(h,dw),
             "sel_h" = sel_h)
  return(out)
}


