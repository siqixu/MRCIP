# MRCIP：A robust Mendelian randomization method accounting for correlated and idiosyncratic pleiotropy
MRCIP is a Mendelian randomization approach for the inference of the causal effect of an exposure on an outcome of interest based on summary-level GWAS data, which accounts for the correlated pleiotropy and the idiosyncratic pleiotropy.
## Setup
The MRCIP package will use the R package "MVN" for the multivariate normality test of the model residuals. Therefore, make sure that this package has 
been installed, which is available on R CRAN https://CRAN.R-project.org/package=MVN. 

Use the following command in R to install the package:
```
install.packages(pkgs="MVN")  # install the "MVN" package
library(devtools)
install_github("siqixu/MRCIP",ref="main") # install the "MRCIP" package
```
## Usage
The MRCIP handles the idiosyncratic pleiotropy by using a weighting scheme, where the selection of bandwidth h in the kernel function is relevant to the trade-off between robustness and efficiency. Generally, smaller bandwidth h is more sensitive to the idiosyncratic pleiotropy, while larger h leads to higher efficiency. Here, we select h by monitoring the empirical mean down-weighting level over a grid of h values. 
 
The plot_h function provides a scatter plot of the mean down-weighting levels against a grid of bandwidth h, as well as a suggested h value. The suggested h value corresponds to a possible abrupt change in the scatter plot. When no abrupt change is detected in the plot, the h value corresponding to the smallest mean down-weighting level will be used as the suggested h value.
```
plot_h(data,beta=0,mu_gamma=NA,s_gamma=NA,s_alpha=NA,rho=0,h_start=0.01,h_step=0.01,tol_dw=1e-4,tol=1e-8,n_iter=1000,plot_h="TRUE")
```
Given h, the MRCIP function provides the estimations and inferences for the causal effect and the correlated pleiotropy.

```
MRCIP(data,beta=0,mu_gamma=NA,s_gamma=NA,s_alpha=NA,rho=0,h=NA,h_start=0.01,h_step=0.01,tol_dw=1e-4,tol=1e-8,n_iter=1000,MVNtest="FALSE")
```
## Example 
```
library(MRCIP)  # load the MRCIP package

# the input data should be a matrix or data frame consisting of four columns: the 1st(2nd) column contains the estimated genetic effects on the outcome (exposure); 
# the 3rd (4th) column contains the estimated standard errors of the estimated genetic effects on the outcome (exposure).

plot_h(data=example)   # select the bandwidth h

out = MRCIP(data=example) # analyze the data with the MRCIP method. By default, MRCIP uses the suggested h value in the "plot_h" function.

-----------Estimation and inference for the causal effect-----------
out$beta 
     beta_hat    pval 95%CI (ll) 95%CI (ul)
[1,]  0.03491 0.35649   -0.03929     0.1091

# beta_hat: the estimated causal effect is 0.03491.
# pval: the corresponding p value (H0: the causal effect is zero) is 0.0.35649.
# 95%CI: the 95% confidence interval is (-0.03929,0.1091). The 95%CI (ll) and 95%CI (ul) represent the lower limit and upper limit, respectively.

-----------Estimation and inference for correlated pleiotropy index (CPI)-------
out$rho
      rho_hat    pval
[1,] -0.16256 0.05359

# rho_hat: the estimated correlated pleiotropy index (CPI) is -0.16256.
# pval: the corresponding p value (H0: the correlated pleiotropy index is zero) is 0.05359.

```
## Reference
1. Xu S, Fung WK and Liu Z. MRCIP: A robust Mendelian randomization method accounting for correlated and idiosyncratic pleiotropy. Briefings in bioinformatics. 2021.


