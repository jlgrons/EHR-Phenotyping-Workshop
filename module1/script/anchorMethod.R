source("~/Documents/GitHub/heart-failure-phenotyping/Model development/modelFitting.R")

anchor_method <- function(s, x) {
  
  s_label <- as.integer(s >= 1)

  beta <- fit_alasso_bic(s_label, x, family = "binomial", x_standardize = F)$beta_hat
  beta <- beta[-1] #drop intercept
  beta <- as.double(which(beta != 0))

  return(beta)
}