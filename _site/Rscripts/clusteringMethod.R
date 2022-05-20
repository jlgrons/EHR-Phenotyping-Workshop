source("~/Documents/GitHub/heart-failure-phenotyping/Model development/modelFitting.R")

library(mclust)

clustering_method <- function(s, x, n_subset = 2000, p_subset = 1, reps = 200, resample = F) {
  
  pi_s <- guassian_clustering(s)$pi_s
  
  if (resample) {
    beta <- c()
    for (i in c(1:reps)) {
      ind <- sample(1:nrow(x), n_subset)
      
      beta_tmp <- fit_alasso_bic(pi_s[ind], x[ind, ], family = "gaussian")$beta_hat 
      beta_tmp <- beta_tmp[-1] #drop intercept
      beta_tmp <- as.integer(beta_tmp != 0)
      beta <- rbind(beta, beta_tmp)
      
      beta_select <- which(colMeans(beta, na.rm = T) >= 0.5)
    }
  } else {
    # beta <- fit_alasso_approx(pi_s, x, family = "gaussian")
    beta <- fit_alasso_bic(pi_s, x, family = "gaussian")$beta_hat 
    beta <- beta[-1]
    beta_select <- which(beta != 0)
  }

  return(list(beta_select = beta_select, beta = beta))
}

guassian_clustering <- function(s) {
  m <- Mclust(s, G = 2, verbose = F)
  
  # Group with the larger mean is the case group
  case_label <- as.double(which.max(colSums(m$parameters$mean)))
  
  # Save the predicted probability of surrogates as pi_s
  pi_s <- m$z[, case_label]
  
  s_group <- m$classification
  s_group[s_group != case_label] <- 0
  s_group[s_group == case_label] <- 1
  
  return(list(pi_s = pi_s, s_group = s_group))
}

