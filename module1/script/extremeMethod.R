# Extreme method. 
# From the paper below:
# Yu, S., Chakrabortty, A., Liao, K. P., Cai, T., 
# Ananthakrishnan, A. N., Gainer, V. S., Churchill, 
# S. E., Szolovits, P., Murphy, S. N., Kohane, I. S., 
# & Cai, T. (2017). Surrogate-assisted feature extraction 
# for high-throughput phenotyping. Journal of the American 
# Medical Informatics Association : JAMIA, 24(e1), e143–e149. 

extreme_method <- function(s, x, n_subset = 400, p_subset = 1, reps = 200,
                          u_bound = NULL, l_bound = NULL, q = 0.025) {
  
  # If only have 1 surrogate
  if (is.vector(s)) {
    s_label <- get_silver_label(s, u_bound, l_bound, q)
  # If have more than 1 surrogate
  } else {
    s_silver <- c()
    for (i in c(1:ncol(s))) {
      s_tmp <- get_silver_label(s[, i], u_bound[i], l_bound[i], q)
      s_silver <- cbind(s_silver, s_tmp)
    }
    s_silver <- rowMeans(s_silver)
    s_label <- rep(0.5, length(s_silver))
    s_label[s_silver < 0.5] <- 0
    s_label[s_silver > 0.5] <- 1
  }
  
  ind1 <- which(s_label == 1)
  ind0 <- which(s_label == 0)
  
  n1_subset <- floor(length(ind1)*p_subset)
  n0_subset <- floor(length(ind0)*p_subset)
  n_subset <- min(n1_subset, n0_subset, n_subset)
  
  # Weights for Adaptive lasso
  w <- rep(1, 2*n_subset)
  p <- ncol(x)
  lambda_init <- p / (2*n_subset)
  
  beta <- c()
  for (i in c(1:reps)) {
    ind <- c(sample(ind1, n_subset), sample(ind0, n_subset))
    beta_tmp <- fit_alasso_bic(s_label[ind], x[ind, ], family = "binomial", x_standardize = F)$beta_hat 
    beta_tmp <- beta_tmp[-1] #drop intercept
    beta_tmp <- as.integer(beta_tmp != 0)
    beta <- rbind(beta, beta_tmp)
  }
  
  beta_select <- which(colMeans(beta, na.rm = T) >= 0.5)
  
  return(list(beta_select = beta_select, beta_all = beta))
}

# Function to generate silver standard label. 
get_silver_label <- function(s, u_bound, l_bound, q) {
  
  # s: N x 1 surrogates
  # u_bound: upper bound 
  # l_bound: lower bound
  # q: quantile of s to get extreme surrogates

  if (is.null(u_bound)) {
    u_bound <- quantile(s, 1 - q)
  }
  
  if (is.null(l_bound)) {
    l_bound <- quantile(s, q)
  }
  
  s_silver <- rep(0.5, length(s))
  s_silver[s <= l_bound] <- 0
  s_silver[s >= u_bound] <- 1
  
  return(s_silver)
}
