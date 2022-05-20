library(glmnet)
library(glmpath)

# Function to calculate adaptive LASSO
fit_alasso_bic <- function(y, x, family = "binomial", x_standardize = F, offset = NULL, weights = NULL, 
                           init_lambda = NULL, init_opt = "ridge", bic_opt = "modified", bic_factor = 0.1) {
  
  if (is.null(offset)) {
    offset <- rep(0, length(y))
  }
  
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
  
  if (is.null(init_lambda)) {
    init_lambda <- ncol(x) / nrow(x)
  }
  
  # Step 1. Init the fit using ridge or ols estimator
  w_alasso <- init_weights(x, y, family = family, weights = weights, x_standardize = x_standardize, 
                           offset = offset, init_lambda = init_lambda, init_opt = init_opt)
  
  
  # Step 2. Choose the best parameter according to modified BIC 
  fit_tmp <- glmnet(x, y, family = family, weight = weights, offset = offset, 
                    penalty.factor = w_alasso, alpha = 1, 
                    nlambda = 100L, lambda.min.ratio = 1e-3)
  
  # Calculate bic. 
  # n <- length(y)
  # dev <- deviance(fit_tmp)
  # df <- sapply(predict(fit_tmp, type = "nonzero"), length) + 1L
  # 
  # bics <- bic(dev, n, df)
  n <- length(y)
  dev <- deviance(fit_tmp)
  nonzero <- sapply(predict(fit_tmp, type = "nonzero"), length) + 1L
  complexity <- nonzero * log(n)
  best_lambda <- fit_tmp$lambda[which.min(dev + complexity)]
  
  #best_lambda <- fit_tmp$lambda[which.min(bics)]
  beta_hat <- as.double(predict(fit_tmp, s = best_lambda, type = "coefficients"))
  
  # Alternative method. 
  # # Generate a sequence of x for presenting solution paths
  # x_transform <- x / vec2mat(w_alasso, nrow(x))
  # 
  # # The package glmpath provides the solution paths for a range of penalty parameters
  # tmpfit <- glmpath(x_transform, y, family = family, weight = weights, offset = offset, 
  #                  min.lambda = 0, standardize = x_standardize, lambda2 = init_lambda)
  # 
  # # Select a grid of lambda from min(lambda) to max(lambda)
  # min_lambda <- min(tmpfit$lambda)
  # max_lambda <- max(tmpfit$lambda)
  # lambda_seq <- seq(min_lambda, max_lambda, length.out = 500)
  # 
  # # Predict with the sequence of lambda 
  # # Return a matrix for betas, with lambda values as row names
  # beta_all <- predict.glmpath(tmpfit, s = lambda_seq, type = "coefficients", 
  #                             mode = "lambda", offset = offset)
  # beta_all <- beta_all / vec2mat(c(1, w_alasso), nrow(beta_all))
  # df_all <- apply(beta_all[, -1, drop = F] != 0, 1, sum)
  # 
  # # Return a matrix for log likelihood
  # loglik <- predict.glmpath(tmpfit, newx = x_rep, newy = y, s = lambda_seq, type = "loglik", 
  #                 mode = "lambda", offset = offset)
  # loglik <- apply(loglik, 2, sum)
  # bic_lambda <- modified_bic(loglik, length(y), bic_factor, df_all)
  # best_bic_ind <- which.min(bic_lambda) 
  # beta_hat <- beta_all[best_bic_ind, ]
  # lambda_hat <- lambda_seq[best_bic_ind]
  
  return(list(beta_hat = beta_hat, model = fit_tmp, lambda = best_lambda))
}

fit_alasso_approx <- function(y, x, family, x_standardize = F, offset = NULL, weights = NULL, 
                              bic_opt = "modified", bic_factor = 0.1) {
  if (is.null(offset)) {
    offset <- rep(0, length(y))
  }
  
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
  
  # Step 1. Init the fit using ridge or ols estimator.
  init_fit <- glm(y ~ x, offset = offset, weight = weights)
  init_coef <- as.vector(coef(init_fit))
  
  
  # Weights for alasso. 
  w_alasso <- 1 / abs(init_coef)
  
  #ddl is the second derivative of the log likelihood.
  ddl <- solve(summary(init_fit)$cov.unscaled)
  ddl_half <- svd(ddl)
  Xtilde <- ddl_half$u %*% diag(sqrt(ddl_half$d)) %*% t(ddl_half$u)
  Ytilde <- Xtilde %*% init_coef
  
  # Step 2: Fit alasso. Choose the best parameter according to BIC.
  fit_tmp <- glmnet(x = Xtilde, y = Ytilde, penalty.factor = w_alasso, alpha = 1, 
                    nlambda = 100L, lambda.min.ratio = 1e-3)
  
  # Calculate BIC. 
  n <- length(Ytilde)
  dev <- deviance(fit_tmp)
  df <- sapply(predict(fit_tmp, type = "nonzero"), length) + 1L
  
  if (bic_opt == "modified") {
    bics <- bic(dev, n, df)
  } else {
    bics <- modified_bic(dev, n, bic_factor, df)
  }
  
  # Choose the best parameter according to BIC.
  best_lambda <- fit_tmp$lambda[which.min(bics)]
  beta_hat <- as.double(predict(fit_tmp, s = best_lambda, type = "coefficients"))
}

# Init the fit using ridge or ols estimator. 
init_weights <- function(x, y, family, weights, x_standardize, offset, init_lambda, init_opt = "ridge") {
  if (init_opt == "ridge") {
    init_fit <- glmnet(x = x, y = y, family = family, weights = weights, offset = offset, 
                       alpha = 0, standardize = x_standardize, lambda = init_lambda)
  } else if (init_opt == "unpenalized") {
    init_fit <- glm(y ~ x, family = family, weights = weights, offset = offset)
  } else {
    stop("Please re-specify whether to use ridge or unpenalized GLM in initial fit.")
  }
  init_coef <- as.vector(coef(init_fit))
  
  # Weights for aLASSO
  w_alasso <- 1 / abs(init_coef[-1])
  return(w_alasso)
}

bic <- function(dev, n, df) {
  bic <- dev + df*log(n)
  return(bic)
}

# Modified BIC in (Minnier et al. 2011), which replace log(n) by min{n^0.1, log(n)}
# to avoid over shrinkage. 
modified_bic <- function(dev, n, bic_factor, df) {
  mod_bic <- dev + min(n^bic_factor, log(n))*df
  return(mod_bic)
}

# Vector to matrix. 
vec2mat <- function(vc, dm){
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}

# Compute the AUC. 

cal_AUC = function(data){
  dd = data[,1]; xx = data[,2]; n0 = sum(1-dd); n1 = sum(dd)
  x0 = xx[dd == 0]; x1 = xx[dd == 1]
  sum((sum.I(x0, "<=", x1)+sum.I(x0,"<",x1))/2)/(n0*n1)
}

# Computes sums efficiently based on ranks
sum.I <-function(yy, FUN, Yi, Vi=NULL, ties.method = "first") 
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  pos <- rank(c(yy,Yi), ties.method = ties.method)[1:length(yy)]-rank(yy,ties.method = ties.method)
  if (substring(FUN,2,2) == "=") pos <- length(Yi)-pos
  if (!is.null(Vi)) {
    if(substring(FUN,2,2) == "=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}


# Function to validate model.
validate_model <- function(train_y, test_y, train_y_hat, test_y_hat) {
  train_roc <- roc(train_y, train_y_hat)
  test_roc <- roc(test_y, test_y_hat)
  
  train_auc <- auc(train_roc)
  train_auc_ci <- ci(train_auc, method = "bootstrap", boot.n = 500, boot.stratified = F)
  train_auc_ci <- round(train_auc_ci, 3)
  train_auc <- train_auc_ci[2]
  train_auc_ci_len <- as.numeric(abs(train_auc_ci[3] - train_auc_ci[1]))
  train_auc_ci <- paste0("(", train_auc_ci[1], ",", train_auc_ci[3], ")")
  
  test_auc <- auc(test_roc)
  test_auc_ci <- ci(test_auc, method = "bootstrap", boot.n = 500, boot.stratified = F)
  test_auc_ci <- round(test_auc_ci, 3)
  test_auc <- test_auc_ci[2]
  test_auc_ci_len <- as.numeric(test_auc_ci[3] - test_auc_ci[1])
  test_auc_ci <- paste0("(", test_auc_ci[1], ",", test_auc_ci[3], ")")
  
  return(list(train_AUC = train_auc,
              train_CI = train_auc_ci,
              train_CI_length = train_auc_ci_len,
              test_AUC = test_auc, 
              test_CI = test_auc_ci,
              test_CI_length = test_auc_ci_len))
}
