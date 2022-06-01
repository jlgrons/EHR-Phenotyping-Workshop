library(mclust)
library(glmnet)
library(glmpath)

#################################################################
################### Feature selection methods ###################
#################################################################

# Absolute rank correlation method.
# From the paper below:
# Sheng Yu, Katherine P Liao, Stanley Y Shaw, Vivian S Gainer,
# Susanne E Churchill, Peter Szolovits, Shawn N Murphy,
# Isaac S. Kohane, Tianxi Cai,
# Toward high-throughput phenotyping: unbiased automated feature
# extraction and selection from knowledge sources, Journal of the
# American Medical Informatics Association, Volume 22, Issue 5,
# September 2015, Pages 993–1000.

rankCor <- function(s, x, threshold = 0.15) {

  # s: N x 1 surrogates
  # x: N x p features
  # threshold: an constant for significant correlation

  rank_cor <- abs(cor(s, x, method = "spearman"))
  ind <- which(rank_cor > threshold)
  return(ind)
}

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

  n1_subset <- floor(length(ind1) * p_subset)
  n0_subset <- floor(length(ind0) * p_subset)
  n_subset <- min(n1_subset, n0_subset, n_subset)

  # Weights for Adaptive lasso
  w <- rep(1, 2 * n_subset)
  p <- ncol(x)
  lambda_init <- p / (2 * n_subset)

  beta <- c()
  for (i in c(1:reps)) {
    ind <- c(sample(ind1, n_subset), sample(ind0, n_subset))
    beta_tmp <- fit_alasso_bic(s_label[ind], x[ind, ], family = "binomial", x_standardize = F)$beta_hat
    beta_tmp <- beta_tmp[-1] # drop intercept
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

clustering_method <- function(s, x, mclust_model = NULL) {
  s <- as.matrix(s)
  N <- nrow(s)
  tmpind <- 1:N

  # get pi_S via clustering or just use S if univariate/non-binary
  if (ncol(s) > 1) {
    mclustfit <- Mclust(s, G = 2, modelNames = mclust_model)
    pi_s <- ProbD.S(s, par = mclustfit$par)
    beta <- fit_alasso_approx(pi_s, x, family = "gaussian")
  } else {
    pi_s <- s[, 1]
    beta <- fit_alasso_bic(pi_s, x, family = "gaussian")$beta_hat
  }

  beta <- beta[-1]
  beta_select <- which(beta != 0)

  return(list(beta_select = beta_select, beta = beta, cluster_model = mclustfit))
}

# Computes probabilities after clustering
ProbD.S <- function(Si, par) {
  par.list <- list(
    pro = par$pro, mu = matrix(par$mean, ncol = 2),
    var = list(1, 1)
  )
  Si <- as.matrix(Si)
  k1 <- which.max(apply(par.list$mu, 2, mean))
  k0 <- setdiff(1:2, k1)
  if (ncol(Si) == 1) {
    sig2 <- par$variance$sigmasq
    if (length(sig2) == 1) {
      sig2 <- rep(sig2, 2)
    }
    par.list$var <- as.list(sig2)
  } else {
    for (kk in 1:2) {
      par.list$var[[kk]] <- par$variance$sigma[, , kk]
    }
  }
  tmp1 <- dmvnorm(Si, mean = par.list$mu[, k1], sigma = as.matrix(par.list$var[[k1]])) *
    par$pro[k1]
  tmp0 <- dmvnorm(Si, mean = par.list$mu[, k0], sigma = as.matrix(par.list$var[[k0]])) *
    par$pro[k0]
  tmp1 / (tmp1 + tmp0)
}

#################################################################
################### Model fitting methods #######################
#################################################################

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
  w_alasso <- init_weights(x, y,
    family = family, weights = weights, x_standardize = x_standardize,
    offset = offset, init_lambda = init_lambda, init_opt = init_opt
  )


  # Step 2. Choose the best parameter according to modified BIC
  fit_tmp <- glmnet(x, y,
    family = family, weight = weights, offset = offset,
    penalty.factor = w_alasso, alpha = 1,
    nlambda = 100L, lambda.min.ratio = 1e-3
  )

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

  # best_lambda <- fit_tmp$lambda[which.min(bics)]
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

  # ddl is the second derivative of the log likelihood.
  ddl <- solve(summary(init_fit)$cov.unscaled)
  ddl_half <- svd(ddl)
  Xtilde <- ddl_half$u %*% diag(sqrt(ddl_half$d)) %*% t(ddl_half$u)
  Ytilde <- Xtilde %*% init_coef

  # Step 2: Fit alasso. Choose the best parameter according to BIC.
  fit_tmp <- glmnet(
    x = Xtilde, y = Ytilde, penalty.factor = w_alasso, alpha = 1,
    nlambda = 100L, lambda.min.ratio = 1e-3
  )

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
    init_fit <- glmnet(
      x = x, y = y, family = family, weights = weights, offset = offset,
      alpha = 0, standardize = x_standardize, lambda = init_lambda
    )
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
  bic <- dev + df * log(n)
  return(bic)
}

# Modified BIC in (Minnier et al. 2011), which replace log(n) by min{n^0.1, log(n)}
# to avoid over shrinkage.
modified_bic <- function(dev, n, bic_factor, df) {
  mod_bic <- dev + min(n^bic_factor, log(n)) * df
  return(mod_bic)
}

# Vector to matrix.
vec2mat <- function(vc, dm) {
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}

# Compute the AUC.

cal_AUC <- function(data) {
  dd <- data[, 1]
  xx <- data[, 2]
  n0 <- sum(1 - dd)
  n1 <- sum(dd)
  x0 <- xx[dd == 0]
  x1 <- xx[dd == 1]
  sum((sum.I(x0, "<=", x1) + sum.I(x0, "<", x1)) / 2) / (n0 * n1)
}

# Computes sums efficiently based on ranks
sum.I <- function(yy, FUN, Yi, Vi = NULL, ties.method = "first") {
  if (FUN == "<" | FUN == ">=") {
    yy <- -yy
    Yi <- -Yi
  }
  pos <- rank(c(yy, Yi), ties.method = ties.method)[1:length(yy)] - rank(yy, ties.method = ties.method)
  if (substring(FUN, 2, 2) == "=") pos <- length(Yi) - pos
  if (!is.null(Vi)) {
    if (substring(FUN, 2, 2) == "=") tmpind <- order(-Yi) else tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind, , drop = F], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  } else {
    return(pos)
  }
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
  test_auc_ci <- ci(test_auc, method = "bootstrap", boot.n = 1000, boot.stratified = F)
  test_auc_ci <- round(test_auc_ci, 3)
  test_auc <- test_auc_ci[2]
  test_auc_ci_len <- as.numeric(test_auc_ci[3] - test_auc_ci[1])
  test_auc_ci <- paste0("(", test_auc_ci[1], ",", test_auc_ci[3], ")")

  return(list(
    train_AUC = train_auc,
    train_CI = train_auc_ci,
    train_CI_length = train_auc_ci_len,
    test_AUC = test_auc,
    test_CI = test_auc_ci,
    test_CI_length = test_auc_ci_len
  ))
}



## roc and auc ----
ROC.Est.FUN <- function(Di, yyi, yy0, fpr0 = NULL, wgti = NULL, yes.smooth = F) {
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- out.F1 <- NULL
  if (is.null(wgti)) {
    wgti <- rep(1, length(Di))
  }
  yyi <- as.matrix(yyi)
  pp <- ncol(as.matrix(yyi))
  mu0 <- sum(wgti * (1 - Di)) / sum(wgti)
  mu1 <- 1 - mu0
  for (k in 1:pp)
  {
    yy <- yy0
    if (!is.null(fpr0)) {
      tpr.all <- S.FUN(yyi[, k], Yi = yyi[, k], Di * wgti, yes.smooth = yes.smooth)
      fpr.all <- S.FUN(yyi[, k], Yi = yyi[, k], (1 - Di) * wgti, yes.smooth = yes.smooth)
      TPR <- approx(c(0, fpr.all, 1), c(0, tpr.all, 1), fpr0, method = "linear", rule = 2)$y
      TPR <- c(S.FUN(yy0, Yi = yyi[, k], Di * wgti, yes.smooth = yes.smooth), TPR)
      yy <- c(yy, Sinv.FUN(fpr0, Yi = yyi[, k], (1 - Di) * wgti, yes.smooth = yes.smooth))
      FPR <- S.FUN(yy, Yi = yyi[, k], (1 - Di) * wgti, yes.smooth = yes.smooth)
    } else {
      TPR <- S.FUN(yy, Yi = yyi[, k], Di * wgti, yes.smooth = yes.smooth)
      FPR <- S.FUN(yy, Yi = yyi[, k], (1 - Di) * wgti, yes.smooth = yes.smooth)
    }
    out.yy <- cbind(out.yy, yy)
    out.pp <- cbind(out.pp, S.FUN(yy, Yi = yyi[, k], wgti, yes.smooth = yes.smooth))
    out.TPR <- cbind(out.TPR, TPR)
    out.FPR <- cbind(out.FPR, FPR)
    PPV <- 1 / (1 + FPR * mu0 / (TPR * mu1))
    NPV <- 1 / (1 + (1 - TPR) * mu1 / ((1 - FPR) * mu0))
    out.PPV <- cbind(out.PPV, PPV)
    out.NPV <- cbind(out.NPV, NPV)
    AUC <- sum(S.FUN(yyi[, k], Yi = yyi[, k], Di * wgti, yes.smooth = yes.smooth) * (1 - Di) * wgti) / sum((1 - Di) * wgti)
    F1 <- 2 * TPR * PPV / (TPR + PPV)
    F1[is.na(F1) == 1] <- 0
    out.AUC <- c(out.AUC, AUC)
    out.F1 <- c(out.F1, F1)
  }
  out <- c(out.AUC, out.yy, out.pp, out.FPR, out.TPR, out.PPV, out.NPV, out.F1)
  out
}

get_roc <- function(y_true, y_score, subject_weight = NULL) {
  junk <- ROC.Est.FUN(y_true, y_score, yy0 = 0.5, fpr0 = seq(0, 1, 0.01), wgti = subject_weight, yes.smooth = F)[-1]
  df <- matrix(junk, ncol = 7)[-1, ]
  colnames(df) <- c("cutoff", "pos.rate", "FPR", "TPR", "PPV", "NPV", "F1")
  return(df)
}
get_auc <- function(y_true, y_score, subject_weight = NULL) {
  ROC.Est.FUN(y_true, y_score, yy0 = 0.5, fpr0 = seq(0, 1, 0.01), wgti = subject_weight, yes.smooth = F)[1]
}

S.FUN <- function(yy, Yi, Di, yes.smooth = F) {
  if (yes.smooth) {
    Y1i <- Yi[Di == 1]
    n1 <- sum(Di)
    bw <- bw.nrd(Y1i) / n1^0.6
    c(t(rep(1 / n1, n1)) %*% pnorm((Y1i - VTM(yy, n1)) / bw))
  } else {
    return((sum.I(yy, "<", Yi, Vi = Di) + sum.I(yy, "<=", Yi, Vi = Di)) / sum(Di) / 2)
  }
  ## sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

sum.I <- function(yy, FUN, Yi, Vi = NULL)
                  ## sum_i I(yy FUN Yi)Vi
# Vi weight
{
  if (FUN == "<" | FUN == ">=") {
    yy <- -yy
    Yi <- -Yi
  }
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy, Yi), ties.method = "f")[1:length(yy)] - rank(yy, ties.method = "f")
  if (substring(FUN, 2, 2) == "=") pos <- length(Yi) - pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if (substring(FUN, 2, 2) == "=") tmpind <- order(-Yi) else tmpind <- order(Yi)
    ## Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind, , drop = F], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  } else {
    return(pos)
  }
}
Sinv.FUN <- function(uu, Yi, Di, yes.smooth = F) {
  yy0 <- unique(sort(Yi, decreasing = T))
  ss0 <- S.FUN(yy0, Yi, Di, yes.smooth = yes.smooth)
  return(approx(ss0[!duplicated(ss0)], yy0[!duplicated(ss0)], uu, method = "linear", f = 0, rule = 2)$y)
}

VTM <- function(vc, dm) {
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}

balanced_cv_fold <- function(kfolds,
                             n = NULL, y = NULL, seed = NULL) {
  stopifnot(!is.null(n) || !is.null(y))
  stopifnot(kfolds >= 2)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(y)) {
    n <- length(y)
  } else {
    n <- as.integer(n)
    y <- runif(n)
  }
  kfolds <- as.integer(kfolds)

  quo <- n %/% kfolds # quotient
  rem <- n %% kfolds # remainder
  offset <- rep.int(
    seq_len(quo),
    c(rep_len(kfolds + 1L, rem), rep_len(kfolds, quo - rem))
  )
  u <- runif(n, 0, 1) + offset

  i <- sort.list(y) # get sorting indexes
  i <- i[sort.list(u)] # permute indexes only within similar y values
  fold <- integer(n)
  fold[i] <- rep_len(seq_len(kfolds), n)

  fold
}


get_lambda_sequence <- function(x, y, foldid = NULL,
                                nlambda = 50L, lambda_max_multiplier = 1.5,
                                lambda_min_ratio = 1 / length(y)) {
  if (is.null(foldid)) {
    lambda_max <- max(abs(crossprod(x, y - mean(y)) / length(y)))
  } else {
    lambda_max <- max(sapply(unique(foldid), function(f) {
      i <- which(foldid != f)
      x1 <- x[i, ]
      y1 <- y[i]
      max(abs(crossprod(x1, y1 - mean(y1)) / length(y1)))
    }))
  }
  lambda_max <- lambda_max * lambda_max_multiplier
  lambda_min <- pmax(lambda_max * lambda_min_ratio, 1e-6)
  lambda <- exp(seq(log(lambda_max), log(lambda_min),
    length.out = nlambda
  ))
  lambda
}


get_estimate <- function(z, y, family, l1_fraction, lambda, factor, exclude,
                         tuning, foldid, df_multiplier) {
  if (tuning == "cv") {
    model <- cv.glmnet(
      z, y,
      foldid = foldid,
      family = family, alpha = l1_fraction,
      standardize = FALSE, intercept = TRUE,
      lambda = lambda * mean(factor),
      penalty.factor = factor / mean(factor), exclude = exclude
    )
    beta <- as.double(coef(model, s = "lambda.min"))
    tuning <- data.frame(
      lambda = model$lambda, nonzero = model$nzero + 1L,
      cvm = model$cvm, cvsd = model$cvsd
    )
  } else if (tuning == "ic") {
    model <- glmnet(
      z, y,
      family = family, alpha = l1_fraction,
      standardize = FALSE, intercept = TRUE,
      lambda = lambda * mean(factor),
      penalty.factor = factor / mean(factor), exclude = exclude
    )
    nonzero <- sapply(predict(model, type = "nonzero"), length) + 1L
    ic <- deviance(model) + df_multiplier * nonzero
    beta <- as.double(coef(model, s = model$lambda[which.min(ic)]))
    tuning <- data.frame(
      lambda = model$lambda, nonzero = nonzero,
      ic = ic
    )
  } else {
    stop("Invalid value for 'tuning'")
  }
  list(beta = beta, tuning = tuning)
}


lasso_fit <- function(x, y, family = "gaussian", l1_fraction = 1.0,
                      factor = 1.0, exclude = integer(0L),
                      tuning = c("cv", "ic"), kfolds = 10L, foldid = 123456L,
                      df_multiplier = "log(n)",
                      ...) {
  tuning <- match.arg(tuning)
  if (tuning == "cv") {
    if (length(foldid) == 1L) {
      foldid <- balanced_cv_fold(kfolds, y = y, seed = foldid)
    } else if (length(foldid) != nrow(x)) {
      foldid <- rep_len(foldid, nrow(x))
    }
  } else {
    foldid <- NULL
    df_multiplier <- eval(
      parse(text = df_multiplier), list(n = length(y))
    )
  }

  if (length(factor) != ncol(x)) {
    factor <- rep_len(factor, ncol(x))
  }
  z <- scale(x)
  z[is.na(z)] <- 0.0
  lambda <- get_lambda_sequence(z, y, foldid)

  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude,
    tuning, foldid, df_multiplier
  )
  beta <- result$beta

  intercept <- beta[1L]
  slope <- beta[-1L]
  slope <- ifelse(
    attr(z, "scaled:scale") > 1e-8,
    slope / attr(z, "scaled:scale"), 0.0
  )
  intercept <- intercept - sum(attr(z, "scaled:center") * slope)
  beta <- c(intercept, slope)

  beta <- structure(beta, foldid = foldid, tuning = result$tuning)
  beta
}


adaptive_lasso_fit <- function(x, y, family = "gaussian", l1_fraction = 1.0, power = 1.0,
                               factor = 1.0, exclude = integer(0L),
                               tuning = c("cv", "ic"), kfolds = 10L, foldid = 123456L,
                               df_multiplier = "log(n)",
                               ...) {
  tuning <- match.arg(tuning)
  if (tuning == "cv") {
    if (length(foldid) == 1L) {
      foldid <- balanced_cv_fold(kfolds, y = y, seed = foldid)
    } else if (length(foldid) != nrow(x)) {
      foldid <- rep_len(foldid, nrow(x))
    }
  } else {
    foldid <- NULL
    df_multiplier <- eval(
      parse(text = df_multiplier), list(n = length(y))
    )
  }

  if (length(factor) != ncol(x)) {
    factor <- rep_len(factor, ncol(x))
  }
  z <- scale(x)
  z[is.na(z)] <- 0.0
  lambda <- get_lambda_sequence(z, y, foldid)

  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude,
    tuning, foldid, df_multiplier
  )
  beta <- result$beta

  # drop intercept
  if (sum(abs(beta[-1L]) > 1e-6) == 0L) {
    beta[-1L] <- 0.0
    beta <- structure(beta, foldid = foldid, tuning = result$tuning)
    return(beta)
  }

  factor <- 1.0 / abs(beta[-1L])**power
  factor <- factor / min(factor)
  threshold <- min(1e6, max(factor) * 0.99)
  exclude <- which(factor > threshold)
  factor[exclude] <- threshold

  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude,
    tuning, foldid, df_multiplier
  )
  beta <- result$beta

  intercept <- beta[1L]
  slope <- beta[-1L]
  slope <- ifelse(
    attr(z, "scaled:scale") > 1e-8,
    slope / attr(z, "scaled:scale"), 0.0
  )
  intercept <- intercept - sum(attr(z, "scaled:center") * slope)
  beta <- c(intercept, slope)

  beta <- structure(beta, foldid = foldid, tuning = result$tuning)
  return(beta)
}


linear_model_predict <- function(beta, x, probability = FALSE) {
  y <- drop(x %*% beta[-1L]) + beta[1L]
  if (probability) {
    y <- plogis(y)
  }
  y
}

validate_supervised <- function(dat, nsim, n.train = c(50, 70, 90)) {
  temp <- parallel::mclapply(1:nsim, FUN = function(n) {
    set.seed(12345+n)
    id.x <- lapply(n.train, function(n) sample(dat$patient_id, size = n))
    id.y <- lapply(id.x, function(i) {
      sample(dat$patient_id[which(!(dat$patient_id %in% i))])
    })

    lasso <- sapply(1:3, function(i) {
      auc_roc(
        actuals = ehr_data[id.y[[i]], ]$label,
        preds = linear_model_predict(
          beta =
            lasso_fit(
              x = ehr_data[id.x[[i]], 3:ncol(ehr_data)],
              y = ehr_data[id.x[[i]], ]$label,
              family = "binomial",
              tuning = "cv"
            ),
          x = as.matrix(ehr_data[id.y[[i]], 3:ncol(ehr_data)]),
          probability = TRUE
        )
      )
    })
    alasso <- sapply(1:3, function(i) {
      auc_roc(
        actuals = ehr_data[id.y[[i]], ]$label,
        preds = linear_model_predict(
          beta =
            adaptive_lasso_fit(
              x = ehr_data[id.x[[i]], 3:ncol(ehr_data)],
              y = ehr_data[id.x[[i]], ]$label,
              family = "binomial",
              tuning = "cv"
            ),
          x = as.matrix(ehr_data[id.y[[i]], 3:ncol(ehr_data)]),
          probability = TRUE
        )
      )
    })
    # ss <- sapply(1:3, function(i) {
    #   auc_roc(
    #     actuals = ehr_data[id.y[[i]], ]$label,
    #     preds = twostep_pred(
    #       train_data = ehr_data[id.x[[i]], ],
    #       test_data = ehr_data[id.y[[i]], ],
    #       X = ehr_data[, 6:ncol(ehr_data)],
    #       S = sicdnlp,
    #       health_count = health_count,
    #       beta.step1 = beta.step1
    #     )
    #   )
    # })

    #c(lasso, alasso, ss)
    c(lasso, alasso)
  })
  temp <- do.call(rbind.data.frame, temp)
  colnames(temp) <- outer(paste0("n=", c(50, 70, 90)), 
                          c("LASSO", "ALASSO"), paste, sep = ",")
  return(temp)
}

# 2 Step AUC
twostep_pred <- function(train_data, test_data, X, S, health_count, beta.step1) {

  # linear predictor without intercept
  bhatx <- linear_model_predict(beta = beta.step1, x = as.matrix(X))
  # Step 2
  idx <- train_data$patient_id
  idy <- test_data$patient_id
  model2ssl_step2 <- glm(
    train_data$label ~ bhatx[idx] + S[idx] + health_count[idx]
  )
  beta_step2 <- coef(model2ssl_step2)
  # recover beta
  beta <- beta_step2[2] * beta.step1
  mu <- beta_step2[1] +
    as.numeric(as.matrix(X[test_data$patient_id, ]) %*% beta[-1]) +
    as.numeric(beta_step2[3] %*% S[test_data$patient_id]) +
    as.numeric(beta_step2[4] %*% health_count[test_data$patient_id])
  # expit
  y_hat.ss <- plogis(mu)
}

validate_ss <- function(dat, nsim, n.train = c(50, 70, 90), bhatx) {
  temp <- parallel::mclapply(1:nsim, FUN = function(n) {
    set.seed(Sys.time())
    id.x <- lapply(n.train, function(n) sample(dat$patient_id, size = n))
    id.y <- lapply(id.x, function(i) {
      sample(dat$patient_id[which(!(dat$patient_id %in% i))])
    })
    
    lasso <- sapply(1:3, function(i) {
      auc_roc(
        actuals = ehr_data[id.y[[i]], ]$label,
        preds = linear_model_predict(
          beta =
            lasso_fit(
              x = ehr_data[id.x[[i]], 3:ncol(ehr_data)],
              y = ehr_data[id.x[[i]], ]$label,
              family = "binomial",
              tuning = "cv"
            ),
          x = as.matrix(ehr_data[id.y[[i]], 3:ncol(ehr_data)]),
          probability = TRUE
        )
      )
    })
    alasso <- sapply(1:3, function(i) {
      auc_roc(
        actuals = ehr_data[id.y[[i]], ]$label,
        preds = linear_model_predict(
          beta =
            adaptive_lasso_fit(
              x = ehr_data[id.x[[i]], 3:ncol(ehr_data)],
              y = ehr_data[id.x[[i]], ]$label,
              family = "binomial",
              tuning = "cv"
            ),
          x = as.matrix(ehr_data[id.y[[i]], 3:ncol(ehr_data)]),
          probability = TRUE
        )
      )
    })
    # ss <- sapply(1:3, function(i) {
    #   auc_roc(
    #     actuals = ehr_data[id.y[[i]], ]$label,
    #     preds = twostep_pred(
    #       train_data = ehr_data[id.x[[i]], ],
    #       test_data = ehr_data[id.y[[i]], ],
    #       X = ehr_data[, 6:ncol(ehr_data)],
    #       S = sicdnlp,
    #       health_count = health_count,
    #       beta.step1 = beta.step1
    #     )
    #   )
    # })
    
    #c(lasso, alasso, ss)
    c(lasso, alasso)
  })
  temp <- do.call(rbind.data.frame, temp)
  colnames(temp) <- outer(paste0("n=", c(50, 70, 90)), 
                          c("LASSO", "ALASSO"), paste, sep = ",")
  return(temp)
}