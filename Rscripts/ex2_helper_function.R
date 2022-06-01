source("../Rscripts/helper_function.R")

validate <- function(dat, nsim, n.train = c(50, 70, 90)){
  temp <- mclapply(1:n.sim, FUN = function(n) {
  set.seed(Sys.time())
  id.x <- lapply(n.train, function(n) sample(dat$patient_id, size = n))
  id.y <- lapply(id.x, function(i) 
    sample(labeled_data$patient_id[which(!(labeled_data$patient_id %in% i))]))
  
  lasso <- sapply(1:3, function(i){
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
  alasso <- sapply(1:3, function(i){
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
  
  c(lasso, alasso)
  })
  temp <- do.call(rbind.data.frame, temp)
  colnames(temp) <- outer(paste0("n=",c(50, 70, 90)), c("LASSO", "ALASSO"), paste, sep=",")
  return(temp)
  }

# Inverse logit function
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

# LASSO AUC
lasso_pred <- function(train_data, test_data) {
  model <- cv.glmnet(
    x = as.matrix(train_data[, 3:ncol(ehr_data)]),
    y = train_data$label,
    family = "binomial",
    alpha = 1 # default, LASSO
  )
  # prediction on testing set
  y_hat <- predict(model,
    newx = as.matrix(test_data[, 3:ncol(ehr_data)]),
    s = "lambda.min", type = "response"
  )
  return(as.numeric(y_hat))
}

lasso_auc <- function(train_data, test_data){
  return(auc_roc(test_data$label, lasso_pred(train_data, test_data)))
}

# ALASSO AUC
alasso_pred <- function(train_data, test_data) {
  model_alasso <- fit_alasso_bic(
    y = train_data$label,
    x = as.matrix(train_data[, 3:ncol(ehr_data)])
  )
  y_hat <- expit(cbind(1, as.matrix(test_data[, 3:ncol(ehr_data)]))
  %*% model_alasso$beta_hat) # Inverse Logit
  return(as.numeric(y_hat))
}

alasso_auc <- function(train_data, test_data){
  return(auc_roc(test_data$label, alasso_pred(train_data, test_data)))
}
# PheCAP AUC
phe_pred <- function(train_data, test_data) {
  model <- fit_alasso_bic(
    y = train_data$label,
    x = train_data[, c(other_feature, SAFE_feature)]
  )
  # prediction on testing set
  y_hat <- expit(cbind(
    1,
    as.matrix(test_data[, c(other_feature, SAFE_feature)])
  )
  %*% model$beta_hat) # Inverse Logit
  return(as.numeric(y_hat))
}

# 2 Step AUC
twostep_pred <- function(train_data, test_data, X, S) {
  # Step 1
  beta.step1 <- alasso_fit(
    y = S, # surrogate
    x = X, # all X
    family = "gaussian"
  )
  # linear predictor without intercept
  bhatx <- linear_model_predict(beta = beta.step1, x = X)
  # Step 2
  idx <- train_data$patient_id
  idy <- test_data$patient_id
  model2ssl_step2 <- glm(
    train_data$label ~ bhatx[idx] + S[idx]
  )
  beta_step2 <- coef(model2ssl_step2)
  # recover beta
  beta <- beta_step2[2] * beta.step1
  return(linear_model_predict(beta = beta., x = test_data[,3:ncol(test_data)],
                              probability = TRUE))
}



balanced_cv_fold <- function(kfolds, 
                             n = NULL, y = NULL, seed = NULL)
{
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
  
  quo <- n %/% kfolds  # quotient
  rem <- n %% kfolds   # remainder
  offset <- rep.int(seq_len(quo), 
                    c(rep_len(kfolds + 1L, rem), rep_len(kfolds, quo - rem)))
  u <- runif(n, 0, 1) + offset
  
  i <- sort.list(y)     # get sorting indexes
  i <- i[sort.list(u)]  # permute indexes only within similar y values
  fold <- integer(n)
  fold[i] <- rep_len(seq_len(kfolds), n)
  
  fold
}


get_lambda_sequence <- function(
  x, y, foldid = NULL,
  nlambda = 50L, lambda_max_multiplier = 1.5,
  lambda_min_ratio = 1 / length(y))
{
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
                    length.out = nlambda))
  lambda
}


get_estimate <- function(
  z, y, family, l1_fraction, lambda, factor, exclude, 
  tuning, foldid, df_multiplier)
{
  if (tuning == "cv") {
    model <- cv.glmnet(
      z, y, foldid = foldid,
      family = family, alpha = l1_fraction, 
      standardize = FALSE, intercept = TRUE,
      lambda = lambda * mean(factor),
      penalty.factor = factor / mean(factor), exclude = exclude)
    beta <- as.double(coef(model, s = "lambda.min"))
    tuning <- data.frame(
      lambda = model$lambda, nonzero = model$nzero + 1L,
      cvm = model$cvm, cvsd = model$cvsd)
  } else if (tuning == "ic") {
    model <- glmnet(
      z, y, 
      family = family, alpha = l1_fraction, 
      standardize = FALSE, intercept = TRUE,
      lambda = lambda * mean(factor),
      penalty.factor = factor / mean(factor), exclude = exclude)
    nonzero <- sapply(predict(model, type = "nonzero"), length) + 1L
    ic <- deviance(model) + df_multiplier * nonzero
    beta <- as.double(coef(model, s = model$lambda[which.min(ic)]))
    tuning <- data.frame(
      lambda = model$lambda, nonzero = nonzero,
      ic = ic)
  } else {
    stop("Invalid value for 'tuning'")
  }
  list(beta = beta, tuning = tuning)
}


lasso_fit <- function(
  x, y, family = "gaussian", l1_fraction = 1.0, 
  factor = 1.0, exclude = integer(0L), 
  tuning = c("cv", "ic"), kfolds = 10L, foldid = 123456L, 
  df_multiplier = "log(n)",
  ...)
{
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
      parse(text = df_multiplier), list(n = length(y)))
  }
  
  if (length(factor) != ncol(x)) {
    factor <- rep_len(factor, ncol(x))
  }
  z <- scale(x)
  z[is.na(z)] <- 0.0
  lambda <- get_lambda_sequence(z, y, foldid)
  
  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude, 
    tuning, foldid, df_multiplier)
  beta <- result$beta
  
  intercept <- beta[1L]
  slope <- beta[-1L]
  slope <- ifelse(
    attr(z, "scaled:scale") > 1e-8, 
    slope / attr(z, "scaled:scale"), 0.0)
  intercept <- intercept - sum(attr(z, "scaled:center") * slope)
  beta <- c(intercept, slope)
  
  beta <- structure(beta, foldid = foldid, tuning = result$tuning)
  beta
}


adaptive_lasso_fit <- function(
  x, y, family = "gaussian", l1_fraction = 1.0, power = 1.0, 
  factor = 1.0, exclude = integer(0L),
  tuning = c("cv", "ic"), kfolds = 10L, foldid = 123456L, 
  df_multiplier = "log(n)",
  ...)
{
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
      parse(text = df_multiplier), list(n = length(y)))
  }
  
  if (length(factor) != ncol(x)) {
    factor <- rep_len(factor, ncol(x))
  }
  z <- scale(x)
  z[is.na(z)] <- 0.0
  lambda <- get_lambda_sequence(z, y, foldid)
  
  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude, 
    tuning, foldid, df_multiplier)
  beta <- result$beta
  
  # drop intercept
  if (sum(abs(beta[-1L]) > 1e-6) == 0L) {
    beta[-1L] <- 0.0
    beta <- structure(beta, foldid = foldid, tuning = result$tuning)
    return(beta)
  }
  
  factor <- 1.0 / abs(beta[-1L]) ** power
  factor <- factor / min(factor)
  threshold <- min(1e6, max(factor) * 0.99)
  exclude <- which(factor > threshold)
  factor[exclude] <- threshold
  
  result <- get_estimate(
    z, y, family, l1_fraction, lambda, factor, exclude, 
    tuning, foldid, df_multiplier)
  beta <- result$beta
  
  intercept <- beta[1L]
  slope <- beta[-1L]
  slope <- ifelse(
    attr(z, "scaled:scale") > 1e-8, 
    slope / attr(z, "scaled:scale"), 0.0)
  intercept <- intercept - sum(attr(z, "scaled:center") * slope)
  beta <- c(intercept, slope)
  
  beta <- structure(beta, foldid = foldid, tuning = result$tuning)
  return(beta)
}


linear_model_predict <- function(beta, x, probability = FALSE)
{
  y <- drop(x %*% beta[-1L]) + beta[1L]
  if (probability) {
    y <- plogis(y)
  }
  y
}


assess_quality_real <- function(data, beta_true, beta_hat)
{
  link_true <- as.double(data$x %*% beta_true[-1L]) + beta_true[1L]
  response_true <- plogis(link_true)
  
  link_hat <- as.double(data$x %*% beta_hat[-1L]) + beta_hat[1L]
  response_hat <- plogis(link_hat)
  
  c(auc = get_auc(data$y, response_hat),
    excess_risk = 0.5 * (get_deviance(data$y, response_hat) - 
                           get_deviance(data$y, response_true)),
    mse_link = mean((link_hat - link_true) ** 2),
    mse_response = mean((data$y - response_hat) ** 2),
    bias_response = mean(response_hat - response_true),
    var_response = var(response_hat - response_true),
    var_ratio_response = var(response_hat) / var(response_true),
    l1_error = sum(abs(beta_hat - beta_true)),
    l2_error = sqrt(sum(beta_hat - beta_true) ** 2),
    l0_norm = sum(abs(beta_hat) > 1e-8),
    l1_norm = sum(abs(beta_hat)))
}


adjust_note_count <- function(u, note_count)
{
  u <- u - mean(u)
  note_count <- note_count - mean(note_count)
  gamma <- sum(u * note_count) / sum(note_count ** 2)
  u - gamma * note_count
}


clean_data <- function(dt, count_normalization)
{
  dt <- dt[!is.na(note_count), ]
  dt[, surrogate := icd_main + nlp_main]
  dt[, icd_main := NULL]
  dt[, nlp_main := NULL]
  
  xn <- grep("^(surrogate|note_count|C[0-9]{7})", names(dt), 
             ignore.case = TRUE, value = TRUE)
  for (j in xn) {
    set(dt, j = j, value = log1p(dt[[j]]))
  }
  
  nonzero <- dt[!is.na(label), sapply(.SD, function(z) sum(z > 0.5))]
  nonzero <- names(nonzero)[nonzero >= 1L]
  
  xn <- intersect(xn, nonzero)
  if (count_normalization) {
    for (j in setdiff(xn, "note_count")) {
      set(dt, j = j, value = adjust_note_count(dt[[j]], dt[["note_count"]]))
    }
  }
  
  y <- dt[, label]
  x <- as.matrix(
    dt[, c("surrogate", setdiff(xn, "surrogate")), with = FALSE])
  
  label_table <- table(y, exclude = NULL)
  surrogate_auc <- get_auc(y[!is.na(y)], x[!is.na(y), 1L])
  
  list(y = y, x = x,
       label_table = label_table, surrogate_auc = surrogate_auc)
}


fit_full_labeled_data <- function(
  x, y, 
  method = c("supervised_lasso", "supervised_adaptive_lasso", 
             "surrogate_only", "perfect_surrogate", 
             "naive_semisupervised", "adaptive_semisupervised",
             "prior_lasso_support", "prior_lasso_value"))
{
  method <- match.arg(method)
  compute_estimate <- get(method)
  
  beta <- compute_estimate(
    x, y, surrogate_index = 1L, kfolds = 10L, fold_seed = 1L)
  beta <- as.double(beta)
  
  beta
}




train_and_evaluate <- function(
  x, y, 
  method = c('U_lasso', "supervised_lasso", "supervised_adaptive_lasso", 
             "surrogate_only", "perfect_surrogate", 
             "naive_semisupervised", "adaptive_semisupervised",
             "prior_lasso_support", "prior_lasso_value", 'U_lasso_SS'),
  beta_good = NULL, num_replications = 5, num_bootstrap = 20, train_size = 50)
{
  method <- match.arg(method)
  compute_estimate <- get(method)
  
  labeled <- which(!is.na(y))
  unlabeled <- which(is.na(y))
  alpha <- double(ncol(x))
  if (method %in% c("perfect_surrogate", 
                    "naive_semisupervised", "adaptive_semisupervised",
                    "prior_lasso_support", "prior_lasso_value")) {
    alpha <- adaptive_lasso_fit(
      x[, -1L], x[, 1L],
      family = "gaussian", tuning = "ic",
      # kfolds = kfolds, foldid = fold_seed,
      df_multiplier = "log(n)")
  }
  if (method %in% c('U_lasso', 'U_lasso_SS')){
    alpha <- U_lasso_direction(x, y, 1L)
  }
  
  names(alpha) <- c("intercept", colnames(x)[-1L])
  
  auc_x_alpha <- get_auc(
    y[labeled],
    drop(x[labeled, -1L] %*% alpha[-1L]))
  auc_s <- get_auc(
    y[labeled],
    x[labeled, 1L])
  auc_list <- c(auc_x_alpha = auc_x_alpha, auc_s = auc_s)
  
  n <- length(labeled)
  
  if (is.null(beta_good)) {
    beta_good <- double(ncol(x) + 1L)
  }
  
  result <- lapply(seq_len(num_replications), function(k) {
    #print(k)
    set.seed(1234L + k)
    sample_id <- sample(labeled, length(labeled), replace = F)
    fold_size <- as.integer(0.25 * length(labeled))
    
    eva_mat <- c()
    for (fold in 1:4) {
      i_test <- sample_id[(1 + (fold - 1) * fold_size): (fold * fold_size)]
      i_train_all <- setdiff(sample_id, i_test)
      
      eva_mat_fold <- c()
      for (j in 1:num_bootstrap) {
        set.seed(66666 + j + 1000 * k)
        i_train <- sample(i_train_all, train_size, replace = F)
        i_train <- c(i_train, unlabeled)
        beta <- compute_estimate(
          x[i_train, ], y[i_train], surrogate_index = 1L,
          kfolds = 10L, fold_seed = 5678L + k, alpha = alpha)
        
        test_data <- list(x = x[i_test, ], y = y[i_test])
        evaluate <- list(quality = assess_quality_real(test_data, beta_good, beta), 
                         beta = beta)
        eva_mat_fold <- rbind(eva_mat_fold, evaluate$quality)
      }
      eva_mat_fold <- as.data.frame(eva_mat_fold)
      eva_mat <- rbind(eva_mat,
                       c(mean(eva_mat_fold$auc), sd(eva_mat_fold$auc),
                         mean(eva_mat_fold$mse_response), sd(eva_mat_fold$mse_response)))
    }
    colMeans(eva_mat)
  })
  
  metric <- rbindlist(lapply(result, function(r) as.data.table(rbind((r)))))
  metric <- colMeans(metric)
  print(method)
  return(metric)
}