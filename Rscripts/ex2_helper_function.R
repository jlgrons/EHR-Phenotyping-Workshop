source("../Rscripts/helper_function.R")

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

# 2Step AUC
twostep_pred <- function(train_data, test_data) {
  # Step 1
  model2ssl_step1 <- fit_alasso_bic(
    y = sicdnlp, # surrogate
    x = x, # all X
    family = "gaussian"
  )
  # linear predictor without intercept
  bhatx <- x %*% model2ssl_step1$beta_hat[-1]
  # Step 2
  idx <- train_data$patient_id
  idy <- test_data$patient_id
  model2ssl_step2 <- glm(
    train_data$label ~ bhatx[idx] + sicdnlp[idx] + health_count[idx]
  )
  beta_step2 <- coef(model2ssl_step2)
  # recover beta
  beta <- beta_step2[2] * model2ssl_step1$beta_hat[-1]
  mu <- beta_step2[1] +
    as.numeric(x[idy, ] %*% beta) +
    as.numeric(beta_step2[3] %*% sicdnlp[idy]) + as.numeric(beta_step2[4]
    %*% health_count[idy])
  y_hat <- expit(mu)
  return(y_hat)
}


# K fold
k_fold <- function(dat, k = 5, seed = 129) {
  set.seed(seed)
  idx <- labeled_data$patient_id %>% caret::createFolds(k)

  auc <- sapply(idx, function(i) {
    train_data <- dat[-i, ]
    test_data <- dat[i, ]
    c(
      auc_roc(lasso_pred(train_data, test_data), test_data$label),
      auc_roc(alasso_pred(train_data, test_data), test_data$label),
      auc_roc(phe_pred(train_data, test_data), test_data$label),
      auc_roc(twostep_pred(train_data, test_data), test_data$label)
    )
  })
  out <- apply(auc, 1, mean)
  names(out) <- c("LASSO", "ALASSO", "PheCAP", "Two Step SS")
  return(out)
}