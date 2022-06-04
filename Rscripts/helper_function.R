library(glmnet)
library(glmpath)

#################################################################
################### Plotting  ###################################
#################################################################

plot_roc <- function(rocs, legend, n_total_method = 7, method_index) {
  
  label_text <- c()
  for (i in c(1:length(rocs))) {
    label_text <- c(label_text, paste0(legend[i], ", AUC=", round(auc(rocs[[i]]), 3)))
  }
  
  ggroc(rocs, legacy.axes = TRUE) +
    scale_colour_manual(values = scales::hue_pal()(n_total_method)[method_index], 
                        labels = label_text) +
    theme(legend.position = "bottom", text = element_text(size = 20)) +
    ggtitle("The operating receiver characteristic (ROC) curve") + 
    guides(color = guide_legend(title = element_blank(), ncol = 2)) + 
    labs(x = "False positive rate (FPR)", y = "True positive rate (TPR)")
}

plot_sims <- function(data, legend, n_total_method = 7, method_index) {
  
  plot_data <- tidyr::gather(data) %>%
    mutate(
      n = gsub(",.*$", "", key),
      method = sub(".*,\\s*", "", key)
    ) 

  plot_data %>%
    ggplot(aes(y = value, color = method)) +
    scale_colour_manual(values = scales::hue_pal()(n_total_method)[method_index],
                        labels = legend) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(. ~ n) +
    ggtitle("Area under the ROC curve (AUC) from 600 simulations") +
    theme(text = element_text(size = 20)) +
    theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks = element_blank())
}

#################################################################
################### ROC Calculation###### #######################
#################################################################
get_roc_parameter <- function(fpr, roc_curve){
  fpr0 <- roc_curve$FPR
  roc_curve %>% filter(FPR == fpr0[which.max(fpr0[fpr0 <= fpr])])
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

#################################################################
################### Model fitting methods #######################
#################################################################
balanced_cv_fold <- function(kfolds, n = NULL, y = NULL, seed = NULL) {
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
    set.seed(1234 + n)
    id.x <- lapply(n.train, function(n) sample(dat$patient_id, size = n))
    id.y <- lapply(id.x, function(i) {
      sample(dat$patient_id[which(!(dat$patient_id %in% i))], 46)
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
    c(lasso, alasso)
  })
  temp <- do.call(rbind.data.frame, temp)
  colnames(temp) <- outer(paste0("n=", c(50, 70, 90)),
    c("LASSO", "ALASSO"), paste,
    sep = ","
  )
  return(temp)
}

# 2 Step AUC
twostep_pred <- function(train_data, test_data, X, S, beta.step1) {

  # linear predictor without intercept
  bhatx <- linear_model_predict(beta = beta.step1, x = as.matrix(X))
  # Step 2
  idx <- train_data$patient_id
  idy <- test_data$patient_id
  model2ssl_step2 <- glm(
    train_data$label ~ bhatx[idx] + S[idx], family = "binomial"
  )
  beta_step2 <- coef(model2ssl_step2)
  # recover beta
  beta <- beta_step2[2] * beta.step1
  mu <- beta_step2[1] +
    as.numeric(as.matrix(X[test_data$patient_id, ]) %*% beta[-1]) +
    as.numeric(beta_step2[3] %*% S[test_data$patient_id])
  # expit
  y_hat.ss <- plogis(mu)
}

validate_phecap <- function(dat, surrogates, feature_selected, nsim, n.train = c(50, 70, 90)) {
  
  surrogate_matrix <- sapply(surrogates, function(surrogate) {
    rowSums(PheCAP::ehr_data[, surrogate$variable_names, drop = FALSE])
  })
  
  colnames(surrogate_matrix) <- sapply(surrogates, function(surrogate) {
    paste0(surrogate$variable_names, collapse = "&")
  })
  
  # Orthogonalize.
  label <- PheCAP::ehr_data$label
  healthcare <- PheCAP::ehr_data$healthcare_utilization
  other_features <- as.matrix(PheCAP::ehr_data[, setdiff(feature_selected$selected, c(colnames(surrogate_matrix))), drop = FALSE])
  other_features <- qr.resid(qr(cbind(1.0, surrogate_matrix, healthcare)), other_features)
  data_transformed <- data.frame(label = label, 
                          surrogate_matrix, 
                          healthcare_utilization = healthcare, 
                          other_features)
  
  temp <- parallel::mclapply(1:nsim, FUN = function(n) {
    set.seed(1234 + n)
    id.x <- lapply(n.train, function(n) sample(dat$patient_id, size = n))
    id.y <- lapply(id.x, function(i) {
      sample(dat$patient_id[which(!(dat$patient_id %in% i))], 46)
    })

    phecap <- sapply(1:3, function(i) {
      mltools::auc_roc(
        actuals = label[id.y[[i]]],
        preds = linear_model_predict(
            beta =
            lasso_fit(
              x = data_transformed[id.x[[i]], 2:ncol(data_transformed)],
              y = label[id.x[[i]]],
              family = "binomial",
              tuning = "cv"
            ),
          x = as.matrix(data_transformed[id.y[[i]], 2:ncol(data_transformed)]),
          probability = TRUE
        )
      )
    })
    phecap
  })

  temp <- do.call(rbind.data.frame, temp)

  colnames(temp) <- outer(paste0("n=", c(50, 70, 90)),
    c("PheCAP"), paste,
    sep = ","
  )
  return(temp)
}

validate_ss <- function(dat, nsim, n.train = c(50, 70, 90), beta, x, S) {
  temp <- parallel::mclapply(1:nsim, FUN = function(n) {
    set.seed(1234 + n)
    id.x <- lapply(n.train, function(n) sample(dat$patient_id, size = n))
    id.y <- lapply(id.x, function(i) {
      sample(dat$patient_id[which(!(dat$patient_id %in% i))], 46)
    })
    ss <- sapply(1:3, function(i) {
      auc_roc(
        actuals = ehr_data[id.y[[i]], ]$label,
        preds = twostep_pred(
          train_data = ehr_data[id.x[[i]], ],
          test_data = ehr_data[id.y[[i]], ],
          X = x,
          S = S,
          beta.step1 = beta
        )
      )
    })
    ss
  })
  temp <- do.call(rbind.data.frame, temp)
  colnames(temp) <- outer(paste0("n=", c(50, 70, 90)),
    c("Two-Step"), paste,
    sep = ","
  )
  return(temp)
}


validate_svmandrf <- function(dat, nsim, n.train = c(50, 70, 90)) {
  temp <- parallel::mclapply(1:nsim, FUN = function(n) {
    set.seed(1234 + n)
    id.x <- lapply(n.train, function(n) sample(dat$patient_id, size = n))
    id.y <- lapply(id.x, function(i) {
      sample(dat$patient_id[which(!(dat$patient_id %in% i))], 46)
    })

    rf <- sapply(1:3, function(i) {
      model <- rfsrc(y ~ .,
        data =
          data.frame(
            y = ehr_data[id.x[[i]], ]$label,
            x = ehr_data[id.x[[i]], 3:ncol(ehr_data)]
          )
      )
      auc_roc(
        actuals = ehr_data[id.y[[i]], ]$label,
        preds = predict(
          model,
          data.frame(x = ehr_data[id.y[[i]], 3:ncol(ehr_data)])
        )$predicted
      )
    })

    svm <- sapply(1:3, function(i) {
      model <- SVMMaj::svmmaj(
        y = ehr_data[id.x[[i]], ]$label,
        X = ehr_data[id.x[[i]], 3:ncol(ehr_data)]
      )
      auc_roc(
        actuals = ehr_data[id.y[[i]], ]$label,
        preds = predict(
          model,
          ehr_data[id.y[[i]], 3:ncol(ehr_data)]
        )
      )
    })

    c(rf, svm)
  })
  temp <- do.call(rbind.data.frame, temp)
  colnames(temp) <- outer(paste0("n=", c(50, 70, 90)),
    c("rf", "svm"), paste,
    sep = ","
  )
  return(temp)
}
