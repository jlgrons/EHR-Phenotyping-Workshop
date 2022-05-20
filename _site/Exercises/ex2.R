library(PheCAP)
library(tidyverse)
library(glmnet)
library(pROC)

source("../Rscripts/rankCor.R")
source("../Rscripts/extremeMethod.R")
source("../Rscripts/clusteringMethod.R")
source("../Rscripts/modelFitting.R")

#-----------------------------------Load data------------------------------------------------#

data(ehr_data)
data <- PhecapData(ehr_data, "healthcare_utilization", "label", 0.4, patient_id = "patient_id")
data

#-----------------------------------Feature selection-----------------------------------------#

# Prepare data for feature selection, transformed. 
sicd <- log(ehr_data$main_ICD + 1)
snlp <- log(ehr_data$main_NLP + 1)
x <- data.matrix(ehr_data %>% select(starts_with("COD") | starts_with("NLP")))
x <- log(x + 1)

# Rank correlation.
AFEP_select <- rankCor(snlp, x, threshold = 0.15)
AFEP_feature <- colnames(x)[AFEP_select]
AFEP_feature

# Tail method.
SAFE_icd <- extreme_method(sicd, x)
SAFE_nlp <- extreme_method(snlp, x)
SAFE_both <- extreme_method(cbind(sicd, snlp), x)
beta <- rbind(SAFE_icd$beta_all, SAFE_nlp$beta_all, SAFE_both$beta_all)
SAFE_select <- which(colMeans(beta, na.rm = T) >= 0.5)
SAFE_feature <- colnames(x)[SAFE_select]
SAFE_feature

# Cluster method.
system.time(Auto <- clustering_method(cbind(sicd, snlp), x))
Auto_select <- Auto$beta_select
Auto_feature <- colnames(x)[Auto_select]
Auto_feature

#-----------------------------------Model development-------------------------------------------#
train_data <- ehr_data %>% filter(patient_id %in% data$training_set)
train_x <- train_data %>% select(starts_with("COD"), starts_with("NLP"), 
                                 starts_with("main"), healthcare_utilization)
train_y <- train_data %>% select(label) %>% pull()

train_x <- log(train_x + 1)

test_data <- ehr_data %>% filter(patient_id %in% data$validation_set) 
test_y <- test_data %>% select(label) %>% pull()

# Supervised - LASSO.
model_sl <- fit_alasso_bic(train_y, train_x, x_standardize = T)
modelselect_sl <- fit_alasso_bic(train_y, train_x[, Auto_feature], x_standardize = T)

# Supervised - Random forest. 
library(randomForestSRC)
model_rf <- rfsrc(y ~., data = data.frame(y = train_y, x = train_x))
#modelselect_rf <- rfsrc(y ~., data = data.frame(y = train_y, x = train_x[, Auto_feature]))

# Semi-supervised. Two-step.
model_twostep <- fit_alasso_bic(train_y, train_x[, Auto_feature], x_standardize = T)

# Semi-supervised. PASS. 


# Weakly-supervised. PheNorm.
library(PheNorm)
data_fit <- data.frame(cbind(train_x, train_y))
model_phenorm <- PheNorm.Prob(
  nm.logS.ori = "main_ICD", 
  nm.utl = "healthcare_utilization", 
  dat = data_fit, 
  corrupt.rate = 0.3, 
  train.size = nrow(data_fit)
)

# Weakly-supervised. MAP. 
library(MAP)
train_x <- train_data %>% select(starts_with("COD"), starts_with("NLP"), 
                                 starts_with("main"), healthcare_utilization)
data_fit <- Matrix(cbind(ICD = train_x$main_ICD, NLP = train_x$main_NLP), sparse = TRUE)
note <- Matrix(train_x$healthcare_utilization, ncol = 1, sparse = TRUE)
model_map <- MAP(mat = data_fit, note = note)











