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


