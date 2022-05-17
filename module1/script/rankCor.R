# Absolute rank correlation method.
# From the paper below:
# Sheng Yu, Katherine P Liao, Stanley Y Shaw, Vivian S Gainer, 
# Susanne E Churchill, Peter Szolovits, Shawn N Murphy, 
# Isaac S. Kohane, Tianxi Cai, 
# Toward high-throughput phenotyping: unbiased automated feature 
# extraction and selection from knowledge sources, Journal of the 
# American Medical Informatics Association, Volume 22, Issue 5, 
# September 2015, Pages 993–1000.

rankCor <- function (s, x, threshold = 0.15) {
  
  # s: N x 1 surrogates
  # x: N x p features
  # threshold: an constant for significant correlation 
  
  rank_cor <- abs(cor(s, x,  method = "pearson"))
  ind <- which(rank_cor > threshold)
  return(ind)
}

