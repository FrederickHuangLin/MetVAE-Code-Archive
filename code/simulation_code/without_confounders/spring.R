pkg_list = c("doParallel", "doRNG", "Hmisc", "tidyverse")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

if (!require(SPRING, quietly = TRUE)) {
  if (!require(devtools, quietly = TRUE)) {
    install.packages("devtools")
    library(devtools)
  }
  devtools::install_github("GraceYoon/SPRING")
}

library(doParallel)
library(doRNG)
library(tidyverse)
library(SPRING)

folder_path = "code/functions/"
r_files = list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  source(file)
}

n = 100
d = c(50, 200, 500)
zero_prop = seq(0, 0.3, 0.1)
iter_num = 100
seed = seq_len(iter_num)

simparams = data.frame(expand.grid(n, d, zero_prop, seed))
colnames(simparams) = c("n", "d", "zero_prop", "seed")
simparams = simparams %>%
  arrange(n, d, zero_prop, seed)
simparams_list = apply(simparams, 1, paste0, collapse = "_")

cl = makeCluster(10)
registerDoParallel(cl)

res_sim = foreach(i = simparams_list, .verbose = TRUE, .combine = rbind, .packages = c("MASS", "SpiecEasi", "SPRING"),
                  .export = c("sim_data")) %dorng% {
  params = strsplit(i, "_")[[1]]
  n = as.numeric(params[1])
  d = as.numeric(params[2])
  zero_prop = as.numeric(params[3])
  seed = as.numeric(params[4])
  cor_pairs = 0.2 * d
  mu = 10:14
  da_prop = 0.1
  
  set.seed(seed)
  sim  =  sim_data(n = n, d = d, cor_pairs = cor_pairs, mu = mu, da_prop = da_prop)
  y = sim$y
  true_cor = sim$cor_matrix
  
  # The log data
  log_y = log(as.matrix(y))
  
  # Add sources of biases
  log_sample_bias = log(runif(n, 1e-3, 1e-1))
  log_feature_bias = log(runif(d, 1e-1, 1))
  log_data = sweep(log_y, 1, log_sample_bias, `+`)
  log_data = sweep(log_data, 2, log_feature_bias, `+`)
  data = exp(log_data)
  
  # Add zeros
  thresholds = apply(data, 2, function(x) quantile(x, zero_prop))
  thresholds = matrix(thresholds, nrow = n, ncol = d, byrow = TRUE)
  data_miss = data
  data_miss[data_miss < thresholds] = 0
  
  # SPRING
  spring_mb = SPRING(data = data_miss, Rmethod = "approx", quantitative = TRUE, 
                     lambdaseq = "data-specific", nlambda = 50, rep.num = 50)
  opt.K = spring_mb$output$stars$opt.index
  est_cor = as.matrix(SpiecEasi::symBeta(spring_mb$output$est$beta[[opt.K]], mode = 'maxabs'))
  rownames(est_cor) = colnames(data_miss)
  colnames(est_cor) = colnames(data_miss)
  est_cor[is.na(est_cor)] = 0
  
  # Summary
  true_idx = (true_cor[lower.tri(true_cor)] != 0)
  est_idx = (est_cor[lower.tri(est_cor)] != 0)
  # TPR
  tpr = sum(est_idx * true_idx)/sum(true_idx)
  # FPR
  fpr = sum(est_idx * (!true_idx))/sum(!true_idx)
  # FDR
  fdr = sum(est_idx * (!true_idx))/sum(est_idx)
  
  c(tpr, fpr, fdr)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "sim_spring_unconfound.csv")