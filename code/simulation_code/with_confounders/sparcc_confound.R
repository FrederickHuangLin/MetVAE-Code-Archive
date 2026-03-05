pkg_list = c("doParallel", "doRNG", "Hmisc", "tidyverse")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

if (!require(SpiecEasi, quietly = TRUE)) {
  if (!require(remotes, quietly = TRUE)) {
    install.packages("remotes")
    library(remotes)
  }
  remotes::install_github("zdk123/SpiecEasi")
}

library(doParallel)
library(doRNG)
library(tidyverse)
library(SpiecEasi)

folder_path = "../../functions"
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

cl = makeCluster(64)
registerDoParallel(cl)

res_sim = foreach(i = simparams_list, .verbose = TRUE, .combine = rbind, .packages = c("MASS", "SpiecEasi"),
                  .export = c("sim_data", "matrix_p_adjust", "p_filter")) %dorng% {
  params = strsplit(i, "_")[[1]]
  n = as.numeric(params[1])
  d = as.numeric(params[2])
  zero_prop = as.numeric(params[3])
  seed = as.numeric(params[4])
  cor_pairs = 0.2 * d
  mu = 10:14
  da_prop = 0.1
  
  set.seed(seed)
  smd = data.frame(x1 = rnorm(n = n),
                   x2 = sample(c("a", "b"), size = n, replace = TRUE))
  
  sim  =  sim_data(n = n, d = d, cor_pairs = cor_pairs, mu = mu, 
                   x = smd, cont_list = c("x1"), cat_list = c("x2"), da_prop = da_prop)
  
  y = sim$y
  x = sim$x
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
  
  # SparCC
  sparcc_values = sparcc(data_miss)$Cor
  colnames(sparcc_values) = colnames(data_miss)
  rownames(sparcc_values) = colnames(data_miss)
  
  boot_res = sparccboot(data_miss, R = 1000)
  p_res = pval.sparccboot(boot_res)
  sparcc_p = diag(0, nrow = d, ncol = d)
  sparcc_p[upper.tri(sparcc_p, diag = FALSE)] = p_res$pvals
  sparcc_p = sparcc_p + t(sparcc_p)
  colnames(sparcc_p) = colnames(data_miss)
  rownames(sparcc_p) = colnames(data_miss)
  sparcc_q = matrix_p_adjust(sparcc_p)
  est_cor = p_filter(sparcc_values, sparcc_q, max_p = 0.05)
  
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

write_csv(data.frame(res_sim), "sim_sparcc_confound.csv")