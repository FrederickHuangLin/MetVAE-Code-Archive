pkg_list = c("doParallel", "doRNG", "Hmisc", "tidyverse")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

if (!require(BiocManager, quietly = TRUE)) {
  install.packages("BiocManager")
}

bioconductor_packages = c("S4Vectors", "TreeSummarizedExperiment", "ANCOMBC")
lapply(bioconductor_packages, function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(package)
  }
})

library(doParallel)
library(doRNG)
library(tidyverse)
library(S4Vectors)
library(TreeSummarizedExperiment)
library(ANCOMBC)

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

cl = makeCluster(10)
registerDoParallel(cl)

res_sim = foreach(i = simparams_list, .verbose = TRUE, .combine = rbind, .packages = c("MASS", "ANCOMBC", "S4Vectors", "TreeSummarizedExperiment"),
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
  
  # SECOM
  feature_table = t(data_miss)
  assays = SimpleList(counts = feature_table)
  smd = DataFrame(smd)
  
  tse = TreeSummarizedExperiment(assays = assays, colData = smd)
  res_linear = secom_linear(data = list(tse), assay_name = "counts",
                            pseudo = 0, prv_cut = 0.1, lib_cut = 0, corr_cut = 0.5, 
                            wins_quant = c(0.05, 0.95), method = "pearson", 
                            soft = FALSE, thresh_len = 20, n_cv = 10, 
                            thresh_hard = 0.3, n_cl = 1)
  
  # P-value filtering
  secom_values = res_linear$corr
  secom_p = res_linear$corr_p
  secom_q = matrix_p_adjust(secom_p, method = "BH")
  secom_cooccur = res_linear$mat_cooccur
  feature_names = colnames(secom_values)
  overlap = 10
  est_cor = secom_values
  est_cor = p_filter(est_cor, secom_q, max_p = 0.05)
  est_cor[secom_cooccur < overlap] = 0
  est_cor[is.na(est_cor)] = 0
  
  # Summary
  true_cor = true_cor[feature_names, feature_names]
  true_idx = (true_cor[lower.tri(true_cor)] != 0)
  est_idx = (est_cor[lower.tri(est_cor)] != 0)
  # TPR
  tpr1 = sum(est_idx * true_idx)/sum(true_idx)
  # FPR
  fpr1 = sum(est_idx * (!true_idx))/sum(!true_idx)
  # FDR
  fdr1 = sum(est_idx * (!true_idx))/sum(est_idx)
  if (is.nan(fdr1)) {fdr1 = 0}
  
  # Thresholding
  secom_values = res_linear$corr_th
  secom_cooccur = res_linear$mat_cooccur
  feature_names = colnames(secom_values)
  overlap = 10
  est_cor = secom_values
  est_cor[secom_cooccur < overlap] = 0
  est_cor[is.na(est_cor)] = 0
  
  # Summary
  true_cor = true_cor[feature_names, feature_names]
  true_idx = (true_cor[lower.tri(true_cor)] != 0)
  est_idx = (est_cor[lower.tri(est_cor)] != 0)
  # TPR
  tpr2 = sum(est_idx * true_idx)/sum(true_idx)
  # FPR
  fpr2 = sum(est_idx * (!true_idx))/sum(!true_idx)
  # FDR
  fdr2 = sum(est_idx * (!true_idx))/sum(est_idx)
  if (is.nan(fdr2)) {fdr2 = 0}
  
  c(tpr1, fpr1, fdr1, tpr2, fpr2, fdr2)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "sim_secom_confound.csv")