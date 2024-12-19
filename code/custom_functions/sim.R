library(MASS)

sim_data = function(n, d, cor_pairs = 0, mu, sigma = 1, x = NULL, cont_list = NULL, cat_list = NULL, da_prop = 0.1) {
  sample_names = paste0("s", seq_len(n))
  feature_names = paste0("f", seq_len(d))
  
  # Simulate the log intensities
  # Mean vector
  mu_vector = sample(mu, size = d, replace = TRUE)
  # Standard deviation vector
  sd_vector = rep(sigma, d)
  
  # Correlation matrix
  cor_matrix = diag(x = 1, nrow = d, ncol = d)
  
  # Check the maximum number of correlated pairs
  max_cor_pairs = (d * d - d) / 2
  if (cor_pairs > max_cor_pairs) {
    stop(paste0("The maximum number of correlated pairs is: ", max_cor_pairs, ". Please reduce the number of correlated pairs."))
  }
  
  if (cor_pairs != 0) {
    # Generate pairs of indices
    idx1 = seq.int(from = 1, by = 2, length.out = cor_pairs)
    idx2 = seq.int(from = 2, by = 2, length.out = cor_pairs)
    
    # Generate correlation values
    cor_pairs_matrix = cbind(idx1, idx2, 
                             sample(c(-0.7, -0.6, -0.5, 0.5, 0.6, 0.7), 
                                    size = cor_pairs, replace = TRUE))
    
    # Insert correlation values into the correlation matrix
    for (i in seq_len(cor_pairs)) {
      row_index = cor_pairs_matrix[i, 1]
      col_index = cor_pairs_matrix[i, 2]
      corr_value = cor_pairs_matrix[i, 3]
      
      cor_matrix[row_index, col_index] = corr_value
      cor_matrix[col_index, row_index] = corr_value
    }
  }
  colnames(cor_matrix) = feature_names
  rownames(cor_matrix) = feature_names
  
  # The covariance matrix
  cov_matrix = cor_matrix * (sd_vector %o% sd_vector)
  
  # Simulate log abundances
  log_y = mvrnorm(n = n, mu = mu_vector, Sigma = cov_matrix)
  
  # Clean the sample meta data
  if (!is.null(x)) {
    # Select continuous confounders
    x_cont = if (!is.null(cont_list)) x[, cont_list, drop = FALSE] else NULL
    # Select and process categorical confounders
    x_cat = if (!is.null(cat_list)) model.matrix(~ ., data = x[, cat_list, drop = FALSE])[, -1, drop = FALSE] else NULL
    # Combine continuous and categorical confounders
    x = if (!is.null(x_cont) & !is.null(x_cat)) {
      cbind(x_cont, x_cat)
    } else if (!is.null(x_cont)) {
      x_cont
    } else if (!is.null(x_cat)) {
      x_cat
    } else {
      stop('At least one of `cont_list` and `cat_list` should be not NULL.')
    }
  }
  
  # Simulate the effect sizes of confounders
  es = c(0, -2, -1, 1, 2)
  if (!is.null(x)) {
    p = ncol(x)
    # Generate effect sizes for the correlated pairs
    beta1 = matrix(sample(es, size = 2 * cor_pairs * p, replace = TRUE, prob = c(0, 0.25, 0.25, 0.25, 0.25)), nrow = 2 * cor_pairs)
    # Generate effect sizes for the uncorrelated pairs
    beta2 = matrix(sample(es, size = (d - 2 * cor_pairs) * p, replace = TRUE, prob = c(1 - da_prop, rep(da_prop / 4, 4))), nrow = d - 2 * cor_pairs)
    # Combine effect sizes
    beta = rbind(beta1, beta2)
    # Update log abundances
    log_y = log_y + as.matrix(x) %*% t(beta)
  } else {
    beta = NULL
  }
  
  # Transform back to the original scale
  y = exp(log_y)
  y = as.data.frame(y)
  rownames(y) = sample_names
  colnames(y) = feature_names
  
  # Outputs
  list(cor_matrix = cor_matrix, y = y, x = x, beta = beta)
}
