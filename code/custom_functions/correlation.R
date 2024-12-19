compute_lrv = function(data) {
  n = nrow(data)
  log_data = log(data)
  log_data[is.infinite(log_data)] = NA
  
  # CLR transformation by features
  shift = apply(log_data, 2, function(x) mean(x, na.rm = TRUE))
  clr_data = t(t(log_data) - shift)
  
  # lrv is proportional to the squared Euclidean distance in CLR coordinates.
  clr_dist = as.matrix(dist(t(clr_data)))
  lrv = clr_dist^2/(n-1)
  rownames(lrv) = colnames(data)
  colnames(lrv) = colnames(data)
  return(lrv)
}

compute_correlation = function(data) {
  n = nrow(data)
  d = ncol(data)
  
  log_data = log(data)
  log_data[is.infinite(log_data)] = NA
  
  # CLR transformation by samples
  shift = apply(log_data, 1, function(x) mean(x, na.rm = TRUE))
  clr_data = log_data - shift
  
  # Compute logratio variance (LRV)
  lrv = compute_lrv(data = data)
  
  # Calculate variance of log abundance for each feature
  clr_var = apply(clr_data, 2, function(x) var(x, na.rm = T))
  sum_log_var = sum(clr_var) * d / (d-1)
  log_var = (clr_var - 1/d^2 * sum_log_var) * d / (d-2) 
  
  # Calculate correlation between log abundances for each feature pair
  log_var_matrix1 = matrix(log_var, nrow = d, ncol = d, byrow = FALSE) 
  log_var_matrix2 = matrix(log_var, nrow = d, ncol = d, byrow = TRUE) 
  log_std_prod = sqrt(log_var_matrix1 * log_var_matrix2)
  rho = (lrv - log_var_matrix1 - log_var_matrix2)/ (-2 * log_std_prod)
  rownames(rho) = colnames(data)
  colnames(rho) = colnames(data)
  
  return(rho)
}

correlation_test = function(data, meta, num_sim, p_adj_method, shift=NULL) {
  n = nrow(data)
  d = ncol(data)
  sample_name = rownames(data)
  feature_name = colnames(data)
  
  if (is.null(meta)) {
    p = 1 #the intercept term
  } else {
    p = ncol(meta) + 1
  }
  
  # CLR transformation by samples
  log_data = log(data)
  log_data[is.infinite(log_data)] = NA
  if (is.null(shift)) {
    shift = apply(log_data, 1, function(x) mean(x, na.rm = TRUE))
  }
  clr_data = log_data - shift
  
  # Parameter estimation in CLR coordinates
  th = apply(data, 2, function(x) min(x[x != 0]))
  th = matrix(th, nrow = n, ncol = d, byrow = TRUE)
  clr_th = log(th) - shift
  
  init_params = apply(clr_data, 2, function(y) {
    if (is.null(meta)) {
      df = data.frame(y = y)
      lm_fit = lm(y ~ 1, data = df)
    } else {
      df = data.frame(y = y, meta)
      lm_fit = lm(y ~ ., data = df)
    }
    summ_fit = summary(lm_fit)
    estimates = c(coef(lm_fit), log_sigma = log(summ_fit$sigma))
    return(estimates)
  })
  init_params[is.na(init_params)] = 1
  
  set.seed(123)
  clr_params = lapply(seq_len(d), function(i) {
    censor_normal_mle(init_params = init_params[, i], meta = meta, 
                      data = clr_data[, i], th = clr_th[, i])
  })
  clr_params = Reduce('cbind', clr_params)
  clr_log_sd = clr_params[p + 1, ]
  clr_sd = exp(clr_log_sd)
  clr_mean = clr_params[1, ]
  
  # Deconfound data
  if (!is.null(meta)) {
    clr_coef = clr_params[2:p, , drop = FALSE]
    X = as.matrix(meta)
    clr_data = clr_data - X %*% clr_coef
  }
  
  # Zero percentages for each feature
  zero_perc = apply(data, 2, function(x) sum(x == 0)/n)
  
  # Hypothesis testing
  rho_list = vector("list", length = num_sim)
  for (i in seq_len(num_sim)) {
    set.seed(i)
    
    # Impute the missing values
    random_clr_data = matrix(rnorm(n * d,
                                   rep(clr_mean, each = n),
                                   rep(clr_sd, each = n)),
                             nrow = n, ncol = d)
    fill_clr_values = lapply(seq_len(d), function(i) {
      values = random_clr_data[, i]
      cutoff = quantile(random_clr_data[, i], zero_perc[i])
      return(values[values < cutoff])
    })
    impute_clr_data = lapply(seq_len(d), function(i) {
      orig_values = clr_data[, i]
      fill_values = fill_clr_values[[i]]
      impute_values = orig_values
      impute_values[is.na(impute_values)] = fill_values
      return(impute_values)
    })
    impute_clr_data = Reduce('cbind', impute_clr_data)
    impute_log_data = impute_clr_data + shift
    colnames(impute_log_data) = feature_name
    rownames(impute_log_data) = sample_name
    impute_data = exp(impute_log_data)
    
    # Calculate correlations using the current imputed data
    rho = compute_correlation(data = impute_data)
    rho_list[[i]] = rho
  }
  
  rho_mean = Reduce(function(x, y) {
    (replace(x, is.na(x), 0) + replace(y, is.na(y), 0))
  }, rho_list) / num_sim
  
  # Statistical inference by Fisher z-transformation
  z = 0.5 * log((1 + rho_mean) / (1 - rho_mean))
  se = 1 / sqrt(n - 3)
  z_score = z / se
  diag(z_score) = 0
  p_value = 2 * pmin(1 - pnorm(abs(z_score)), pnorm(abs(z_score)))
  diag(p_value) = 0
  q_value = matrix_p_adjust(p_value, method = p_adj_method)
  
  outputs = list(shift = shift, clr_th = clr_th, clr_params = clr_params, 
                 th = th, init_params = init_params,
                 estimate = rho_mean, p_value = p_value, q_value = q_value)
  return(outputs)
}
