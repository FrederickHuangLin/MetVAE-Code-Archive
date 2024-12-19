#==================================Compute VLR==================================
compute_vlr = function(data) {
  n = nrow(data)
  log_data = log(data)
  log_data[is.infinite(log_data)] = NA
  
  # CLR transformation by features
  shift = apply(log_data, 2, function(x) mean(x, na.rm = TRUE))
  clr_data = t(t(log_data) - shift)
  
  # vlr is proportional to the squared Euclidean distance in CLR coordinates.
  clr_dist = as.matrix(dist(t(clr_data)))
  vlr_mat = clr_dist^2/(n-1)
  rownames(vlr_mat) = colnames(data)
  colnames(vlr_mat) = colnames(data)
  return(vlr_mat)
}

norm_vlr = function(t) {
  d = nrow(t)
  t_colmed = apply(t, 2, median)
  t_med = median(t)
  
  colmed_matrix1 = matrix(t_colmed, nrow = d, ncol = d, byrow = FALSE)
  colmed_matrix2 = matrix(t_colmed, nrow = d, ncol = d, byrow = TRUE)
  
  if (any(colmed_matrix1 < 0)) {
    colmed_matrix1 = ifelse(colmed_matrix1 < 0, -colmed_matrix1, colmed_matrix1)
    colmed_matrix2 = ifelse(colmed_matrix2 < 0, -colmed_matrix2, colmed_matrix2)
  }
  
  norm_t = (t - colmed_matrix1 - colmed_matrix2 + t_med) / sqrt(colmed_matrix1 * colmed_matrix2)
  
  return(norm_t)
}

#============================VLR hypothesis testing=============================
by_sample_permute = function(data) {
  n = nrow(data)
  d = ncol(data)
  permuted_data = apply(data, 2, sample)
  permuted_data = matrix(permuted_data, nrow = n, ncol = d)
  return(permuted_data)
}

es_direction = function(p_below, p_above) {
  d = nrow(p_below)
  direct = matrix("increase", nrow = d, ncol = d) # increased co-occurrence
  direct[p_below > p_above] = "decrease" # decreased co-occurrence
  diag(direct) = "unchanged"
  return(direct)
}

vlr_test = function(data, meta, num_perm, p_adj_method) {
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
  shift = apply(log_data, 1, function(x) mean(x, na.rm = TRUE))
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
  vlr_list = vector("list", length = num_perm)
  vlr_below = vector("list", length = num_perm)
  vlr_above = vector("list", length = num_perm)
  for (i in seq_len(num_perm)) {
    set.seed(i)
    
    if (any(zero_perc != 0)) {
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
    } else {
      impute_data = data
    }
    
    # Calculate the VLR using the current imputed data
    vlr_values = compute_vlr(data = impute_data)
    vlr_values_norm = norm_vlr(vlr_values)
    
    # Calculate the permuted VLR using the current imputed data
    perm_data = by_sample_permute(data = impute_data)
    vlr_perm = compute_vlr(data = perm_data)
    vlr_perm_norm = norm_vlr(vlr_perm)
    
    vlr_list[[i]] = vlr_values_norm
    vlr_below[[i]] = (vlr_perm_norm <= vlr_values_norm)
    vlr_above[[i]] = (vlr_perm_norm >= vlr_values_norm)
  }
  
  vlr_mean = Reduce(function(x, y) {
    (replace(x, is.na(x), 0) + replace(y, is.na(y), 0))
  }, vlr_list) / num_perm
  
  vlr_p_below = Reduce(function(x, y) {
    (replace(x, is.na(x), 0) + replace(y, is.na(y), 0))
  }, vlr_below) / num_perm
  
  vlr_p_above = Reduce(function(x, y) {
    (replace(x, is.na(x), 0) + replace(y, is.na(y), 0))
  }, vlr_above) / num_perm
  
  vlr_p = 2 * pmin(vlr_p_below, vlr_p_above)
  diag(vlr_p) = 0
  vlr_q = matrix_p_adjust(vlr_p, method = p_adj_method)
  es_direct = es_direction(vlr_p_below, vlr_p_above)
  
  outputs = list(estimate = vlr_mean, direct = es_direct, 
                 p_value = vlr_p, q_value = vlr_q)
  return(outputs)
}
