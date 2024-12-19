#=====================P-value filtering and adjustment==========================

p_filter = function(mat, mat_p, max_p, impute_value = 0) {
  mat_filter = mat
  mat_filter[mat_p > max_p] = impute_value
  return(mat_filter)
}

matrix_p_adjust = function(p_matrix, method = "BH") {
  p_vector = p_matrix[lower.tri(p_matrix)]
  q_vector = p.adjust(p_vector, method = method)
  q_matrix = diag(0, nrow = nrow(p_matrix))
  q_matrix[lower.tri(q_matrix)] = q_vector
  q_matrix[upper.tri(q_matrix)] = t(q_matrix)[upper.tri(q_matrix)]
  return(q_matrix)
}

#===========Estimate parameters for the censored Normal distribution============
log_likelihood = function(params, meta, data, th) {
  
  b = matrix(params[1:(length(params) - 1)], ncol = 1)
  sigma = exp(params[length(params)])
  
  if (is.null(meta)) {
    X = matrix(1, nrow = length(data), ncol = 1)
  } else {
    X = as.matrix(cbind(1, meta))
  }
  
  mu = X %*% b
  e = data - mu
  
  uncensored_idx = which(!is.na(data))
  censored_idx = which(is.na(data))
  
  ll_uncensored = - log(sigma) * length(uncensored_idx)-
    1 / (2 * sigma^2) * t(e[uncensored_idx, ]) %*% e[uncensored_idx, ]
  
  ll_censored = sum(pnorm((th[censored_idx] - mu[censored_idx, ])/sigma,
                          lower.tail = TRUE, log.p = TRUE))
  ll = as.numeric(ll_uncensored) + as.numeric(ll_censored)
  return(-ll)
}

gradient = function(params, meta, data, th) {
  
  b = matrix(params[1:(length(params) - 1)], ncol = 1)
  sigma = exp(params[length(params)])
  
  if (is.null(meta)) {
    X = matrix(1, nrow = length(data), ncol = 1)
  } else {
    X = as.matrix(cbind(1, meta))
  }
  
  mu = X %*% b
  e = data - mu
  
  uncensored_idx = which(!is.na(data))
  censored_idx = which(is.na(data))
  
  grad_b = -1 /(sigma^2) * t(X[uncensored_idx, ]) %*% e[uncensored_idx, ] +
    1/sigma * t(X[censored_idx, ]) %*% (dnorm((th[censored_idx] - mu[censored_idx, ])/sigma)/
                                          pnorm((th[censored_idx] - mu[censored_idx, ])/sigma, lower.tail = TRUE))
  grad_log_sigma = length(uncensored_idx) -
    (1/sigma^2) * t(e[uncensored_idx, ]) %*% e[uncensored_idx, ] +
    (1/sigma) * t(dnorm((th[censored_idx] - mu[censored_idx, ])/sigma)/
                    pnorm((th[censored_idx] - mu[censored_idx, ])/sigma, lower.tail = TRUE)) %*%
    (th[censored_idx] - mu[censored_idx, ])
  return(c(as.vector(grad_b), as.numeric(grad_log_sigma)))
}

censor_normal_mle = function(init_params, meta, data, th, fn = log_likelihood,
                             gr = gradient, method = "BFGS") {
  suppressWarnings(optim_results <- try(optim(par = init_params, fn = fn, gr = gr,
                                              meta = meta, data = data, th = th,
                                              method = method),
                                        silent = TRUE))
  if (inherits(optim_results, "try-error")) {
    outputs = init_params
  } else {
    outputs = optim_results$par
  }
  return(outputs)
}