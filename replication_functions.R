library(dplyr)
library(tidyr)
library(scales)
library(glmnet)

# Two-sided p-value: 2 * Phi(-|estimate| / bootstrap SE).
pvalue_boot_se <- function(boot_draws, point_est) {
  x <- boot_draws[is.finite(boot_draws)]
  if (!length(x) || is.na(point_est)) return(NA_real_)
  se <- stats::sd(x)
  if (!is.finite(se)) return(NA_real_)
  if (se <= 0) {
    return(if (abs(point_est) < sqrt(.Machine$double.eps)) 1 else 0)
  }
  2 * stats::pnorm(-abs(point_est / se))
}

# Counterfactual rows: assign named columns by position (same pattern as simulation imp).
set_cols <- function(data, named_values) {
  out <- data
  for (nm in names(named_values)) out[[nm]] <- named_values[[nm]]
  out
}

#' Core "Proof of Concept" Heterogeneity Analysis
#'
#' Called by \code{bootstrap_poc()}.
poc <- function(df, Y, R, X, Zvars, Wvars, Wprime, w_fix) {
  safe_coef <- function(model_fit, term_name) {
    v <- coef(model_fit)[term_name]
    if (is.na(v) || !length(v)) 0 else unname(v)
  }
  rx <- paste(R, X, sep = ":")

  tau <- safe_coef(lm(as.formula(paste0(Y, " ~ ", R, " * ", X)), data = df), rx)

  fit_cond <- lm(
    as.formula(paste0(Y, " ~ ", R, " * ", X, " * (", paste(Zvars, collapse = " + "), ")")),
    data = df
  )
  zeta_z <- safe_coef(fit_cond, rx)
  for (zv in Zvars) {
    z_term <- paste(R, X, zv, sep = ":")
    zeta_z <- zeta_z + safe_coef(fit_cond, z_term) * mean(df[[zv]], na.rm = TRUE)
  }
  delta_z <- tau - zeta_z

  fit_out <- lm(
    as.formula(
      paste0(
        Y, " ~ ", R, " * ", X, " * (", paste(Zvars, collapse = " + "), ") + ",
        paste(Wvars, collapse = " + "), " + ", Wprime, " + ",
        X, ":", Wprime, " + ", R, ":", X, ":", Wprime
      )
    ),
    data = df
  )
  gamma_rx <- safe_coef(fit_out, rx)
  gamma_rxwprime <- safe_coef(fit_out, paste(R, X, Wprime, sep = ":"))
  zeta_static_poc <- gamma_rx + gamma_rxwprime * w_fix
  for (zv in Zvars) {
    z_term <- paste(R, X, zv, sep = ":")
    zeta_static_poc <- zeta_static_poc + safe_coef(fit_out, z_term) * mean(df[[zv]], na.rm = TRUE)
  }
  delta_static_poc <- zeta_z - zeta_static_poc

  fit_wp <- lm(
    as.formula(paste0(Wprime, " ~ ", paste(Zvars, collapse = " + "))),
    data = df[df[[R]] == 0, , drop = FALSE]
  )
  Ewprime_r0_z <- mean(as.numeric(predict(fit_wp, newdata = df)), na.rm = TRUE)
  zeta_stoch_poc <- gamma_rx + gamma_rxwprime * Ewprime_r0_z
  for (zv in Zvars) {
    z_term <- paste(R, X, zv, sep = ":")
    zeta_stoch_poc <- zeta_stoch_poc + safe_coef(fit_out, z_term) * mean(df[[zv]], na.rm = TRUE)
  }
  delta_stoch_poc <- zeta_z - zeta_stoch_poc

  list(
    zeta_z = zeta_z, delta_z = delta_z,
    zeta_static_poc = zeta_static_poc, delta_static_poc = delta_static_poc,
    zeta_stoch_poc = zeta_stoch_poc, delta_stoch_poc = delta_stoch_poc
  )
}

#' @param B Number of bootstrap replications.
#' @param ... Passed to \code{poc()}.
#' @return Tibble: estimand, estimate (mean of \code{B} bootstrap draws), std_error, p_value.
bootstrap_poc <- function(B = 1000, ...) {
  args <- list(...)
  df <- args$df
  if (is.null(df)) stop("A data frame must be provided with argument name 'df'.", call. = FALSE)

  point_estimates <- unlist(do.call(poc, args))
  est_names <- names(point_estimates)

  boot_results <- replicate(B, {
    idx <- sample.int(nrow(df), nrow(df), replace = TRUE)
    boot_args <- args
    boot_args$df <- df[idx, , drop = FALSE]
    tryCatch({
      vals <- unlist(do.call(poc, boot_args))
      vals[est_names]
    }, error = function(e) {
      setNames(rep(NA_real_, length(est_names)), est_names)
    })
  }, simplify = TRUE)

  boot_estimates_df <- as_tibble(t(boot_results))
  names(boot_estimates_df) <- est_names

  summary_stats <- tibble(
    estimand = est_names,
    estimate = vapply(est_names, function(nm) mean(boot_estimates_df[[nm]], na.rm = TRUE), numeric(1)),
    std_error = vapply(est_names, function(nm) sd(boot_estimates_df[[nm]], na.rm = TRUE), numeric(1)),
    p_value = vapply(est_names, function(nm) {
      m <- mean(boot_estimates_df[[nm]], na.rm = TRUE)
      pvalue_boot_se(boot_estimates_df[[nm]], m)
    }, numeric(1))
  )

  summary_stats %>%
    mutate(p_value = pvalue(p_value, accuracy = 0.001, add_p = FALSE)) %>%
    select(estimand, estimate, std_error, p_value)
}

#' Imputation-based heterogeneity (full linear models; simulation-aligned).
#'
#' @return Named list of six estimands (\code{zeta_z_imp}, \code{delta_z_imp}, static/stoch variants).
imputation <- function(df, Y, R, X, Zvars, Wvars, Wprime, w_fix = 0) {
  rx <- paste(R, X, sep = ":")
  tau <- unname(coef(lm(as.formula(paste0(Y, " ~ ", R, " * ", X)), data = df))[rx])

  form_cond <- as.formula(paste0(Y, " ~ ", R, " * ", X, " * (", paste(Zvars, collapse = " + "), ")"))
  fit_cond <- lm(form_cond, data = df)
  new_11 <- set_cols(df, setNames(list(1, 1), c(R, X)))
  new_10 <- set_cols(df, setNames(list(1, 0), c(R, X)))
  new_01 <- set_cols(df, setNames(list(0, 1), c(R, X)))
  new_00 <- set_cols(df, setNames(list(0, 0), c(R, X)))
  mu11 <- predict(fit_cond, newdata = new_11)
  mu10 <- predict(fit_cond, newdata = new_10)
  mu01 <- predict(fit_cond, newdata = new_01)
  mu00 <- predict(fit_cond, newdata = new_00)
  zeta_z_imp <- mean((mu11 - mu10) - (mu01 - mu00))
  delta_z_imp <- tau - zeta_z_imp

  Zc_vars <- paste0(Zvars, "_c")
  for (i in seq_along(Zvars)) df[[Zc_vars[i]]] <- df[[Zvars[i]]] - mean(df[[Zvars[i]]], na.rm = TRUE)

  form_out <- as.formula(
    paste0(
      Y, " ~ ", R, " * ", X, " * (", paste(Zc_vars, collapse = " + "), ") + ",
      paste(Wvars, collapse = " + "), " + ", Wprime, " + ",
      X, ":", Wprime, " + ", R, ":", X, ":", Wprime
    )
  )
  fit_out <- lm(form_out, data = df)
  new_x1 <- set_cols(df, setNames(list(1, w_fix), c(X, Wprime)))
  new_x0 <- set_cols(df, setNames(list(0, w_fix), c(X, Wprime)))
  mu_x1 <- predict(fit_out, newdata = new_x1)
  mu_x0 <- predict(fit_out, newdata = new_x0)
  mu_w <- mu_x1 - mu_x0
  dat_mu <- cbind(df, mu_w = mu_w)
  form_mu <- as.formula(paste0("mu_w ~ ", R, " + ", paste(Zc_vars, collapse = " + ")))
  zeta_static_imp <- unname(coef(lm(form_mu, data = dat_mu))[R])
  delta_static_imp <- zeta_z_imp - zeta_static_imp

  fit_wp <- lm(as.formula(paste0(Wprime, " ~ ", paste(Zvars, collapse = " + "))), data = df[df[[R]] == 0, , drop = FALSE])
  df_r1 <- df[df[[R]] == 1, , drop = FALSE]
  tilde_W <- as.numeric(predict(fit_wp, newdata = df_r1))
  new_wx1 <- set_cols(df_r1, setNames(list(1, tilde_W), c(X, Wprime)))
  new_wx0 <- set_cols(df_r1, setNames(list(0, tilde_W), c(X, Wprime)))
  mu1_t <- predict(fit_out, newdata = new_wx1)
  mu0_t <- predict(fit_out, newdata = new_wx0)
  new_wx1$mu_t <- mu1_t - mu0_t
  theta_w1_z <- unname(coef(lm(as.formula(paste0("mu_t ~ ", paste(Zc_vars, collapse = " + "))), data = new_wx1))[1])
  zeta_stoch_imp <- theta_w1_z - mean((mu01 - mu00))
  delta_stoch_imp <- zeta_z_imp - zeta_stoch_imp

  list(
    zeta_z_imp = zeta_z_imp, delta_z_imp = delta_z_imp,
    zeta_static_imp = zeta_static_imp, delta_static_imp = delta_static_imp,
    zeta_stoch_imp = zeta_stoch_imp, delta_stoch_imp = delta_stoch_imp
  )
}

#' @param B Number of bootstrap replications.
#' @param ... Passed to \code{imputation()}.
#' @return Tibble: estimand, estimate (mean of \code{B} bootstrap draws), std_error, p_value.
bootstrap_imputation <- function(B = 1000, ...) {
  args <- list(...)
  df <- args$df
  if (is.null(df)) stop("A data frame must be provided with argument name 'df'.", call. = FALSE)

  point_estimates <- unlist(do.call(imputation, args))

  boot_results <- replicate(B, {
    boot_indices <- sample.int(nrow(df), nrow(df), replace = TRUE)
    boot_args <- args
    boot_args$df <- df[boot_indices, , drop = FALSE]
    tryCatch(
      unlist(do.call(imputation, boot_args)),
      error = function(e) rep(NA_real_, length(point_estimates))
    )
  }, simplify = TRUE)

  boot_estimates_df <- as_tibble(t(boot_results))
  est_names <- names(point_estimates)
  names(boot_estimates_df) <- est_names

  summary_stats <- tibble(
    estimand = est_names,
    estimate = vapply(est_names, function(nm) mean(boot_estimates_df[[nm]], na.rm = TRUE), numeric(1)),
    std_error = vapply(est_names, function(nm) sd(boot_estimates_df[[nm]], na.rm = TRUE), numeric(1)),
    p_value = vapply(est_names, function(nm) {
      x <- boot_estimates_df[[nm]]
      pvalue_boot_se(x, mean(x, na.rm = TRUE))
    }, numeric(1))
  )

  summary_stats %>%
    mutate(p_value = pvalue(p_value, accuracy = 0.001, add_p = FALSE)) %>%
    select(estimand, estimate, std_error, p_value)
}

# Partially penalized post-double-selection:
# force in all main effects (no ":"), penalize only interactions (contains ":").
pds_predict_hybrid <- function(X, y, d = NULL, lambda_rule = c("1se", "min"),
                               main_idx = NULL, interaction_idx = NULL) {
  stopifnot(nrow(X) == length(y))
  lambda_rule <- match.arg(lambda_rule)
  s <- if (lambda_rule == "1se") "lambda.1se" else "lambda.min"
  p <- ncol(X)
  if (is.null(main_idx)) main_idx <- integer(0)
  if (is.null(interaction_idx)) interaction_idx <- integer(0)
  main_idx <- sort(unique(as.integer(main_idx)))
  main_idx <- main_idx[main_idx >= 1L & main_idx <= p]
  interaction_idx <- sort(unique(as.integer(interaction_idx)))
  interaction_idx <- interaction_idx[interaction_idx >= 1L & interaction_idx <= p]
  if (!length(main_idx) && !length(interaction_idx)) {
    stop("pds_predict_hybrid: no candidate columns supplied.", call. = FALSE)
  }
  cand <- sort(unique(c(main_idx, interaction_idx)))
  Xcand <- X[, cand, drop = FALSE]
  pf <- rep(1, length(cand))
  if (length(main_idx)) {
    pf[match(main_idx, cand)] <- 0
  }
  pick_local <- function(target) {
    fit <- cv.glmnet(Xcand, target, alpha = 1, family = "gaussian", standardize = TRUE,
                     nfolds = 10, penalty.factor = pf)
    cc <- as.numeric(coef(fit, s = s))[-1]
    sel_main_local <- if (length(main_idx)) match(main_idx, cand) else integer(0)
    sel_int_local <- if (length(interaction_idx)) {
      loc <- match(interaction_idx, cand)
      loc[abs(cc[loc]) > 1e-10]
    } else {
      integer(0)
    }
    sort(unique(c(sel_main_local, sel_int_local)))
  }
  sel_y_local <- pick_local(y)
  if (is.null(d)) {
    sel_local <- sel_y_local
  } else {
    stopifnot(nrow(X) == length(d))
    sel_d_local <- pick_local(d)
    sel_local <- union(sel_y_local, sel_d_local)
  }
  sel <- cand[sel_local]
  Xs <- X[, sel, drop = FALSE]
  Xdm <- cbind(1, Xs)
  fit <- lm.fit(x = Xdm, y = y)
  function(Xnew) {
    nc <- ncol(Xnew)
    Xfeat <- if (nc == p) {
      Xnew[, sel, drop = FALSE]
    } else if (nc == p + 1L) {
      Xnew[, sel + 1L, drop = FALSE]
    } else {
      stop("pds_predict_hybrid: Xnew must have ", p, " cols (no intercept) or ", p + 1L,
           " (with intercept); got ", nc, call. = FALSE)
    }
    Xn <- cbind(1, Xfeat)
    drop(Xn %*% fit$coefficients)
  }
}

#' Hybrid LASSO heterogeneity: all main effects forced in; interactions selected by LASSO.
#'
#' Uses penalty.factor = 0 for main-effect columns (no ":"), 1 for interaction columns (":").
#' For the outcome model, double LASSO unions interaction selections from Y and Wprime.
lasso_analysis_hybrid <- function(df, Y, R, X, Zvars, Wvars, Wprime, w_fix = 0,
                                  lambda_rule = c("1se", "min"),
                                  lasso_type = c("double", "single")) {
  lambda_rule <- match.arg(lambda_rule)
  lasso_type <- match.arg(lasso_type)
  safe_lm_coef <- function(formula_obj, data_obj, term = NULL, index = NULL) {
    fit <- tryCatch(lm(formula_obj, data = data_obj), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    cf <- coef(fit)
    if (!is.null(term)) {
      val <- cf[term]
      if (!length(val) || is.na(val)) return(NA_real_)
      return(unname(val))
    }
    if (!is.null(index)) {
      if (length(cf) < index || is.na(cf[index])) return(NA_real_)
      return(unname(cf[index]))
    }
    NA_real_
  }
  mm_train <- function(form, data) {
    M <- model.matrix(form, data = data)
    M[, -1L, drop = FALSE]
  }
  split_idx <- function(cols) {
    list(
      main = which(!grepl(":", cols, fixed = TRUE)),
      inter = which(grepl(":", cols, fixed = TRUE))
    )
  }

  rx <- paste(R, X, sep = ":")
  tau <- safe_lm_coef(as.formula(paste0(Y, " ~ ", R, " * ", X)), df, term = rx)

  form_cond <- as.formula(paste0(Y, " ~ ", R, " * ", X, " * (", paste(Zvars, collapse = " + "), ")"))
  Xmat <- mm_train(form_cond, df)
  idx_cond <- split_idx(colnames(Xmat))
  yvec <- df[[Y]]
  d_wprime <- as.numeric(df[[Wprime]])
  pred_tau <- pds_predict_hybrid(Xmat, yvec, d = NULL, lambda_rule = lambda_rule,
                                 main_idx = idx_cond$main, interaction_idx = idx_cond$inter)

  new_11 <- set_cols(df, setNames(list(1, 1), c(R, X)))
  new_10 <- set_cols(df, setNames(list(1, 0), c(R, X)))
  new_01 <- set_cols(df, setNames(list(0, 1), c(R, X)))
  new_00 <- set_cols(df, setNames(list(0, 0), c(R, X)))
  mu11 <- pred_tau(model.matrix(form_cond, data = new_11))
  mu10 <- pred_tau(model.matrix(form_cond, data = new_10))
  mu01 <- pred_tau(model.matrix(form_cond, data = new_01))
  mu00 <- pred_tau(model.matrix(form_cond, data = new_00))
  zeta_z_lasso <- mean((mu11 - mu10) - (mu01 - mu00))
  delta_z_lasso <- tau - zeta_z_lasso

  Zc_vars <- paste0(Zvars, "_c")
  for (i in seq_along(Zvars)) df[[Zc_vars[i]]] <- df[[Zvars[i]]] - mean(df[[Zvars[i]]], na.rm = TRUE)

  form_out <- as.formula(paste0(Y, " ~ ", R, " * ", X, " * ", Wprime, " * (", paste(c(Zvars, Wvars), collapse = " + "), ")"))
  X_out <- mm_train(form_out, df)
  idx_out <- split_idx(colnames(X_out))
  d_out <- if (lasso_type == "double") d_wprime else NULL
  pred_out <- pds_predict_hybrid(X_out, yvec, d = d_out, lambda_rule = lambda_rule,
                                 main_idx = idx_out$main, interaction_idx = idx_out$inter)

  new_x1 <- set_cols(df, setNames(list(1, w_fix), c(X, Wprime)))
  new_x0 <- set_cols(df, setNames(list(0, w_fix), c(X, Wprime)))
  mu_x1 <- pred_out(model.matrix(form_out, data = new_x1))
  mu_x0 <- pred_out(model.matrix(form_out, data = new_x0))
  mu_w <- mu_x1 - mu_x0
  dat_mu <- cbind(df, mu_w = mu_w)
  form_mu <- as.formula(paste0("mu_w ~ ", R, " + ", paste(Zc_vars, collapse = " + ")))
  zeta_static_lasso <- safe_lm_coef(form_mu, dat_mu, term = R)
  delta_static_lasso <- zeta_z_lasso - zeta_static_lasso

  fit_wp <- tryCatch(
    lm(as.formula(paste0(Wprime, " ~ ", paste(Zvars, collapse = " + "))), data = df[df[[R]] == 0, , drop = FALSE]),
    error = function(e) NULL
  )
  df_r1 <- df[df[[R]] == 1, , drop = FALSE]
  tilde_W <- if (!is.null(fit_wp)) as.numeric(predict(fit_wp, newdata = df_r1)) else rep(mean(df[[Wprime]], na.rm = TRUE), nrow(df_r1))
  new_wx1 <- set_cols(df_r1, setNames(list(1, tilde_W), c(X, Wprime)))
  new_wx0 <- set_cols(df_r1, setNames(list(0, tilde_W), c(X, Wprime)))
  mu1_t <- pred_out(model.matrix(form_out, data = new_wx1))
  mu0_t <- pred_out(model.matrix(form_out, data = new_wx0))
  new_wx1$mu_t <- mu1_t - mu0_t
  theta_w1_z <- safe_lm_coef(as.formula(paste0("mu_t ~ ", paste(Zc_vars, collapse = " + "))), new_wx1, index = 1)
  zeta_stoch_lasso <- theta_w1_z - mean((mu01 - mu00))
  delta_stoch_lasso <- zeta_z_lasso - zeta_stoch_lasso

  list(
    zeta_z_lasso = zeta_z_lasso, delta_z_lasso = delta_z_lasso,
    zeta_static_lasso = zeta_static_lasso, delta_static_lasso = delta_static_lasso,
    zeta_stoch_lasso = zeta_stoch_lasso, delta_stoch_lasso = delta_stoch_lasso
  )
}

# Hybrid selection diagnostics:
# counts only selected interactions; main effects are always retained by design.
lasso_selection_sizes_hybrid <- function(df, Y, R, X, Zvars, Wvars, Wprime,
                                         lambda_rule = c("1se", "min"),
                                         lasso_type = c("double", "single")) {
  lambda_rule <- match.arg(lambda_rule)
  lasso_type <- match.arg(lasso_type)
  s <- if (lambda_rule == "1se") "lambda.1se" else "lambda.min"
  mm_train <- function(form, data) {
    M <- model.matrix(form, data = data)
    M[, -1L, drop = FALSE]
  }
  pick_sel <- function(Xmat, yvec, dvec) {
    cols <- colnames(Xmat)
    main_idx <- which(!grepl(":", cols, fixed = TRUE))
    int_idx <- which(grepl(":", cols, fixed = TRUE))
    cand <- sort(unique(c(main_idx, int_idx)))
    Xcand <- Xmat[, cand, drop = FALSE]
    pf <- rep(1, length(cand))
    if (length(main_idx)) pf[match(main_idx, cand)] <- 0
    pick_int <- function(target) {
      fit <- cv.glmnet(Xcand, target, alpha = 1, family = "gaussian", standardize = TRUE,
                       nfolds = 10, penalty.factor = pf)
      cc <- as.numeric(coef(fit, s = s))[-1]
      if (!length(int_idx)) return(integer(0))
      loc <- match(int_idx, cand)
      int_idx[abs(cc[loc]) > 1e-10]
    }
    sel_y <- pick_int(yvec)
    if (is.null(dvec)) {
      sel_d <- integer(0)
      sel <- sel_y
    } else {
      sel_d <- pick_int(dvec)
      sel <- union(sel_y, sel_d)
    }
    c(
      n_main_forced = length(main_idx),
      n_int_sel_y = length(sel_y),
      n_int_sel_d = length(sel_d),
      n_int_sel_union = length(sel)
    )
  }

  yvec <- df[[Y]]
  d_wprime <- as.numeric(df[[Wprime]])
  form_cond <- as.formula(paste0(Y, " ~ ", R, " * ", X, " * (", paste(Zvars, collapse = " + "), ")"))
  form_out <- as.formula(paste0(Y, " ~ ", R, " * ", X, " * ", Wprime, " * (", paste(c(Zvars, Wvars), collapse = " + "), ")"))
  X_cond <- mm_train(form_cond, df)
  X_out <- mm_train(form_out, df)
  d_out <- if (lasso_type == "double") d_wprime else NULL
  cond_sizes <- pick_sel(X_cond, yvec, NULL)
  out_sizes <- pick_sel(X_out, yvec, d_out)

  tibble(
    n_cond_main_forced = unname(cond_sizes["n_main_forced"]),
    n_sel_cond_y = unname(cond_sizes["n_int_sel_y"]),
    n_sel_cond_d = unname(cond_sizes["n_int_sel_d"]),
    n_sel_cond_union = unname(cond_sizes["n_int_sel_union"]),
    n_out_main_forced = unname(out_sizes["n_main_forced"]),
    n_sel_out_y = unname(out_sizes["n_int_sel_y"]),
    n_sel_out_d = unname(out_sizes["n_int_sel_d"]),
    n_sel_out_union = unname(out_sizes["n_int_sel_union"])
  )
}

#' In-sample selection diagnostics for hybrid LASSO specifications.
#'
#' Builds one row per combination of \code{lasso_type} and \code{lambda_rule}.
#'
#' @param lasso_types Character vector passed through, typically \code{c("double", "single")}.
#' @return A tibble sorted by \code{lasso_type} and \code{lambda_rule}.
collect_lasso_diagnostics <- function(df, Y, R, X, Zvars, Wvars, Wprime,
                                      lasso_types = c("double", "single")) {
  lasso_types <- match.arg(lasso_types, c("double", "single"), several.ok = TRUE)
  rules <- c("min", "1se")
  rows <- list()
  k <- 0L
  for (rule in rules) {
    for (lt in lasso_types) {
      k <- k + 1L
      sz <- lasso_selection_sizes_hybrid(df, Y, R, X, Zvars, Wvars, Wprime,
                                         lambda_rule = rule, lasso_type = lt)
      rows[[k]] <- tibble(lasso_type = lt, lambda_rule = rule) %>% bind_cols(sz)
    }
  }
  out <- bind_rows(rows)
  out %>% arrange(lasso_type, match(lambda_rule, c("min", "1se")))
}

#' @param B Number of bootstrap replications.
#' @param hybrid Deprecated; hybrid LASSO is always used.
#' @param ... Passed to \code{lasso_analysis_hybrid()}.
#' @return Tibble: estimand, estimate (mean of \code{B} bootstrap draws), std_error (bootstrap SD), p_value.
bootstrap_lasso <- function(B = 1000, hybrid = TRUE, ...) {
  if (!isTRUE(hybrid)) warning("Standard LASSO was removed; using hybrid LASSO.", call. = FALSE)
  model_args <- list(...)
  df <- model_args$df
  if (is.null(df)) stop("A data frame must be provided with the argument name 'df'.", call. = FALSE)

  point_estimates <- unlist(do.call(lasso_analysis_hybrid, model_args))

  boot_results <- replicate(B, {
    boot_indices <- sample.int(nrow(df), nrow(df), replace = TRUE)
    boot_args <- model_args
    boot_args$df <- df[boot_indices, , drop = FALSE]
    tryCatch(
      unlist(do.call(lasso_analysis_hybrid, boot_args)),
      error = function(e) rep(NA_real_, length(point_estimates))
    )
  }, simplify = TRUE)
  boot_estimates_df <- as_tibble(t(boot_results))
  est_names <- names(point_estimates)
  names(boot_estimates_df) <- est_names

  summary_stats <- tibble(
    estimand = est_names,
    estimate = vapply(est_names, function(nm) mean(boot_estimates_df[[nm]], na.rm = TRUE), numeric(1)),
    std_error = vapply(est_names, function(nm) sd(boot_estimates_df[[nm]], na.rm = TRUE), numeric(1)),
    p_value = vapply(est_names, function(nm) {
      x <- boot_estimates_df[[nm]]
      pvalue_boot_se(x, mean(x, na.rm = TRUE))
    }, numeric(1))
  )

  summary_stats %>%
    mutate(p_value = pvalue(p_value, accuracy = 0.001, add_p = FALSE)) %>%
    select(estimand, estimate, std_error, p_value)
}

#' Compare hybrid LASSO diagnostics for lambda.min vs lambda.1se on identical bootstrap samples.
#'
#' @param B Number of bootstrap replications.
#' @param lasso_types Character vector: \code{"double"}, \code{"single"}, or both.
#' @param ... Passed to \code{lasso_analysis_hybrid()} and \code{lasso_selection_sizes_hybrid()}.
#'   If \code{lasso_type} is included, it is ignored in favor of \code{lasso_types}.
#' @return A list with point estimates, bootstrap draws, and bootstrap summary by
#'   \code{lambda_rule} and \code{lasso_type}.
bootstrap_lasso_diagnostics <- function(B = 200, lasso_types = c("double", "single"), ...) {
  model_args <- list(...)
  df <- model_args$df
  if (is.null(df)) stop("A data frame must be provided with the argument name 'df'.", call. = FALSE)
  model_args$lasso_type <- NULL
  lasso_types <- match.arg(lasso_types, c("double", "single"), several.ok = TRUE)
  sel_args_base <- model_args
  sel_args_base$w_fix <- NULL

  rules <- c("min", "1se")
  point <- bind_rows(lapply(lasso_types, function(lt) {
    bind_rows(lapply(rules, function(rule) {
      args_rule <- modifyList(model_args, list(lambda_rule = rule, lasso_type = lt))
      est <- unlist(do.call(lasso_analysis_hybrid, args_rule))
      sizes_rule <- modifyList(sel_args_base, list(lambda_rule = rule, lasso_type = lt))
      sizes <- do.call(lasso_selection_sizes_hybrid, sizes_rule)
      tibble(lasso_type = lt, lambda_rule = rule) %>%
        bind_cols(tibble::as_tibble_row(as.list(est)), sizes)
    }))
  }))

  boot <- bind_rows(lapply(seq_len(B), function(rep_i) {
    idx <- sample.int(nrow(df), nrow(df), replace = TRUE)
    boot_df <- df[idx, , drop = FALSE]
    bind_rows(lapply(lasso_types, function(lt) {
      bind_rows(lapply(rules, function(rule) {
        args_rule <- modifyList(model_args, list(df = boot_df, lambda_rule = rule, lasso_type = lt))
        sizes_rule <- modifyList(sel_args_base, list(df = boot_df, lambda_rule = rule, lasso_type = lt))
        est <- tryCatch(unlist(do.call(lasso_analysis_hybrid, args_rule)), error = function(e) NULL)
        sizes <- tryCatch(do.call(lasso_selection_sizes_hybrid, sizes_rule), error = function(e) NULL)
        if (is.null(est) || is.null(sizes)) return(NULL)
        tibble(rep = rep_i, lasso_type = lt, lambda_rule = rule) %>%
          bind_cols(tibble::as_tibble_row(as.list(est)), sizes)
      }))
    }))
  }))

  ref_args <- modifyList(model_args, list(lambda_rule = "min", lasso_type = lasso_types[1]))
  est_cols <- names(unlist(do.call(lasso_analysis_hybrid, ref_args)))
  size_cols <- c("n_cond_main_forced", "n_sel_cond_y", "n_sel_cond_d", "n_sel_cond_union",
                 "n_out_main_forced", "n_sel_out_y", "n_sel_out_d", "n_sel_out_union")
  summary <- boot %>%
    group_by(lasso_type, lambda_rule) %>%
    summarise(
      across(all_of(c(est_cols, size_cols)), list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))),
      .groups = "drop"
    )

  list(point = point, bootstrap = boot, summary = summary)
}
