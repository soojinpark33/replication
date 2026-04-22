library(glmnet)

# --- POC (linear models; simulation DGP) ---
poc <- function(df) {
  fit_Y0 <- lm(Y ~ R * X * Z, data = df)
  cf0 <- coef(fit_Y0)
  phi_rx <- cf0["R:X"]
  phi_rxz <- cf0["R:X:Z"]
  tau_z_poc <- phi_rx + phi_rxz * mean(df$Z)

  fit_Y <- lm(Y ~ R * X * Z + W + Wprime + X:Wprime + R:X:Wprime, data = df)
  coefs <- coef(fit_Y)
  gamma_rx <- coefs["R:X"]
  gamma_rxwprime <- coefs["R:X:Wprime"]
  gamma_rxz <- coefs["R:X:Z"]
  w_prime <- 1
  zeta_static_poc <- gamma_rx + gamma_rxwprime * w_prime + gamma_rxz * mean(df$Z)

  sum_z_Ewprime_pz <- {
    terms <- vapply(sort(unique(df$Z)), function(z) {
      pz <- mean(df$Z == z)
      idx <- df$R == 0 & df$Z == z
      if (!any(idx)) NA_real_ else pz * mean(df$Wprime[idx])
    }, numeric(1))
    if (anyNA(terms)) NA_real_ else sum(terms)
  }
  zeta_stoch_poc <- gamma_rx + gamma_rxwprime * sum_z_Ewprime_pz + gamma_rxz * mean(df$Z)

  c(
    tau_rxz = unname(as.numeric(tau_z_poc)),
    zeta_static_poc = unname(as.numeric(zeta_static_poc)),
    zeta_stoch_poc = unname(as.numeric(zeta_stoch_poc))
  )
}

Ewprime_r0_z <- function(df) {
  df_r0 <- df[df$R == 0, , drop = FALSE]
  if (!nrow(df_r0)) return(rep(mean(df$Wprime), nrow(df)))
  m0 <- mean(df_r0$Wprime)
  mz <- tapply(df_r0$Wprime, df_r0$Z, mean)
  out <- unname(mz[as.character(df$Z)])
  out[is.na(out)] <- m0
  out
}

# --- Imputation (full linear models) ---
imp <- function(df) {
  fit_cond <- lm(Y ~ R * X * Z, data = df)
  new_11 <- transform(df, R = 1, X = 1)
  new_10 <- transform(df, R = 1, X = 0)
  new_01 <- transform(df, R = 0, X = 1)
  new_00 <- transform(df, R = 0, X = 0)
  mu11 <- predict(fit_cond, newdata = new_11)
  mu10 <- predict(fit_cond, newdata = new_10)
  mu01 <- predict(fit_cond, newdata = new_01)
  mu00 <- predict(fit_cond, newdata = new_00)
  tau_z_imp <- mean((mu11 - mu10) - (mu01 - mu00))

  df$Zc <- df$Z - mean(df$Z)
  fit_out <- lm(Y ~ R * X * Zc + W + Wprime + X:Wprime + R:X:Wprime, data = df)
  w_fix <- 1
  newdata_x1 <- transform(df, X = 1, Wprime = w_fix)
  newdata_x0 <- transform(df, X = 0, Wprime = w_fix)
  mu_x1 <- predict(fit_out, newdata = newdata_x1)
  mu_x0 <- predict(fit_out, newdata = newdata_x0)
  mu_w <- mu_x1 - mu_x0
  zeta_static_imp <- unname(coef(lm(mu_w ~ R + Zc, data = df))["R"])

  tilde_W <- Ewprime_r0_z(df)
  newdata_wx1 <- subset(transform(df, X = 1, Wprime = tilde_W), R == 1)
  newdata_wx0 <- subset(transform(df, X = 0, Wprime = tilde_W), R == 1)
  mu1_t <- predict(fit_out, newdata = newdata_wx1)
  mu0_t <- predict(fit_out, newdata = newdata_wx0)
  newdata_wx1$mu_t <- mu1_t - mu0_t
  theta_w1_z <- unname(coef(lm(mu_t ~ Zc, data = newdata_wx1))[1])
  zeta_stoch_imp <- theta_w1_z - mean((mu01 - mu00))

  c(
    tau_z_imp = tau_z_imp,
    zeta_static_imp = zeta_static_imp,
    zeta_stoch_imp = zeta_stoch_imp
  )
}

# Hybrid post-double-selection:
# force-in main effects (no ":"), penalize interactions (contains ":").
pds_predict <- function(X, y, d = NULL, lambda_rule = c("min", "1se"),
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
  cand <- sort(unique(c(main_idx, interaction_idx)))
  if (!length(cand)) stop("pds_predict: no candidate columns supplied.", call. = FALSE)
  Xcand <- X[, cand, drop = FALSE]
  pf <- rep(1, length(cand))
  if (length(main_idx)) pf[match(main_idx, cand)] <- 0

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
    Xn <- cbind(1, Xnew[, sel, drop = FALSE])
    drop(Xn %*% fit$coefficients)
  }
}

lasso <- function(df, lambda_rule = c("min", "1se")) {
  lambda_rule <- match.arg(lambda_rule)
  form <- as.formula("Y ~ R * X * Z")
  mm <- function(d) model.matrix(form, data = d)[, -1, drop = FALSE]

  Xmat <- mm(df)
  yvec <- df$Y
  d_wprime <- as.numeric(df$Wprime)
  idx_tau_main <- which(!grepl(":", colnames(Xmat), fixed = TRUE))
  idx_tau_inter <- which(grepl(":", colnames(Xmat), fixed = TRUE))
  pred_tau <- pds_predict(Xmat, yvec, d = NULL, lambda_rule = lambda_rule,
                          main_idx = idx_tau_main, interaction_idx = idx_tau_inter)

  new_11 <- transform(df, R = 1, X = 1)
  new_10 <- transform(df, R = 1, X = 0)
  new_01 <- transform(df, R = 0, X = 1)
  new_00 <- transform(df, R = 0, X = 0)
  mu11 <- pred_tau(mm(new_11))
  mu10 <- pred_tau(mm(new_10))
  mu01 <- pred_tau(mm(new_01))
  mu00 <- pred_tau(mm(new_00))
  tau_z_imp <- mean((mu11 - mu10) - (mu01 - mu00))

  df$Zc <- df$Z - mean(df$Z)
  form1 <- as.formula("Y ~ R * X * Z * W * Wprime")
  X_out <- model.matrix(form1, data = df)[, -1, drop = FALSE]
  y_out <- df$Y
  idx_out_main <- which(!grepl(":", colnames(X_out), fixed = TRUE))
  idx_out_inter <- which(grepl(":", colnames(X_out), fixed = TRUE))
  pred_out <- pds_predict(X_out, y_out, d_wprime, lambda_rule = lambda_rule,
                          main_idx = idx_out_main, interaction_idx = idx_out_inter)

  w_fix <- 1
  newdata_x1 <- transform(df, X = 1, Wprime = w_fix)
  newdata_x0 <- transform(df, X = 0, Wprime = w_fix)
  mu_x1 <- pred_out(model.matrix(form1, data = newdata_x1)[, -1, drop = FALSE])
  mu_x0 <- pred_out(model.matrix(form1, data = newdata_x0)[, -1, drop = FALSE])
  mu_w <- mu_x1 - mu_x0
  zeta_static_imp <- unname(coef(lm(mu_w ~ R + Zc, data = df))["R"])

  tilde_W <- Ewprime_r0_z(df)
  newdata_wx1 <- subset(transform(df, X = 1, Wprime = tilde_W), R == 1)
  newdata_wx0 <- subset(transform(df, X = 0, Wprime = tilde_W), R == 1)
  mu1_t <- pred_out(model.matrix(form1, data = newdata_wx1)[, -1, drop = FALSE])
  mu0_t <- pred_out(model.matrix(form1, data = newdata_wx0)[, -1, drop = FALSE])
  newdata_wx1$mu_t <- mu1_t - mu0_t
  theta_w1_z <- unname(coef(lm(mu_t ~ Zc, data = newdata_wx1))[1])
  zeta_stoch_imp <- theta_w1_z - mean((mu01 - mu00))

  c(
    tau_z_imp = tau_z_imp,
    zeta_static_imp = zeta_static_imp,
    zeta_stoch_imp = zeta_stoch_imp
  )
}
