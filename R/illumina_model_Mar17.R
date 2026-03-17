# ==============================================================================
# 0. Utility functions
# ==============================================================================

# Safe Beta simulator
rbeta_safe <- function(n, alpha, beta, eps = 1e-12) {
  alpha <- pmax(alpha, eps)
  beta  <- pmax(beta, eps)
  y <- stats::rbeta(n, shape1 = alpha, shape2 = beta)
  pmin(pmax(y, eps), 1 - eps)
}

# Delta-method CI helper on the logit scale
delta_ci_logit <- function(eta, X, V, conf_level = 0.95) {
  z_alpha <- qnorm(0.5 + conf_level / 2)
  n <- nrow(X)
  out <- matrix(NA_real_, n, 2)
  colnames(out) <- c("lower", "upper")

  for (j in seq_len(n)) {
    xj <- matrix(X[j, ], ncol = 1)
    vj <- drop(t(xj) %*% V %*% xj)
    if (!is.finite(vj) || vj < 0) {
      out[j, ] <- c(NA_real_, NA_real_)
      next
    }
    se <- sqrt(vj)
    ci_eta <- eta[j] + c(-1, 1) * z_alpha * se
    out[j, ] <- plogis(ci_eta)
  }
  out
}

# ==============================================================================
# 1. Likelihood Functions
# ==============================================================================

# Female-only Beta likelihood
beta_loglik_y <- function(par, y, z, w, eps = 1e-12) {
  pz <- ncol(z)
  pw <- ncol(w)

  a <- par[1:pz]
  b <- par[(pz + 1):(pz + pw)]
  kappa <- par[pz + pw + 1]

  nu <- plogis(drop(z %*% a))
  pi <- plogis(drop(w %*% b))

  q2  <- nu * (1 - pi)
  q1  <- nu * pi * (1 - nu) + (1 - nu) * (1 - pi)
  eta <- 0.5 * q1 + q2
  eta <- pmin(pmax(eta, eps), 1 - eps)

  y <- pmin(pmax(y, eps), 1 - eps)

  alpha <- kappa * eta
  beta  <- kappa * (1 - eta)

  sum(dbeta(y, shape1 = alpha, shape2 = beta, log = TRUE))
}

# ------------------------------------------------------------------------------
# Combined likelihood for females and males
# Females: Beta with XCI mixture mean eta_f
# Males:   Beta with mean nu_m
# ------------------------------------------------------------------------------
neg_beta_loglik_combined <- function(par, y, z, w, fem_idx, eps = 1e-12) {
  pz <- ncol(z)
  pw <- ncol(w)

  a <- par[1:pz]
  b <- par[(pz + 1):(pz + pw)]
  kappa <- par[pz + pw + 1]

  fem_idx  <- as.logical(fem_idx)
  male_idx <- !fem_idx

  y <- pmin(pmax(y, eps), 1 - eps)

  ll_fem <- 0
  if (any(fem_idx)) {
    z_f <- z[fem_idx, , drop = FALSE]
    w_f <- w[fem_idx, , drop = FALSE]
    y_f <- y[fem_idx]

    nu_f <- plogis(drop(z_f %*% a))
    pi_f <- plogis(drop(w_f %*% b))

    q2  <- nu_f * (1 - pi_f)
    q1  <- nu_f * pi_f * (1 - nu_f) + (1 - nu_f) * (1 - pi_f)
    eta_f <- 0.5 * q1 + q2
    eta_f <- pmin(pmax(eta_f, eps), 1 - eps)

    alpha_f <- kappa * eta_f
    beta_f  <- kappa * (1 - eta_f)

    ll_fem <- sum(dbeta(y_f, shape1 = alpha_f, shape2 = beta_f, log = TRUE))
  }

  ll_male <- 0
  if (any(male_idx)) {
    z_m <- z[male_idx, , drop = FALSE]
    y_m <- y[male_idx]

    nu_m <- plogis(drop(z_m %*% a))
    eta_m <- pmin(pmax(nu_m, eps), 1 - eps)

    alpha_m <- kappa * eta_m
    beta_m  <- kappa * (1 - eta_m)

    ll_male <- sum(dbeta(y_m, shape1 = alpha_m, shape2 = beta_m, log = TRUE))
  }

  -(ll_fem + ll_male)
}

# ==============================================================================
# 2. Parametric bootstrap for Beta model
# ==============================================================================

bootstrap_fit_and_ci_beta <- function(
    par_hat, y, z, w, fem_idx,
    theta_start,
    conf_level = 0.95,
    B = 300,
    seed = 1,
    eps = 1e-12,
    control = list(maxit = 5000, factr = 1e7),
    verbose = FALSE
) {
  set.seed(seed)

  fem_idx  <- as.logical(fem_idx)
  male_idx <- !fem_idx

  pz <- ncol(z)
  pw <- ncol(w)
  p_all <- pz + pw + 1

  a_hat     <- par_hat[1:pz]
  b_hat     <- par_hat[(pz + 1):(pz + pw)]
  kappa_hat <- par_hat[pz + pw + 1]

  nu0 <- plogis(drop(z %*% a_hat))
  pi0 <- plogis(drop(w %*% b_hat))

  q2  <- nu0 * (1 - pi0)
  q1  <- nu0 * pi0 * (1 - nu0) + (1 - nu0) * (1 - pi0)
  eta_f <- 0.5 * q1 + q2
  eta_f <- pmin(pmax(eta_f, eps), 1 - eps)

  eta_m <- pmin(pmax(nu0, eps), 1 - eps)

  nu_boot  <- matrix(NA_real_, nrow(z), B)
  pi_boot  <- matrix(NA_real_, nrow(w), B)
  par_boot <- matrix(NA_real_, nrow = B, ncol = p_all)

  n_ok <- 0L

  for (b in seq_len(B)) {
    y_star <- numeric(length(y))

    if (any(fem_idx)) {
      idx <- which(fem_idx)
      alpha_f <- kappa_hat * eta_f[idx]
      beta_f  <- kappa_hat * (1 - eta_f[idx])
      y_star[idx] <- rbeta_safe(length(idx), alpha_f, beta_f, eps = eps)
    }

    if (any(male_idx)) {
      idx <- which(male_idx)
      alpha_m <- kappa_hat * eta_m[idx]
      beta_m  <- kappa_hat * (1 - eta_m[idx])
      y_star[idx] <- rbeta_safe(length(idx), alpha_m, beta_m, eps = eps)
    }

    fit_b <- tryCatch(
      optim(
        par       = theta_start,
        fn        = neg_beta_loglik_combined,
        y         = y_star,
        z         = z,
        w         = w,
        fem_idx   = fem_idx,
        eps       = eps,
        method    = "L-BFGS-B",
        lower     = c(rep(-Inf, pz), rep(-Inf, pw), 1e-8),
        upper     = c(rep( Inf, pz), rep( Inf, pw), Inf),
        control   = control
      ),
      error = function(e) NULL
    )

    if (is.null(fit_b) || is.null(fit_b$par) || any(!is.finite(fit_b$par))) {
      if (verbose) message("Bootstrap refit failed at b = ", b)
      next
    }

    n_ok <- n_ok + 1L
    par_b <- fit_b$par
    par_boot[n_ok, ] <- par_b

    a_b <- par_b[1:pz]
    b_b <- par_b[(pz + 1):(pz + pw)]

    nu_boot[, n_ok] <- plogis(drop(z %*% a_b))
    pi_boot[, n_ok] <- plogis(drop(w %*% b_b))
  }

  if (n_ok < max(50, 0.2 * B)) {
    warning("Too few successful bootstrap refits (n_ok = ", n_ok, "). Results may be unstable.")
  }

  if (n_ok == 0L) {
    return(list(
      nu_ci_boot = matrix(NA_real_, nrow(z), 2, dimnames = list(NULL, c("lower", "upper"))),
      pi_ci_boot = matrix(NA_real_, nrow(w), 2, dimnames = list(NULL, c("lower", "upper"))),
      par_boot   = NULL,
      vc_boot    = NULL,
      se_boot    = NULL,
      B          = B,
      n_ok       = n_ok,
      seed       = seed
    ))
  }

  par_boot <- par_boot[seq_len(n_ok), , drop = FALSE]
  nu_boot  <- nu_boot[, seq_len(n_ok), drop = FALSE]
  pi_boot  <- pi_boot[, seq_len(n_ok), drop = FALSE]

  alpha_q <- (1 - conf_level) / 2
  probs <- c(alpha_q, 1 - alpha_q)

  nu_ci <- t(apply(nu_boot, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 10) return(c(NA_real_, NA_real_))
    quantile(v, probs = probs, names = FALSE, type = 8)
  }))

  pi_ci <- t(apply(pi_boot, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 10) return(c(NA_real_, NA_real_))
    quantile(v, probs = probs, names = FALSE, type = 8)
  }))

  colnames(nu_ci) <- c("lower", "upper")
  colnames(pi_ci) <- c("lower", "upper")

  vc_boot <- if (n_ok >= 2) stats::cov(par_boot) else NULL
  se_boot <- if (!is.null(vc_boot)) sqrt(pmax(diag(vc_boot), 0)) else NULL

  list(
    nu_ci_boot = nu_ci,
    pi_ci_boot = pi_ci,
    par_boot   = par_boot,
    vc_boot    = vc_boot,
    se_boot    = se_boot,
    B          = B,
    n_ok       = n_ok,
    seed       = seed
  )
}

# ==============================================================================
# 3. Main fitting function for Illumina beta values
# ==============================================================================

#' Fit the LVE-X model for Illumina methylation beta values at a single CpG
#'
#' Fits the LVE-X Beta model at one CpG/site using Illumina-array methylation
#' beta values (continuous outcomes in (0,1)). The function estimates
#' covariate-dependent parameters for the active-allele methylation level
#' \eqn{\nu} (Xa) and the inactive-allele methylation level \eqn{\pi} (Xi),
#' while allowing female samples to follow the LVE-X mixture mean structure and
#' male samples to follow a simpler Beta model with mean \eqn{\nu}.
#'
#' Confidence intervals are obtained using either:
#' \itemize{
#'   \item a Hessian-based delta method, or
#'   \item a parametric bootstrap fallback if Hessian-based inference fails.
#' }
#'
#' @param i Integer. Column index of the CpG/site to fit, i.e. which column of `BetaMat`.
#' @param BetaMat Numeric matrix (n x p). Illumina methylation beta values;
#'   rows are samples and columns are CpG sites. Values should lie in 0 and 1.
#' @param family Character. Model family; currently only `"Beta"` is supported.
#' @param z Optional numeric matrix/data.frame (n x q). Covariates for
#'   \eqn{\nu} (Xa). If `NULL`, an intercept-only model is used.
#' @param w Optional numeric matrix/data.frame (n x r). Covariates for
#'   \eqn{\pi} (Xi). If `NULL`, an intercept-only model is used.
#' @param fem_idx Integer or logical vector indicating female samples.
#'   Must be supplied.
#' @param theta_start Optional numeric vector of starting values for optimization.
#'   Its required length is `ncol(z) + ncol(w) + 1` after intercepts are added
#'   internally. If `NULL`, defaults to zeros for regression coefficients and 10
#'   for `kappa`.
#' @param conf_level Numeric in (0,1). Confidence level for intervals.
#'   Default is 0.95.
#' @param eps Numeric. Small value used to clip beta values and probabilities
#'   away from 0 and 1 to avoid numerical issues. Default is `1e-12`.
#' @param bootstrap_on_fail Logical. If `TRUE`, performs parametric bootstrap
#'   when Hessian-based inference fails. Default is `TRUE`.
#' @param bootstrap_B Integer. Number of bootstrap replicates. Default is 300.
#' @param bootstrap_seed Integer. Random seed for bootstrap reproducibility.
#'   Default is 1.
#' @param bootstrap_verbose Logical. If `TRUE`, prints bootstrap messages.
#'   Default is `FALSE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{point_estimates}{Named vector of parameter estimates.}
#'   \item{standard_errors}{Named vector of standard errors, if available.}
#'   \item{covariance_matrix}{Variance-covariance matrix of the parameter estimates, if available.}
#'   \item{pi_hat}{Fitted sample-level \eqn{\pi} values.}
#'   \item{nu_hat}{Fitted sample-level \eqn{\nu} values.}
#'   \item{pi_ci}{Confidence interval matrix for \eqn{\pi}, with columns `lower` and `upper`.}
#'   \item{nu_ci}{Confidence interval matrix for \eqn{\nu}, with columns `lower` and `upper`.}
#'   \item{conf_level}{Confidence level used.}
#'   \item{ci_method}{Method used for intervals: `"delta-hessian"`, `"parametric-bootstrap"`, or `"none"`.}
#'   \item{bootstrap_details}{Bootstrap output if bootstrap was used, including parameter draws.}
#' }
#'
#' @details
#' Internally, the function:
#' \itemize{
#'   \item adds intercept columns to `z` and `w`,
#'   \item optimizes the negative log-likelihood using `stats::optim()` with
#'   method `"L-BFGS-B"` and constraint `kappa > 0`,
#'   \item attempts Hessian-based inference using `numDeriv::hessian()`,
#'   \item falls back to parametric bootstrap if Hessian-based inference is unavailable or unstable.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 20; p <- 5
#' BetaMat <- matrix(runif(n * p, 0.05, 0.95), n, p)
#' fem_idx <- sample(c(TRUE, FALSE), n, replace = TRUE)
#' z <- data.frame(age = rnorm(n))
#' w <- data.frame(batch = rnorm(n))
#'
#' fit <- LVEX_fit_illumina_beta(
#'   i = 2,
#'   BetaMat = BetaMat,
#'   z = z,
#'   w = w,
#'   fem_idx = fem_idx
#' )
#' str(fit)
#' }
#'
#' @export
LVEX_fit_illumina_beta <- function(
    i,
    BetaMat,
    family = "Beta",
    z = NULL,
    w = NULL,
    fem_idx = NULL,
    theta_start = NULL,
    conf_level = 0.95,
    eps = 1e-12,
    bootstrap_on_fail = TRUE,
    bootstrap_B = 300,
    bootstrap_seed = 1,
    bootstrap_verbose = FALSE
) {
  if (is.null(fem_idx)) stop("Error: argument 'fem_idx' must be provided.")
  if (family != "Beta") stop("Error: 'family' type is not correct.")

  z_name <- if (!is.null(z)) colnames(z) else NULL
  w_name <- if (!is.null(w)) colnames(w) else NULL

  n_rows <- nrow(BetaMat)

  ## z
  if (is.null(z)) {
    z <- matrix(1, nrow = n_rows, ncol = 1L)
    colnames(z) <- "(Intercept)"
  } else {
    z <- as.matrix(z)
    if (!is.numeric(z)) stop("`z` must be numeric.")
    z <- cbind(1, z)
  }

  ## w
  if (is.null(w)) {
    w <- matrix(1, nrow = n_rows, ncol = 1L)
    colnames(w) <- "(Intercept)"
  } else {
    w <- as.matrix(w)
    if (!is.numeric(w)) stop("`w` must be numeric.")
    w <- cbind(1, w)
  }

  pz <- ncol(z)
  pw <- ncol(w)

  if (is.null(theta_start)) {
    theta_start <- c(rep(0, pz), rep(0, pw), 10)
  }

  y_i <- BetaMat[, i]
  y_i <- pmin(pmax(y_i, eps), 1 - eps)

  fit <- optim(
    par       = theta_start,
    fn        = neg_beta_loglik_combined,
    y         = y_i,
    z         = z,
    w         = w,
    fem_idx   = fem_idx,
    eps       = eps,
    method    = "L-BFGS-B",
    lower     = c(rep(-Inf, pz), rep(-Inf, pw), 1e-8),
    upper     = c(rep( Inf, pz), rep( Inf, pw), Inf),
    control   = list(maxit = 5000, factr = 1e7)
  )

  a_names <- if (!is.null(z_name) && any(nzchar(z_name))) {
    paste0("a_", c("Intercept", z_name))
  } else {
    paste0("a", seq_len(pz) - 1L)
  }

  b_names <- if (!is.null(w_name) && any(nzchar(w_name))) {
    paste0("b_", c("Intercept", w_name))
  } else {
    paste0("b", seq_len(pw) - 1L)
  }

  est <- setNames(fit$par, c(a_names, b_names, "kappa"))

  a_hat <- fit$par[1:pz]
  b_hat <- fit$par[(pz + 1):(pz + pw)]
  eta_nu <- drop(z %*% a_hat)
  eta_pi <- drop(w %*% b_hat)
  nu_hat <- plogis(eta_nu)
  pi_hat <- plogis(eta_pi)

  se <- vc <- NULL
  nu_ci <- pi_ci <- matrix(NA_real_, nrow = nrow(z), ncol = 2)
  colnames(nu_ci) <- colnames(pi_ci) <- c("lower", "upper")

  numeric_ok <- FALSE

  if (requireNamespace("numDeriv", quietly = TRUE)) {
    obj <- function(p) {
      neg_beta_loglik_combined(
        par = p,
        y = y_i,
        z = z,
        w = w,
        fem_idx = fem_idx,
        eps = eps
      )
    }

    H <- tryCatch(numDeriv::hessian(obj, fit$par), error = function(e) NULL)
    if (!is.null(H)) {
      H <- 0.5 * (H + t(H))
      vc_try <- tryCatch(solve(H), error = function(e) NULL)

      if (!is.null(vc_try) && all(is.finite(vc_try))) {
        vc <- vc_try
        se <- sqrt(pmax(diag(vc), 0))
        names(se) <- names(est)

        Va <- vc[1:pz, 1:pz, drop = FALSE]
        Vb <- vc[(pz + 1):(pz + pw), (pz + 1):(pz + pw), drop = FALSE]

        nu_ci <- delta_ci_logit(eta_nu, z, Va, conf_level = conf_level)
        pi_ci <- delta_ci_logit(eta_pi, w, Vb, conf_level = conf_level)

        numeric_ok <- (any(is.finite(nu_ci)) && any(is.finite(pi_ci)))
      }
    }
  }

  bootstrap_used <- FALSE
  boot_out <- NULL

  if (!numeric_ok && isTRUE(bootstrap_on_fail)) {
    boot_out <- bootstrap_fit_and_ci_beta(
      par_hat      = fit$par,
      y            = y_i,
      z            = z,
      w            = w,
      fem_idx      = fem_idx,
      theta_start  = fit$par,
      conf_level   = conf_level,
      B            = bootstrap_B,
      seed         = bootstrap_seed,
      eps          = eps,
      verbose      = bootstrap_verbose
    )

    nu_ci <- boot_out$nu_ci_boot
    pi_ci <- boot_out$pi_ci_boot

    vc <- boot_out$vc_boot
    se <- boot_out$se_boot
    if (!is.null(se)) names(se) <- names(est)

    bootstrap_used <- TRUE
  }

  list(
    point_estimates    = round(est, 4),
    standard_errors    = if (!is.null(se)) round(se, 4) else NULL,
    covariance_matrix  = vc,
    pi_hat             = as.numeric(pi_hat),
    nu_hat             = as.numeric(nu_hat),
    pi_ci              = pi_ci,
    nu_ci              = nu_ci,
    conf_level         = conf_level,
    ci_method          = if (numeric_ok) "delta-hessian" else if (bootstrap_used) "parametric-bootstrap" else "none",
    bootstrap_details  = boot_out
  )
}

# ==============================================================================
# 4. Prediction function for Illumina beta values
# ==============================================================================

#' Predict nu (Xa) and pi (Xi) for new samples from a fitted Illumina LVE-X Beta model
#'
#' Generates point predictions and confidence intervals for the active-allele
#' methylation level \eqn{\nu} (Xa) and inactive-allele methylation level
#' \eqn{\pi} (Xi) for new covariate values, using a fitted object returned by
#' `LVEX_fit_illumina_beta()`.
#'
#' The function prefers bootstrap-based prediction intervals if bootstrap
#' parameter draws are available; otherwise it falls back to delta-method
#' intervals using the fitted covariance matrix.
#'
#' @param fit A fitted model object returned by `LVEX_fit_illumina_beta()`.
#'   It must contain `point_estimates`, and ideally `covariance_matrix` and/or
#'   `bootstrap_details$par_boot`.
#' @param new_z Numeric matrix/data.frame. Covariates for \eqn{\nu} (Xa) on new
#'   samples. Column names must match those used in model fitting, excluding the
#'   intercept, which is added automatically.
#' @param new_w Numeric matrix/data.frame. Covariates for \eqn{\pi} (Xi) on new
#'   samples. Column names must match those used in model fitting, excluding the
#'   intercept, which is added automatically.
#' @param conf_level Optional numeric in (0,1). Confidence level for prediction
#'   intervals. If `NULL`, uses `fit$conf_level` if available, otherwise 0.95.
#'
#' @return A data frame with one row per new sample and columns:
#' \describe{
#'   \item{nu_hat}{Point prediction for \eqn{\nu}.}
#'   \item{nu_low}{Lower confidence bound for \eqn{\nu}.}
#'   \item{nu_high}{Upper confidence bound for \eqn{\nu}.}
#'   \item{pi_hat}{Point prediction for \eqn{\pi}.}
#'   \item{pi_low}{Lower confidence bound for \eqn{\pi}.}
#'   \item{pi_high}{Upper confidence bound for \eqn{\pi}.}
#'   \item{ci_method}{Method used for prediction intervals: `"bootstrap-predict"`, `"delta-hessian"`, or `"none"`.}
#' }
#'
#' @details
#' The coefficient-to-covariate mapping is inferred from `fit$point_estimates`:
#' parameters starting with `a_` are used for the \eqn{\nu} model, and those
#' starting with `b_` are used for the \eqn{\pi} model.
#'
#' Bootstrap prediction intervals are formed by transforming bootstrap parameter
#' draws through the logistic link and taking empirical quantiles. If bootstrap
#' draws are not available, the function uses delta-method intervals based on the
#' covariance matrix.
#'
#' @examples
#' \dontrun{
#' new_z <- data.frame(age = c(50, 60))
#' new_w <- data.frame(batch = c(0.1, -0.2))
#' pred <- LVEX_predict_illumina_beta(fit, new_z = new_z, new_w = new_w)
#' pred
#' }
#'
#' @export
LVEX_predict_illumina_beta <- function(fit, new_z, new_w, conf_level = NULL) {
  if (is.null(conf_level)) {
    conf_level <- if (!is.null(fit$conf_level)) fit$conf_level else 0.95
  }

  est <- fit$point_estimates
  vc  <- fit$covariance_matrix

  a_idx <- grepl("^a_", names(est))
  b_idx <- grepl("^b_", names(est))

  if (!any(a_idx)) stop("No 'a_' parameters found in fit.")
  if (!any(b_idx)) stop("No 'b_' parameters found in fit.")

  a_hat <- est[a_idx]
  b_hat <- est[b_idx]

  z_terms <- sub("^a_", "", names(a_hat))
  w_terms <- sub("^b_", "", names(b_hat))

  z_cov_names <- z_terms[-1]
  w_cov_names <- w_terms[-1]

  new_z <- as.matrix(new_z)
  new_w <- as.matrix(new_w)

  if (!all(z_cov_names %in% colnames(new_z))) {
    stop("new_z is missing covariates: ",
         paste(setdiff(z_cov_names, colnames(new_z)), collapse = ", "))
  }
  if (!all(w_cov_names %in% colnames(new_w))) {
    stop("new_w is missing covariates: ",
         paste(setdiff(w_cov_names, colnames(new_w)), collapse = ", "))
  }

  Z <- cbind(Intercept = 1, new_z[, z_cov_names, drop = FALSE])
  W <- cbind(Intercept = 1, new_w[, w_cov_names, drop = FALSE])
  colnames(Z) <- z_terms
  colnames(W) <- w_terms

  eta_nu <- drop(Z %*% a_hat)
  eta_pi <- drop(W %*% b_hat)
  nu_hat <- plogis(eta_nu)
  pi_hat <- plogis(eta_pi)

  alpha_q <- (1 - conf_level) / 2
  probs <- c(alpha_q, 1 - alpha_q)

  boot <- fit$bootstrap_details
  if (!is.null(boot) && !is.null(boot$par_boot) && nrow(boot$par_boot) >= 10) {
    par_boot <- boot$par_boot

    pz <- ncol(Z)
    pw <- ncol(W)

    nu_draws <- matrix(NA_real_, nrow(Z), nrow(par_boot))
    pi_draws <- matrix(NA_real_, nrow(W), nrow(par_boot))

    for (b in seq_len(nrow(par_boot))) {
      a_b <- par_boot[b, 1:pz]
      b_b <- par_boot[b, (pz + 1):(pz + pw)]
      nu_draws[, b] <- plogis(drop(Z %*% a_b))
      pi_draws[, b] <- plogis(drop(W %*% b_b))
    }

    nu_ci <- t(apply(nu_draws, 1, function(v) {
      quantile(v, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
    }))
    pi_ci <- t(apply(pi_draws, 1, function(v) {
      quantile(v, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
    }))
    colnames(nu_ci) <- colnames(pi_ci) <- c("lower", "upper")

    return(data.frame(
      nu_hat    = nu_hat,
      nu_low    = nu_ci[, 1],
      nu_high   = nu_ci[, 2],
      pi_hat    = pi_hat,
      pi_low    = pi_ci[, 1],
      pi_high   = pi_ci[, 2],
      ci_method = "bootstrap-predict"
    ))
  }

  nu_ci <- pi_ci <- matrix(NA_real_, nrow(Z), 2)
  colnames(nu_ci) <- colnames(pi_ci) <- c("lower", "upper")

  if (!is.null(vc)) {
    Va <- vc[a_idx, a_idx, drop = FALSE]
    Vb <- vc[b_idx, b_idx, drop = FALSE]

    nu_ci <- delta_ci_logit(eta_nu, Z, Va, conf_level = conf_level)
    pi_ci <- delta_ci_logit(eta_pi, W, Vb, conf_level = conf_level)

    return(data.frame(
      nu_hat    = nu_hat,
      nu_low    = nu_ci[, 1],
      nu_high   = nu_ci[, 2],
      pi_hat    = pi_hat,
      pi_low    = pi_ci[, 1],
      pi_high   = pi_ci[, 2],
      ci_method = "delta-hessian"
    ))
  }

  warning("No bootstrap draws and no covariance matrix available; returning NA intervals.")
  data.frame(
    nu_hat    = nu_hat,
    nu_low    = NA_real_,
    nu_high   = NA_real_,
    pi_hat    = pi_hat,
    pi_low    = NA_real_,
    pi_high   = NA_real_,
    ci_method = "none"
  )
}
