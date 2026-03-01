# ============================================================
# 0) Utility: simulate Beta-Binomial via Beta + Binomial
# ============================================================
rbbinom_beta <- function(n, size, alpha, beta) {
  p <- stats::rbeta(n, alpha, beta)
  stats::rbinom(n, size = size, prob = p)
}

# ============================================================
# 1) Female likelihood: Beta–Binomial (unchanged)
# ============================================================
bb_loglik <- function(par, r, m, z, w, eps = 1e-12) {
  a <- par[1:ncol(z)]
  b <- par[(ncol(z) + 1):(length(par) - 1)]
  kappa <- par[length(par)]  # >0 by L-BFGS-B lower bound

  nu <- plogis(drop(z %*% a))
  pi <- plogis(drop(w %*% b))

  q2  <- nu * (1 - pi)
  q1  <- nu * pi * (1 - nu) + (1 - nu) * (1 - pi)
  eta <- 0.5 * q1 + q2
  eta <- pmin(pmax(eta, eps), 1 - eps)

  alpha <- kappa * eta
  beta  <- kappa * (1 - eta)

  sum(extraDistr::dbbinom(x = r, size = m, alpha = alpha, beta = beta, log = TRUE))
}

# ============================================================
# 2) Male likelihood: Binomial (unchanged)
# ============================================================
bb_loglik_male <- function(par, r, m, z, eps = 1e-12) {
  a <- par[1:ncol(z)]
  nu <- plogis(drop(z %*% a))
  p  <- pmin(pmax(nu, eps), 1 - eps)
  sum(dbinom(x = r, size = m, prob = p, log = TRUE))
}

# ============================================================
# 3) Negative log-likelihood (unchanged; shared a across sexes)
# ============================================================
neg_bb_loglik <- function(par, r, m, z, w, fem_idx, eps = 1e-12) {
  stopifnot(is.matrix(z), is.matrix(w))
  n  <- NROW(z)
  pz <- NCOL(z)
  pw <- NCOL(w)
  stopifnot(length(r) == n, length(m) == n, NROW(w) == n, length(fem_idx) == n)

  a_f     <- par[1:pz]
  b_f     <- par[(pz + 1):(pz + pw)]
  kappa_f <- par[pz + pw + 1]

  a_m <- par[1:pz]  # shared as in your original code

  fem_idx  <- as.logical(fem_idx)
  male_idx <- !fem_idx

  ll_fem <- 0
  if (any(fem_idx)) {
    par_f <- c(a_f, b_f, kappa_f)
    ll_fem <- bb_loglik(
      par = par_f,
      r   = r[fem_idx],
      m   = m[fem_idx],
      z   = z[fem_idx, , drop = FALSE],
      w   = w[fem_idx, , drop = FALSE],
      eps = eps
    )
  }

  ll_male <- 0
  if (any(male_idx)) {
    ll_male <- bb_loglik_male(
      par = a_m,
      r   = r[male_idx],
      m   = m[male_idx],
      z   = z[male_idx, , drop = FALSE],
      eps = eps
    )
  }

  -(ll_fem + ll_male)
}

# ============================================================
# 4) Delta-method CI helper (numeric way)
# ============================================================
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

# ============================================================
# 5) Parametric bootstrap:
#    returns (i) CI for nu/pi AND (ii) parameter draws + vc/se
# ============================================================
bootstrap_fit_and_ci <- function(
    par_hat, r, m, z, w, fem_idx,
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

  pz <- ncol(z); pw <- ncol(w)
  p_all <- pz + pw + 1

  # unpack fitted params
  a_hat     <- par_hat[1:pz]
  b_hat     <- par_hat[(pz + 1):(pz + pw)]
  kappa_hat <- par_hat[pz + pw + 1]

  # fitted nu/pi/eta
  nu0 <- plogis(drop(z %*% a_hat))
  pi0 <- plogis(drop(w %*% b_hat))

  q2  <- nu0 * (1 - pi0)
  q1  <- nu0 * pi0 * (1 - nu0) + (1 - nu0) * (1 - pi0)
  eta <- 0.5 * q1 + q2
  eta <- pmin(pmax(eta, eps), 1 - eps)

  alpha <- kappa_hat * eta
  beta  <- kappa_hat * (1 - eta)

  # store bootstrap predictions & parameters
  nu_boot  <- matrix(NA_real_, nrow(z), B)
  pi_boot  <- matrix(NA_real_, nrow(w), B)
  par_boot <- matrix(NA_real_, nrow = B, ncol = p_all)

  n_ok <- 0L

  for (b in seq_len(B)) {
    r_star <- numeric(length(r))

    # females: Beta-Binomial
    if (any(fem_idx)) {
      idx <- which(fem_idx)
      r_star[idx] <- rbbinom_beta(
        n     = length(idx),
        size  = m[idx],
        alpha = alpha[idx],
        beta  = beta[idx]
      )
    }

    # males: Binomial with p = nu
    if (any(male_idx)) {
      idx <- which(male_idx)
      p_m <- pmin(pmax(nu0[idx], eps), 1 - eps)
      r_star[idx] <- rbinom(length(idx), size = m[idx], prob = p_m)
    }

    fit_b <- tryCatch(
      optim(
        par     = theta_start,
        fn      = neg_bb_loglik,
        r       = r_star,
        m       = m,
        z       = z,
        w       = w,
        fem_idx = fem_idx,
        eps     = eps,
        method  = "L-BFGS-B",
        lower   = c(rep(-Inf, pz), rep(-Inf, pw), 1e-8),
        upper   = c(rep( Inf, pz), rep( Inf, pw),  Inf),
        control = control
      ),
      error = function(e) NULL
    )

    if (is.null(fit_b) || is.null(fit_b$par) || any(!is.finite(fit_b$par))) {
      if (verbose) message("Bootstrap refit failed at b=", b)
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
    warning("Too few successful bootstrap refits (n_ok=", n_ok, "). Results may be unstable.")
  }

  # trim to successful reps
  par_boot <- par_boot[seq_len(n_ok), , drop = FALSE]
  nu_boot  <- nu_boot[,  seq_len(n_ok), drop = FALSE]
  pi_boot  <- pi_boot[,  seq_len(n_ok), drop = FALSE]

  # percentile CIs for nu/pi
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

  # bootstrap covariance / SE for parameters
  vc_boot <- if (n_ok >= 2) stats::cov(par_boot) else NULL
  se_boot <- if (!is.null(vc_boot)) sqrt(pmax(diag(vc_boot), 0)) else NULL

  list(
    nu_ci_boot = nu_ci,
    pi_ci_boot = pi_ci,
    par_boot   = par_boot,   # draws of (a,b,kappa)
    vc_boot    = vc_boot,
    se_boot    = se_boot,
    B          = B,
    n_ok       = n_ok,
    seed       = seed
  )
}

# ============================================================
# 6) Main fitting function:
#    - delta/Hessian first
#    - bootstrap fallback provides: CIs + SE + VC + draws
# ============================================================
#' Fit the LVE-X model for bisulfite sequencing methylation counts (single CpG)
#'
#' Fits the LVE-X beta-binomial model at a single CpG/site using
#' methylated (MC) and unmethylated (UC) read count matrices from bisulfite sequencing.
#' The function estimates covariate-dependent parameters for the active-allele methylation
#' level \eqn{\nu} (Xa) and the inactive-allele methylation level \eqn{\pi} (Xi), and
#' provides confidence intervals using either a Hessian-based delta method or a
#' parametric bootstrap fallback.
#'
#' @param i Integer. Column index of the CpG/site to fit (i.e., which column of `MC`/`UC`).
#' @param MC Integer matrix (n x p). Methylated read counts; rows are samples, columns are sites.
#' @param UC Integer matrix (n x p). Unmethylated read counts; same dimension/order as `MC`.
#' @param family Character. Model family; currently only `"BetaBinomial"` is supported.
#' @param z Optional numeric matrix/data.frame (n x q). Covariates for \eqn{\nu} (Xa).
#'   If `NULL`, an intercept-only model is used.
#' @param w Optional numeric matrix/data.frame (n x r). Covariates for \eqn{\pi} (Xi).
#'   If `NULL`, an intercept-only model is used.
#' @param fem_idx Integer or logical vector. Indices (or mask) indicating female samples.
#'   Must be provided.
#' @param theta_start Optional numeric vector. Starting values for optimization with length
#'   `ncol(z) + ncol(w) + 1` after adding intercepts internally. If `NULL`, defaults to
#'   zeros for regression coefficients and 10 for `kappa`.
#' @param conf_level Numeric in (0, 1). Confidence level for intervals (default 0.95).
#' @param eps Numeric. Small value used to clip probabilities to avoid 0/1 issues (default 1e-12).
#' @param bootstrap_on_fail Logical. If `TRUE`, uses parametric bootstrap to obtain CIs and
#'   covariance if Hessian-based inference fails (default `TRUE`).
#' @param bootstrap_B Integer. Number of bootstrap replicates (default 300).
#' @param bootstrap_seed Integer. Random seed for bootstrap reproducibility (default 1).
#' @param bootstrap_verbose Logical. If `TRUE`, prints bootstrap progress (default `FALSE`).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{point_estimates}{Named vector of parameter estimates (regression coefficients and `kappa`).}
#'   \item{standard_errors}{Named vector of standard errors if available; otherwise `NULL`.}
#'   \item{covariance_matrix}{Variance-covariance matrix of parameter estimates if available; otherwise `NULL`.}
#'   \item{pi_hat}{Numeric vector of fitted \eqn{\pi} for each sample.}
#'   \item{nu_hat}{Numeric vector of fitted \eqn{\nu} for each sample.}
#'   \item{pi_ci}{Matrix with columns `lower` and `upper` for \eqn{\pi}.}
#'   \item{nu_ci}{Matrix with columns `lower` and `upper` for \eqn{\nu}.}
#'   \item{conf_level}{The confidence level used.}
#'   \item{ci_method}{Character string indicating CI method: `"delta-hessian"`, `"parametric-bootstrap"`, or `"none"`.}
#'   \item{bootstrap_details}{Bootstrap output (if used), including bootstrap parameter draws for downstream prediction.}
#' }
#'
#' @details
#' Internally, the function:
#' \itemize{
#'   \item Adds an intercept column to `z` and `w` (even if you already supplied matrices).
#'   \item Optimizes the negative beta-binomial log-likelihood via `stats::optim()` with
#'   method `"L-BFGS-B"` and constraint `kappa > 0`.
#'   \item Attempts Hessian-based inference using `numDeriv::hessian()` (if installed).
#'   \item Falls back to parametric bootstrap for CIs and covariance if Hessian inversion fails.
#' }
#'
#' @examples
#' # Example (toy) usage:
#' # n samples, p sites
#' # MC/UC are integer count matrices of the same dimension
#' # fem_idx indicates which rows are females
#' # z and w are covariate matrices (rows aligned to MC/UC)
#' #
#' # NOTE: Replace neg_bb_loglik / bootstrap_fit_and_ci with your package implementations.
#' \dontrun{
#' set.seed(1)
#' n <- 20; p <- 5
#' MC <- matrix(rpois(n*p, 5), n, p)
#' UC <- matrix(rpois(n*p, 5), n, p)
#' fem_idx <- sample(c(TRUE, FALSE), n, replace = TRUE)
#' z <- data.frame(age = rnorm(n))
#' w <- data.frame(batch = rnorm(n))
#'
#' out <- LVEX_fit_bisulfite_counts(
#'   i = 2, MC = MC, UC = UC,
#'   z = z, w = w, fem_idx = fem_idx
#' )
#' str(out)
#' }
#'
#' @export
LVEX_fit_bisulfite_counts <- function(
    i, MC, UC,
    family = "BetaBinomial",
    z = NULL, w = NULL,
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
  if (family != "BetaBinomial") stop("Error: 'family' type is not correct.")

  z_name <- if (!is.null(z)) colnames(z) else NULL
  w_name <- if (!is.null(w)) colnames(w) else NULL

  n_rows <- nrow(MC)

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

  pz <- ncol(z); pw <- ncol(w)

  if (is.null(theta_start)) {
    theta_start <- c(rep(0, pz), rep(0, pw), 10)
  }

  r_i <- MC[, i]
  m_i <- MC[, i] + UC[, i]

  fit <- optim(
    par     = theta_start,
    fn      = neg_bb_loglik,
    r       = r_i,
    m       = m_i,
    z       = z,
    w       = w,
    fem_idx = fem_idx,
    eps     = eps,
    method  = "L-BFGS-B",
    lower   = c(rep(-Inf, pz), rep(-Inf, pw), 1e-8),
    upper   = c(rep( Inf, pz), rep( Inf, pw),  Inf),
    control = list(maxit = 5000, factr = 1e7)
  )

  a_names <- if (!is.null(z_name) && any(nzchar(z_name)))
    paste0("a_", c("Intercept", z_name)) else paste0("a", seq_len(pz) - 1L)
  b_names <- if (!is.null(w_name) && any(nzchar(w_name)))
    paste0("b_", c("Intercept", w_name)) else paste0("b", seq_len(pw) - 1L)

  est <- setNames(fit$par, c(a_names, b_names, "kappa"))

  # point predictions
  a_hat <- fit$par[1:pz]
  b_hat <- fit$par[(pz + 1):(pz + pw)]
  eta_nu <- drop(z %*% a_hat)
  eta_pi <- drop(w %*% b_hat)
  nu_hat <- plogis(eta_nu)
  pi_hat <- plogis(eta_pi)

  # --- numeric SE + delta-method CI ---
  se <- vc <- NULL
  nu_ci <- pi_ci <- matrix(NA_real_, nrow = nrow(z), ncol = 2)
  colnames(nu_ci) <- colnames(pi_ci) <- c("lower", "upper")

  numeric_ok <- FALSE

  if (requireNamespace("numDeriv", quietly = TRUE)) {
    obj <- function(p) neg_bb_loglik(p, r = r_i, m = m_i, z = z, w = w, fem_idx = fem_idx, eps = eps)
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

  # --- bootstrap fallback ---
  bootstrap_used <- FALSE
  boot_out <- NULL

  if (!numeric_ok && isTRUE(bootstrap_on_fail)) {
    boot_out <- bootstrap_fit_and_ci(
      par_hat = fit$par,
      r       = r_i,
      m       = m_i,
      z       = z,
      w       = w,
      fem_idx = fem_idx,
      theta_start = fit$par,   # warm-start
      conf_level = conf_level,
      B = bootstrap_B,
      seed = bootstrap_seed,
      eps = eps,
      verbose = bootstrap_verbose
    )

    # Use bootstrap CIs for nu/pi
    nu_ci <- boot_out$nu_ci_boot
    pi_ci <- boot_out$pi_ci_boot

    # And ALSO use bootstrap covariance + SE so prediction can work later
    vc <- boot_out$vc_boot
    se <- boot_out$se_boot
    if (!is.null(se)) names(se) <- names(est)

    bootstrap_used <- TRUE
  }

  list(
    point_estimates    = round(est, 4),
    standard_errors    = if (!is.null(se)) round(se, 4) else NULL,
    covariance_matrix  = vc,                      # now filled from bootstrap if needed
    pi_hat             = as.numeric(pi_hat),
    nu_hat             = as.numeric(nu_hat),
    pi_ci              = pi_ci,
    nu_ci              = nu_ci,
    conf_level         = conf_level,
    ci_method          = if (numeric_ok) "delta-hessian" else if (bootstrap_used) "parametric-bootstrap" else "none",
    bootstrap_details  = boot_out                 # includes par_boot for future prediction
  )
}
# ============================================================
# 7) Predict function:
#    - if bootstrap draws exist -> bootstrap prediction intervals
#    - else if covariance exists -> delta-method intervals
# ============================================================
#' Predict nu (Xa) and pi (Xi) for new samples from a fitted bisulfite-counts LVE-X model
#'
#' Generates point predictions and confidence intervals for the active-allele methylation
#' level \eqn{\nu} (Xa) and the inactive-allele methylation level \eqn{\pi} (Xi) for new
#' covariate values, using a model fitted by `LVEX_fit_bisulfite_counts()`.
#'
#' The function prefers bootstrap-based prediction intervals if the fitted object contains
#' bootstrap parameter draws; otherwise it falls back to delta-method intervals using the
#' fitted covariance matrix.
#'
#' @param fit A fitted model object returned by `LVEX_fit_bisulfite_counts()`.
#'   It must contain `point_estimates`, and ideally `covariance_matrix` and/or
#'   `bootstrap_details$par_boot`.
#' @param new_z Numeric matrix/data.frame. Covariates for \eqn{\nu} (Xa) on new samples,
#'   with column names matching the covariates used in the fitting function (excluding
#'   the intercept, which is added automatically).
#' @param new_w Numeric matrix/data.frame. Covariates for \eqn{\pi} (Xi) on new samples,
#'   with column names matching the covariates used in the fitting function (excluding
#'   the intercept, which is added automatically).
#' @param conf_level Optional numeric in (0, 1). Confidence level for intervals.
#'   If `NULL`, the function uses `fit$conf_level` if available, otherwise defaults to 0.95.
#'
#' @return A data frame with one row per new sample and columns:
#' \describe{
#'   \item{nu_hat}{Point prediction for \eqn{\nu}.}
#'   \item{nu_low}{Lower confidence bound for \eqn{\nu}.}
#'   \item{nu_high}{Upper confidence bound for \eqn{\nu}.}
#'   \item{pi_hat}{Point prediction for \eqn{\pi}.}
#'   \item{pi_low}{Lower confidence bound for \eqn{\pi}.}
#'   \item{pi_high}{Upper confidence bound for \eqn{\pi}.}
#'   \item{ci_method}{Which method was used: `"bootstrap-predict"`, `"delta-hessian"`, or `"none"`.}
#' }
#'
#' @details
#' The mapping between coefficients and covariates is inferred from the names in
#' `fit$point_estimates`:
#' \itemize{
#'   \item Parameters starting with `a_` are used for the \eqn{\nu} model.
#'   \item Parameters starting with `b_` are used for the \eqn{\pi} model.
#' }
#' Therefore, the column names of `new_z` and `new_w` must contain the same covariate
#' names used in fitting (excluding the intercept).
#'
#' Bootstrap prediction intervals are computed by transforming each bootstrap draw of
#' parameters through the logistic link and taking empirical quantiles.
#' Delta-method intervals use `delta_ci_logit()` with the appropriate sub-block of the
#' covariance matrix.
#'
#' @examples
#' \dontrun{
#' # Suppose `fit` is returned by LVEX_fit_bisulfite_counts()
#' new_z <- data.frame(age = c(50, 60))
#' new_w <- data.frame(batch = c(0.1, -0.2))
#' pred <- LVEX_predict_bisulfite_counts(fit, new_z = new_z, new_w = new_w)
#' pred
#' }
#'
#' @export
LVEX_predict_bisulfite_counts <- function(fit, new_z, new_w, conf_level = NULL) {
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

  # point predictions
  eta_nu <- drop(Z %*% a_hat)
  eta_pi <- drop(W %*% b_hat)
  nu_hat <- plogis(eta_nu)
  pi_hat <- plogis(eta_pi)

  alpha_q <- (1 - conf_level) / 2
  probs <- c(alpha_q, 1 - alpha_q)

  # ---- Preferred: bootstrap prediction intervals if available ----
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

    nu_ci <- t(apply(nu_draws, 1, function(v) quantile(v, probs = probs, na.rm = TRUE, names = FALSE, type = 8)))
    pi_ci <- t(apply(pi_draws, 1, function(v) quantile(v, probs = probs, na.rm = TRUE, names = FALSE, type = 8)))
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

  # ---- Fallback: delta-method using covariance matrix ----
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
