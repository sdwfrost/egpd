## Spliced discrete forecasting model:
## additive quantile regression for the bulk + DEGPD tail regression

.spliced_require_qgam <- function() {
  if (!requireNamespace("qgam", quietly = TRUE)) {
    stop(
      "fit_spliced_degpd() requires the 'qgam' package. Install it with install.packages('qgam').",
      call. = FALSE
    )
  }
}

.spliced_degpd_formula <- function(formula, response = NULL) {
  rhs <- formula[[length(formula)]]
  f <- if (is.null(response)) call("~", rhs) else call("~", as.name(response), rhs)
  stats::as.formula(f, env = environment(formula))
}

.spliced_degpd_intercept_only <- function(formula) {
  tt <- terms.formula(formula)
  attr(tt, "intercept") == 1L && !length(attr(tt, "term.labels")) && !length(attr(tt, "offset"))
}

.normalize_spliced_bulk_probs <- function(bulk_probs, threshold_prob, fit_tail) {
  if (!is.numeric(bulk_probs) || length(bulk_probs) == 0L || any(!is.finite(bulk_probs))) {
    stop("'bulk_probs' must be a non-empty numeric vector", call. = FALSE)
  }
  bulk_probs <- sort(unique(as.numeric(bulk_probs)))
  if (any(bulk_probs <= 0 | bulk_probs >= 1)) {
    stop("'bulk_probs' must lie strictly between 0 and 1", call. = FALSE)
  }

  if (!fit_tail) return(bulk_probs)
  if (!is.numeric(threshold_prob) || length(threshold_prob) != 1L || !is.finite(threshold_prob) ||
      threshold_prob <= 0 || threshold_prob >= 1)
    stop("'threshold_prob' must be a single probability in (0, 1)", call. = FALSE)
  c(bulk_probs[bulk_probs < threshold_prob], threshold_prob)
}

.normalize_spliced_tail_formula <- function(tail_formula, response = "exceedance") {
  valid_names <- c("lsigma", "lxi", "lkappa")

  if (is.null(tail_formula)) {
    tail_formula <- list(lsigma = ~1, lxi = ~1, lkappa = ~1)
  } else if (inherits(tail_formula, "formula")) {
    tail_formula <- list(lsigma = tail_formula)
  } else if (!is.list(tail_formula)) {
    stop("'tail_formula' must be NULL, a formula, or a named list of formulae", call. = FALSE)
  }

  if (is.null(names(tail_formula))) {
    if (length(tail_formula) == 1L) {
      names(tail_formula) <- "lsigma"
    } else if (length(tail_formula) == length(valid_names)) {
      names(tail_formula) <- valid_names
    } else {
      stop(
        "'tail_formula' must have names from lsigma, xi/lxi, lkappa, or length 1/3 if unnamed",
        call. = FALSE
      )
    }
  }

  if ("xi" %in% names(tail_formula)) {
    if ("lxi" %in% names(tail_formula)) {
      stop("Specify at most one of 'xi' and 'lxi' in 'tail_formula'", call. = FALSE)
    }
    names(tail_formula)[names(tail_formula) == "xi"] <- "lxi"
  }

  bad_names <- setdiff(names(tail_formula), valid_names)
  if (length(bad_names) > 0L) {
    stop("Unknown tail formula component(s): ", paste(bad_names, collapse = ", "), call. = FALSE)
  }

  for (nm in setdiff(valid_names, names(tail_formula))) {
    tail_formula[[nm]] <- ~1
  }
  tail_formula <- tail_formula[valid_names]

  if (!all(vapply(tail_formula, inherits, what = "formula", logical(1)))) {
    stop("All 'tail_formula' elements must be formulae", call. = FALSE)
  }

  tail_formula$lsigma <- .spliced_degpd_formula(tail_formula$lsigma, response)
  tail_formula$lxi <- .spliced_degpd_formula(tail_formula$lxi)
  tail_formula$lkappa <- .spliced_degpd_formula(tail_formula$lkappa)
  tail_formula
}

.spliced_degpd_has_specials <- function(formula) {
  terms_obj <- terms.formula(formula, specials = c("s", "te", "ti"))
  specials <- attr(terms_obj, "specials")
  if (is.null(specials)) {
    return(FALSE)
  }
  if (length(specials) == 0L) {
    return(FALSE)
  }
  any(unlist(specials, use.names = FALSE) > 0L)
}

.spliced_control_names <- function(x, allowed, label) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyDuplicated(names(x)) ||
                                     any(!names(x) %in% allowed))))
    stop("'", label, "' must be a list with unique names from: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)
}

.fit_spliced_degpd_tail <- function(tail_formula, tail_data, fit_control = list(),
                                    training_data = tail_data) {
  .spliced_control_names(fit_control, c("start", "fix.arg", "optim.method", "control"),
                        "egpd_control")
  if (.spliced_degpd_has_specials(tail_formula$lsigma)) {
    stop("tail_formula$lsigma currently supports parametric terms only", call. = FALSE)
  }
  if (!.spliced_degpd_intercept_only(tail_formula$lxi)) {
    stop("tail_formula$lxi currently supports only ~1", call. = FALSE)
  }
  if (!.spliced_degpd_intercept_only(tail_formula$lkappa)) {
    stop("tail_formula$lkappa currently supports only ~1", call. = FALSE)
  }

  cc <- .gam_complete_rows(list(tail_formula$lsigma), tail_data,
                           required = as.character(tail_formula$lsigma[[2L]]))
  n_used <- sum(cc)
  if (n_used < nrow(tail_data)) {
    warning(sprintf("Omitted %d incomplete tail observation(s); %d exceedances remain",
                    nrow(tail_data) - n_used, n_used), call. = FALSE)
  }
  if (n_used < 10L) {
    stop("Not enough complete exceedances were available to fit the tail model", call. = FALSE)
  }
  mf <- stats::model.frame(tail_formula$lsigma, data = tail_data[cc, , drop = FALSE],
                           na.action = stats::na.fail)
  ## model.frame has recorded the training bases in the terms' predvars.
  tt <- stats::terms(mf)
  scale_terms <- stats::delete.response(tt)
  ## Check rank in the requested factor coding, not one with unsupported
  ## levels silently removed. Numeric bases are still fitted on exceedances.
  for (nm in names(stats::.getXlevels(tt, mf))) {
    i <- match(nm, names(mf))
    v <- eval(attr(tt, "predvars")[[i + 1L]], training_data, environment(tt))
    if (is.character(v)) v <- factor(v)
    if (any(!levels(v) %in% v)) v <- droplevels(v)
    idx <- match(as.character(mf[[nm]]), as.character(v))
    if (anyNA(idx))
      stop("Inconsistent levels for tail predictor '", nm,
           "' on the bulk and exceedance rows; supply it as a column.", call. = FALSE)
    mf[[nm]] <- v[idx]
  }
  scale_xlevels <- stats::.getXlevels(tt, mf)
  X <- stats::model.matrix(tt, data = mf)
  if (ncol(X)) {
    decomp <- qr(X, tol = sqrt(.Machine$double.eps))
    if (decomp$rank < ncol(X)) {
      aliases <- colnames(X)[decomp$pivot[seq.int(decomp$rank + 1L, ncol(X))]]
      stop("Tail scale design is rank deficient on the ", n_used,
           " complete exceedances; coefficients are not identifiable for: ",
           paste(aliases, collapse = ", "), ". Simplify 'tail_formula'.", call. = FALSE)
    }
  }
  scale_contrasts <- attr(X, "contrasts")
  offset <- stats::model.offset(mf)
  if (is.null(offset)) offset <- 0
  y <- stats::model.response(mf)

  fix.arg <- fit_control$fix.arg
  if (is.null(fix.arg)) {
    fix.arg <- list(kappa = 1)
  }
  .spliced_control_names(fix.arg, c("xi", "kappa"), "egpd_control$fix.arg")
  spec <- .distfit_parspec(1L, "degpd")
  .distfit_check_parlist(fix.arg, "fix.arg", spec, "degpd")

  init_fit <- tryCatch(
    distfit(y, type = 1, family = "degpd", fix.arg = fix.arg, hessian = FALSE, optim.method = "BFGS"),
    error = function(e) NULL
  )

  sigma0 <- if (!is.null(init_fit) && "sigma" %in% names(init_fit$estimate)) {
    unname(init_fit$estimate["sigma"])
  } else {
    max(mean(y + 0.5), 1e-2)
  }
  xi0 <- if (!is.null(init_fit) && "xi" %in% names(init_fit$estimate)) {
    unname(init_fit$estimate["xi"])
  } else if (!is.null(fix.arg$xi)) {
    fix.arg$xi
  } else {
    0.2
  }
  kappa0 <- if (!is.null(init_fit) && "kappa" %in% names(init_fit$estimate)) {
    unname(init_fit$estimate["kappa"])
  } else if (!is.null(fix.arg$kappa)) {
    fix.arg$kappa
  } else {
    1
  }

  beta0 <- rep(0, ncol(X))
  if (attr(tt, "intercept") == 1L) beta0[1] <- log(sigma0) - mean(offset)

  start <- fit_control$start
  if (!is.null(start)) {
    .spliced_control_names(start, c("beta", "xi", "kappa"), "egpd_control$start")
    .distfit_check_parlist(start, "start", spec, "degpd")
    if (!is.null(start$beta)) {
      if (!is.numeric(start$beta) || length(start$beta) != ncol(X) || any(!is.finite(start$beta))) {
        stop("tail optimizer 'start$beta' must be finite coefficients of the log-scale design",
             call. = FALSE)
      }
      beta0 <- as.numeric(start$beta)
    }
    if (!is.null(start$xi)) {
      xi0 <- as.numeric(start$xi)
    }
    if (!is.null(start$kappa)) {
      kappa0 <- as.numeric(start$kappa)
    }
  }

  ## Tail shape xi is parameterised on the unconstrained identity scale (xi in R).
  theta0 <- beta0
  if (is.null(fix.arg$xi)) {
    theta0 <- c(theta0, xi0)
  }
  if (is.null(fix.arg$kappa)) {
    theta0 <- c(theta0, log(kappa0))
  }

  unpack <- function(theta) {
    beta <- theta[seq_len(ncol(X))]
    pos <- ncol(X)
    xi <- if (is.null(fix.arg$xi)) {
      pos <- pos + 1L
      theta[pos]
    } else {
      fix.arg$xi
    }
    kappa <- if (is.null(fix.arg$kappa)) {
      pos <- pos + 1L
      exp(theta[pos])
    } else {
      fix.arg$kappa
    }
    list(beta = beta, xi = xi, kappa = kappa)
  }

  nll <- function(theta) {
    pars <- unpack(theta)
    if (!all(is.finite(pars$beta)) || !is.finite(pars$xi) ||
        !is.finite(pars$kappa) || pars$kappa <= 0) {
      return(1e20)
    }
    eta <- drop(X %*% pars$beta) + offset
    if (any(!is.finite(eta)) || any(eta > 700)) {
      return(1e20)
    }
    sigma <- exp(eta)
    if (any(!is.finite(sigma)) || any(sigma <= 0)) {
      return(1e20)
    }
    if (pars$xi < 0 && any(1 + pars$xi * y / sigma <= 0)) {
      return(1e20)
    }
    ldens <- ddiscegpd(y, sigma = sigma, xi = pars$xi, kappa = pars$kappa,
                       type = 1, log = TRUE)
    if (any(!is.finite(ldens))) {
      return(1e20)
    }
    -sum(ldens)
  }

  if (nll(theta0) >= 1e20 && is.null(fix.arg$xi) && (is.null(start) || is.null(start$xi))) {
    theta0[ncol(X) + 1L] <- 0.1
  }

  optim_method <- fit_control$optim.method
  if (is.null(optim_method)) {
    optim_method <- "BFGS"
  }
  optim_control <- fit_control$control
  if (is.null(optim_control)) {
    optim_control <- list()
  }
  default_ctrl <- list(maxit = 500L)
  if (optim_method %in% c("BFGS", "L-BFGS-B") && is.null(optim_control$ndeps)) {
    default_ctrl$ndeps <- rep(1e-5, length(theta0))
  }

  opt <- optim(
    par = theta0,
    fn = nll,
    method = optim_method,
    control = modifyList(default_ctrl, optim_control),
    hessian = FALSE
  )
  if ((opt$value >= 1e20 || opt$convergence != 0L) && length(theta0) >= 2L &&
      identical(optim_method, "BFGS") && is.null(fit_control$optim.method)) {
    opt_nm <- tryCatch(
      optim(
        par = if (opt$value < nll(theta0)) opt$par else theta0,
        fn = nll,
        method = "Nelder-Mead",
        control = modifyList(list(maxit = 2000L), optim_control),
        hessian = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(opt_nm) && opt_nm$value < opt$value) {
      opt <- opt_nm
      optim_method <- "Nelder-Mead"
    }
  }
  ## nll() is 1e20 where some exceedance has zero likelihood, so an optimum there is no fit
  ## (as in distfit()); a fit that did not converge is reported as egpd() reports one.
  if (!is.finite(opt$value) || opt$value >= 1e20) {
    stop("the DEGPD tail fit found no parameter value with a positive likelihood ",
         "from its starting values (an exceedance lies outside the support of every ",
         "distribution visited); supply a feasible `egpd_control$start`.", call. = FALSE)
  }
  converged <- identical(as.integer(opt$convergence), 0L)
  if (!converged) {
    msg <- sprintf(paste0("the DEGPD tail fit did not converge (optim %s, code %d%s): its ",
                          "log-likelihood is not reported and the estimates are unreliable. ",
                          "Try `egpd_control$start`, `egpd_control$optim.method` or a larger ",
                          "`egpd_control$control$maxit`."),
                   optim_method,
                   as.integer(opt$convergence),
                   if (!is.null(opt$message)) paste0(": ", opt$message) else "")
    warning(structure(class = c("egpd_nonconvergence", "warning", "condition"),
                      list(message = msg, call = NULL)))
  }
  pars <- unpack(opt$par)

  structure(
    list(
      formula = stats::formula(scale_terms),
      terms = scale_terms,
      xlevels = scale_xlevels,
      contrasts = scale_contrasts,
      coefficients = setNames(pars$beta, colnames(X)),
      xi = pars$xi,
      kappa = pars$kappa,
      fix.arg = fix.arg,
      nobs = n_used,
      convergence = opt$convergence,
      loglik = if (converged) -opt$value else NA_real_,
      optim = opt
    ),
    class = "spliced_tail_degpd"
  )
}

## mqgam omits raw NAs before evaluating its formula. Bind complete terms
## that resolve a raw NA to private columns so the bulk uses our joint sample.
.spliced_bulk_terms <- function(formula, data) {
  raw <- intersect(.call_vars(formula[[3L]]), names(data))
  if (!anyNA(data[raw])) return(list(formula=formula, data=data))
  ff <- mgcv::interpret.gam(formula)$fake.formula
  environment(ff) <- environment(formula)
  tt <- stats::delete.response(stats::terms(ff))
  vars <- as.list(attr(tt,"variables"))[-1L]
  idx <- which(vapply(vars, function(v) anyNA(data[intersect(.call_vars(v),raw)]), TRUE))
  if (!length(idx)) return(list(formula=formula, data=data))
  mf <- stats::model.frame(tt, data, na.action=stats::na.fail)
  pred <- as.list(attr(stats::terms(mf),"predvars"))[-1L][idx]
  repl <- vector("list",length(idx))
  used <- c(names(data),all.vars(formula))
  for (j in seq_along(idx)) {
    nm <- .aux_response_name(".egpd_bulk", used)
    used <- c(used,nm)
    names(pred)[j] <- nm
    data[[nm]] <- mf[[idx[j]]]
    repl[[j]] <- if (idx[j] %in% attr(tt,"offset")) call("offset", as.name(nm)) else as.name(nm)
  }
  replace_term <- function(x) {
    j <- which(vapply(vars[idx], identical, TRUE, y=x))
    if (length(j)) return(repl[[j[1L]]])
    if (is.call(x)) for (i in seq_along(x)[-1L]) x[[i]] <- replace_term(x[[i]])
    x
  }
  formula[[3L]] <- replace_term(formula[[3L]])
  list(formula=formula, data=data, pred=list(calls=pred, env=environment(formula)))
}

.spliced_bulk_data <- function(bulk_fit, data) {
  data <- .gam_offset_data(bulk_fit[["egpd.offset.formula"]], data)
  pred <- bulk_fit[["egpd.bulk.terms"]]
  if (length(pred) && nrow(data))
    data[names(pred$calls)] <- lapply(pred$calls, eval, envir=data, enclos=pred$env)
  data
}

.spliced_degpd_qdo_matrix <- function(bulk_fit, probs, newdata) {
  if (!nrow(newdata)) return(matrix(numeric(), 0L, length(probs)))
  pred <- qgam::qdo(obj = bulk_fit, qu = probs, fun = .gam_predict_rows,
                    newdata = newdata, type = "link")
  if (is.list(pred)) {
    mat <- do.call(cbind, lapply(pred, function(x) as.numeric(x)))
  } else {
    mat <- as.matrix(pred)
  }
  if (nrow(mat) != nrow(newdata) && ncol(mat) == nrow(newdata)) {
    mat <- t(mat)
  }
  if (nrow(mat) != nrow(newdata) || ncol(mat) != length(probs)) {
    stop("Bulk quantile predictions returned an unexpected shape", call. = FALSE)
  }
  storage.mode(mat) <- "double"
  colnames(mat) <- as.character(probs)
  mat
}

.spliced_degpd_sorted_matrix <- function(mat) {
  if (!nrow(mat) || ncol(mat) == 1L) return(mat)
  out <- t(apply(mat, 1L, sort, na.last = TRUE))
  dimnames(out) <- dimnames(mat)
  out
}

.spliced_degpd_prepare_newdata <- function(object, newdata) {
  if (!missing(newdata)) {
    return(newdata)
  }
  if (is.null(object$data)) {
    stop("Supply 'newdata' or fit with keep_data = TRUE.", call. = FALSE)
  }
  object$data
}

.spliced_degpd_bulk_quantile <- function(probs, values, xout) {
  probs <- as.numeric(probs)
  values <- pmax(as.numeric(values), 0)
  xout <- as.numeric(xout)

  stats::approx(
    x = c(0, probs, 1),
    y = c(0, values, tail(values, 1)),
    xout = xout,
    yleft = 0,
    yright = values[length(values)],
    ties = "ordered"
  )$y
}

.spliced_degpd_bulk_cdf <- function(values, probs, at, upper_prob) {
  values <- pmax(as.numeric(values), 0)
  probs <- as.numeric(probs)
  at <- as.numeric(at)

  values <- c(0, values)
  probs <- c(0, probs)
  knots <- unique(values)
  left <- probs[match(knots, values)]
  right <- probs[length(probs) + 1L - match(knots, rev(values))]
  right[length(right)] <- upper_prob
  interval <- findInterval(at, knots)
  out <- numeric(length(at))
  out[interval == length(knots)] <- upper_prob
  mid <- which(interval > 0 & interval < length(knots))
  j <- interval[mid]
  out[mid] <- right[j] + (left[j + 1L] - right[j]) *
    (at[mid] - knots[j]) / (knots[j + 1L] - knots[j])
  out
}

.spliced_degpd_tail_parameters <- function(object, newdata) {
  if (!nrow(newdata)) return(data.frame(scale = numeric(), shape = numeric(), kappa = numeric()))
  tt <- object$tail_fit$terms
  if (is.null(tt)) tt <- object$tail_fit$formula
  mf <- stats::model.frame(tt, data = newdata, na.action = stats::na.pass,
                           xlev = object$tail_fit$xlevels)
  X <- stats::model.matrix(tt, data = mf, contrasts.arg = object$tail_fit$contrasts)
  offset <- stats::model.offset(mf)
  if (is.null(offset)) offset <- 0
  data.frame(
    scale = exp(drop(X %*% object$tail_fit$coefficients) + offset),
    shape = rep(object$tail_fit$xi, nrow(newdata)),
    kappa = rep(object$tail_fit$kappa, nrow(newdata))
  )
}

## First count modelled by the tail.  Maia et al. (2026, Eq. 3) assign
## y <= floor(q) to the bulk, with alpha_tilde = F_QR(floor(q)), and the rest
## to the DGP.  Using floor(q) + 1 (rather than ceiling(q), which equals q for
## an integer threshold quantile) keeps the fitted exceedances Y - s | Y >= s
## and the predictive tail mass 1 - alpha_tilde = P(Y >= s) consistent.
.spliced_degpd_shift <- function(threshold_values) {
  pmax(floor(as.numeric(threshold_values)) + 1, 0)
}

.spliced_degpd_alpha_tilde <- function(bulk_mat, bulk_probs, threshold_values,
                                       threshold_prob, eps = 1e-12) {
  n <- nrow(bulk_mat)
  out <- numeric(n)
  for (i in seq_len(n)) {
    val <- .spliced_degpd_bulk_cdf(
      values = bulk_mat[i, ],
      probs = bulk_probs,
      at = floor(threshold_values[i]),
      upper_prob = threshold_prob
    )
    out[i] <- min(max(val, 0), 1 - eps)
  }
  out
}

.spliced_degpd_tail_cdf <- function(x, shift, scale, shape, kappa) {
  out <- numeric(length(x))
  idx <- x >= shift
  if (any(idx)) {
    out[idx] <- pdiscegpd(
      q = x[idx] - shift,
      sigma = scale,
      xi = shape,
      kappa = kappa,
      type = 1
    )
  }
  out
}

.spliced_degpd_row_cdf <- function(at, bulk_row, bulk_probs, upper_prob,
                                   has_tail = FALSE, shift = 0,
                                   alpha_tilde = 0, scale = 1,
                                   shape = 0, kappa = 1) {
  at <- as.numeric(at)
  out <- numeric(length(at))
  neg <- at < 0
  if (all(neg)) return(out)
  ok <- !neg
  at_ok <- floor(at[ok])
  if (!has_tail) {
    res <- .spliced_degpd_bulk_cdf(bulk_row, bulk_probs, at_ok, upper_prob = 1)
  } else {
    res <- numeric(length(at_ok))
    below <- at_ok < shift
    if (any(below)) {
      res[below] <- pmin(
        .spliced_degpd_bulk_cdf(bulk_row, bulk_probs, at_ok[below], upper_prob = upper_prob),
        alpha_tilde
      )
    }
    above <- !below
    if (any(above)) {
      res[above] <- alpha_tilde + (1 - alpha_tilde) * pdiscegpd(
        q = at_ok[above] - shift,
        sigma = scale,
        xi = shape,
        kappa = kappa,
        type = 1
      )
    }
  }
  res[res < 0] <- 0
  res[res > 1] <- 1
  out[ok] <- res
  out
}

.spliced_degpd_bisect_count <- function(lo, hi, reached) {
  ## lo fails and hi reaches the target. Above 2^53 even mid + 1 may be mid.
  repeat {
    mid <- floor(lo + (hi - lo) / 2)
    if (mid <= lo || mid >= hi) return(hi)
    if (reached(mid)) hi <- mid else lo <- mid
  }
}

.spliced_degpd_row_quantile <- function(probs, bulk_row, bulk_probs, upper_prob,
                                        has_tail = FALSE, shift = 0,
                                        alpha_tilde = 0, scale = 1,
                                        shape = 0, kappa = 1) {
  cdf_at <- function(k) {
    .spliced_degpd_row_cdf(
      at = k, bulk_row = bulk_row, bulk_probs = bulk_probs,
      upper_prob = upper_prob, has_tail = has_tail, shift = shift,
      alpha_tilde = alpha_tilde, scale = scale, shape = shape, kappa = kappa
    )
  }

  vapply(probs, function(p) {
    if (is.nan(p)) return(NaN)
    if (is.na(p)) return(NA_real_)
    if (p < 0 || p > 1) return(NaN)
    if (p == 0) return(0)

    if (!has_tail || (p <= alpha_tilde && shift > 0)) {
      hi <- if (has_tail) shift else ceiling(max(bulk_row, 0))
      return(.spliced_degpd_bisect_count(-1, hi, function(k) cdf_at(k) >= p))
    }

    if (is.na(scale)) return(NA_real_)
    if (p == 1) {
      return(if (shape < 0) shift + .gpd_last_count(scale, shape) else Inf)
    }

    ## Invert the smaller conditional probability directly. In particular,
    ## never obtain a small survivor by subtracting a rounded tail CDF from 1.
    lower <- p - alpha_tilde < 1 - p
    target <- (if (lower) p - alpha_tilde else 1 - p) / (1 - alpha_tilde)
    exc <- qdiscegpd(target, sigma = scale, xi = shape, kappa = kappa, lower.tail = lower)
    k <- shift + exc
    if (!is.finite(k)) return(k)

    ## Honor isolated, ordinary boundaries on the same scale as the public
    ## CDF (including p exactly F(k)). Near 1 the lower CDF no longer resolves
    ## the survivor accurately; do not use its rounded values to move counts.
    if (k < 2^53 && 1 - p > sqrt(.Machine$double.eps)) {
      near <- k + (-2:2)
      f <- cdf_at(near)
      hit <- which(f[1:3] < p & p <= f[2:4] & f[2:4] < f[3:5])
      if (length(hit)) return(near[hit[1L] + 1L])
    }

    reached <- function(x) {
      v <- pdiscegpd(x, sigma = scale, xi = shape, kappa = kappa, lower.tail = lower)
      if (lower) v >= target else v <= target
    }
    ## Bracket on the exceedance scale before adding the shift. A step of
    ## at least one ulp still moves a double where unit increments cannot.
    step <- max(1, ceiling(exc * .Machine$double.eps))
    if (reached(exc)) {
      hi <- exc
      repeat {
        lo <- max(-1, hi - step)
        if (lo < 0 || !reached(lo)) break
        hi <- lo
        step <- 2 * step
      }
    } else {
      lo <- exc
      repeat {
        hi <- min(lo + step, .Machine$double.xmax)
        if (reached(hi)) break
        if (hi == .Machine$double.xmax) return(Inf)
        lo <- hi
        step <- 2 * step
      }
    }
    shift + .spliced_degpd_bisect_count(lo, hi, reached)
  }, numeric(1))
}

.spliced_degpd_cdf_matrix <- function(object, newdata, at, bulk_mat = NULL, tail_pars = NULL) {
  at <- as.numeric(at)
  if (!length(at) || any(!is.finite(at)) || any(at != floor(at))) {
    stop("'at' must contain finite integer values", call. = FALSE)
  }
  at <- sort(unique(at))

  if (is.null(bulk_mat)) {
    bulk_mat <- .spliced_degpd_sorted_matrix(
      .spliced_degpd_qdo_matrix(object$bulk_fit, object$bulk_probs,
                                .spliced_bulk_data(object$bulk_fit, newdata))
    )
  }

  has_tail <- !is.null(object$tail_fit)
  upper_index <- if (!has_tail) length(object$bulk_probs) else object$threshold_index
  upper_prob <- if (!has_tail) 1 else object$threshold_prob
  threshold_idx <- object$threshold_index
  threshold_values <- unname(bulk_mat[, threshold_idx])
  threshold_shift <- .spliced_degpd_shift(threshold_values)
  bulk_probs <- object$bulk_probs[seq_len(upper_index)]
  bulk_sub <- bulk_mat[, seq_len(upper_index), drop = FALSE]

  alpha_tilde <- if (has_tail) {
    .spliced_degpd_alpha_tilde(bulk_sub, bulk_probs, threshold_values, upper_prob)
  } else {
    rep(1, nrow(bulk_sub))
  }

  if (has_tail && is.null(tail_pars)) {
    tail_pars <- .spliced_degpd_tail_parameters(object, newdata)
  }

  cdf <- matrix(0, nrow = nrow(bulk_sub), ncol = length(at))
  for (i in seq_len(nrow(bulk_sub))) {
    cdf[i, ] <- .spliced_degpd_row_cdf(
      at = at,
      bulk_row = bulk_sub[i, ],
      bulk_probs = bulk_probs,
      upper_prob = upper_prob,
      has_tail = has_tail,
      shift = threshold_shift[i],
      alpha_tilde = alpha_tilde[i],
      scale = if (has_tail) tail_pars$scale[i] else 1,
      shape = if (has_tail) tail_pars$shape[i] else 0,
      kappa = if (has_tail) tail_pars$kappa[i] else 1
    )
  }

  if (nrow(cdf) && ncol(cdf) > 1L) {
    cdf <- matrix(t(apply(cdf, 1, cummax)), nrow = nrow(cdf), ncol = ncol(cdf))
  }
  colnames(cdf) <- paste0("F:", format(at, scientific = FALSE, trim = TRUE))
  rownames(cdf) <- rownames(bulk_mat)
  cdf
}

#' Fit a spliced discrete EGPD forecasting model
#'
#' Fits a hybrid count-distribution model following the X-flexForecast
#' framework of Maia, Castro-Camilo, and Browell (2026): multiple additive
#' quantile regressions describe the bulk of the count distribution, while
#' threshold exceedances are modeled with a type-1 DEGPD tail.
#'
#' @details Bulk quantiles interpolate from probability zero at value zero,
#' with negative fitted knots set to zero for bulk interpolation. Tied
#' knots produce atoms, with the CDF taking the upper probability at a tie.
#' Without a tail fit, all remaining upper probability is assigned to the
#' largest bulk knot, so the CDF reaches one rather than losing mass.
#'
#' Because the interpolation starts at \eqn{(0, 0)}, the predicted
#' \eqn{P(Y = 0)} is zero whenever the lowest fitted bulk quantile is positive
#' (\eqn{P(Y = 0)} is then only as large as the bulk probabilities tied at
#' zero). Observed zeros in such rows receive a log score of \code{-Inf};
#' include a low bulk probability (or use a zero-inflated count model) when
#' zeros are common.
#'
#' Following Eq. (3) of Maia, Castro-Camilo, and Browell (2026), with threshold
#' quantile \eqn{q_i = q_{\alpha^{(T)}, i}} the counts
#' \eqn{y \le \lfloor q_i \rfloor} are described by the bulk, up to
#' \eqn{\tilde{\alpha}^{(T)}_i = \widehat{F}_{QR,i}(\lfloor q_i \rfloor)}, and the
#' counts \eqn{y \ge s_i = \max(\lfloor q_i \rfloor + 1, 0)} by the tail:
#' \eqn{\widehat{F}_i(y) = \tilde{\alpha}^{(T)}_i + (1 - \tilde{\alpha}^{(T)}_i) F_{\mathrm{DEGPD}}(y - s_i; \sigma_i, \xi_i, \kappa_i)},
#' avoiding any artificial probability jump at \eqn{y = s_i}. The tail is fitted
#' to the exceedances \eqn{y - s_i} of the observations with \eqn{y \ge s_i}, so
#' the fitted tail law and the predictive tail mass
#' \eqn{1 - \tilde{\alpha}^{(T)}_i = P(Y \ge s_i)} refer to the same event.
#' Here \eqn{q_i} is the fitted threshold after sorting the row's bulk
#' quantiles to remove crossings, not a clipped interpolation
#' knot. A negative threshold gives \eqn{s_i = 0} and
#' \eqn{\tilde{\alpha}^{(T)}_i = 0}: the tail describes the entire row,
#' including zero. For a non-negative, non-integer \eqn{q_i},
#' \eqn{s_i = \lceil q_i \rceil} as in the paper;
#' for an integer \eqn{q_i} (where the paper's two support ranges overlap at
#' \eqn{q_i}) the count \eqn{q_i} belongs to the bulk. The tail shape
#' \eqn{\xi} is parameterised on the unconstrained identity scale
#' (\eqn{\xi \in \mathbb{R}}), allowing light and bounded tails as well as heavy
#' tails. By default the tail optimizer fixes \code{kappa = 1}, which recovers
#' the paper's discrete GPD tail. Supplying \code{egpd_control$fix.arg}
#' replaces that default: include \code{kappa = 1} explicitly to retain it
#' when fixing \code{xi}. The tail scale supports parametric terms only.
#' Offsets enter the log scale with coefficient one. Incomplete rows of the
#' tail model frame are omitted jointly with a warning. A term that replaces
#' a missing input with a finite value does not by itself exclude the row.
#' \code{n_exceed} counts the complete exceedances; at least ten are required.
#' The scale design must have full column rank on these rows, not just
#' on the bulk sample. In particular a training factor level without a
#' complete exceedance cannot identify a separate scale.
#'
#' The \code{qgam} package is an optional dependency and must be installed to
#' use \code{fit_spliced_degpd()}.
#'
#' @param formula a formula for the bulk additive quantile regression. Its
#'   response expression supplies the counts for both the bulk and tail.
#'   Multiple formula offsets are summed in both fitting and prediction.
#' @param data a data frame.
#' @param bulk_probs probabilities used for the bulk quantile regression fit.
#'   With a tail, only knots below \code{threshold_prob} are retained and
#'   the threshold itself is included once.
#' @param threshold_prob probability threshold where the tail model begins.
#'   Ignored when \code{fit_tail = FALSE}: only \code{bulk_probs} supplies knots.
#' @param fit_tail logical; should a DEGPD tail be fitted above the threshold?
#' @param tail_formula NULL, a formula, or a named list of formulae for the
#'   tail model. If a single formula is supplied, it is used for \code{lsigma}
#'   and the shape and kappa submodels default to \code{~1}. If a list is
#'   supplied, valid names are \code{lsigma}, \code{lxi}, and \code{lkappa};
#'   the latter two support only \code{~1}.
#' @param degpd.args list of DEGPD family arguments for the tail fit. Currently
#'   only \code{m = 1} is supported.
#' @param qgam_control optional control list passed to \code{qgam::mqgam()}.
#' @param egpd_control optional list of tail optimizer options. Supported
#'   entries are \code{start}, \code{fix.arg}, \code{optim.method}, and
#'   \code{control}. The \code{start} list accepts \code{beta} (the vector of
#'   coefficients of log scale in model-matrix order), \code{xi} (shape)
#'   and \code{kappa} (positive power). The \code{fix.arg} list accepts
#'   \code{xi} and \code{kappa}; fixed values take precedence over starts.
#'   Unknown option or parameter names are errors. If the tail optimizer finds no parameter value with a
#'   positive likelihood (for example from an infeasible \code{start}), the fit
#'   stops with an error; if it does not converge, a warning of class
#'   \code{"egpd_nonconvergence"} is raised and the tail log-likelihood
#'   (\code{tail_fit$loglik}) is \code{NA}, as for \code{\link{egpd}}.
#' @param keep_data logical; should the training data be stored in the fitted
#'   object?
#'
#' @return \code{fit_spliced_degpd()} returns an object of class
#'   \code{"spliced_degpd"} containing the fitted bulk quantile model, the
#'   fitted DEGPD tail model, the quantile threshold, and the modeling
#'   metadata.
#'
#'   \code{predict.spliced_degpd()} returns a data frame with one row per row
#'   of \code{newdata} (carrying its row names). Missing bulk predictors
#'   give missing predictions only on those rows. During fitting, incomplete
#'   bulk rows (including missing responses) are omitted with a warning
#'   before selecting exceedances; \code{n} counts
#'   the bulk rows used. An unseen fixed-effect
#'   factor level in the bulk predictor gives a warning and \code{NA} for
#'   that row, without changing estimable rows:
#'   \describe{
#'     \item{\code{type = "quantile"}}{discrete integer quantiles
#'       \eqn{\inf\{k \in \mathbb{N}_0 : \widehat{F}_i(k) \ge p\}}, one
#'       column per \code{prob}. Extreme upper percentiles use the survivor,
#'       even where the reported CDF rounds to one. At probability one the
#'       quantile is the last count on bounded support, otherwise \code{Inf}.
#'       Above \eqn{2^{53}} adjacent counts need not be distinct in double
#'       precision;}
#'     \item{\code{type = "cdf"}, \code{"pmf"}}{predicted CDF values or
#'       probability masses with one column per element of \code{at}, in the
#'       order supplied (duplicates are kept);}
#'     \item{\code{type = "parameter"}}{the tail shift \eqn{s_i}
#'       (\code{threshold}), the fitted threshold quantile, the splice
#'       probability \code{alpha_tilde}, and the fitted tail parameters.
#'       Without a tail fit, the two threshold fields are \code{NA};}
#'     \item{\code{type = "bulk"}}{the sorted fitted bulk quantiles.}
#'   }
#' @examples
#' \dontrun{
#' if (requireNamespace("qgam", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 400
#'   x <- runif(n, -1, 1)
#'   y <- rpois(n, lambda = exp(1 + x)) +
#'     rbinom(n, 1, 0.15) * rdiscegpd(n, sigma = 3, xi = 0.2, kappa = 1.4)
#'   df <- data.frame(y = y, x = x)
#'
#'   fit <- fit_spliced_degpd(
#'     y ~ s(x, k = 6),
#'     data = df,
#'     bulk_probs = c(0.1, 0.5, 0.8, 0.9),
#'     threshold_prob = 0.8,
#'     tail_formula = list(lsigma = ~x, lxi = ~1, lkappa = ~1)
#'   )
#'
#'   predict(fit, newdata = df[1:5, ], type = "quantile", prob = c(0.5, 0.9, 0.99))
#'   predict(fit, newdata = df[1:2, ], type = "pmf", at = 0:20)
#' }
#' }
#' @export
fit_spliced_degpd <- function(formula,
                              data,
                              bulk_probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                              threshold_prob = 0.9,
                              fit_tail = TRUE,
                              tail_formula = NULL,
                              degpd.args = list(m = 1),
                              qgam_control = list(),
                              egpd_control = list(),
                              keep_data = FALSE) {
  .spliced_require_qgam()

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("'formula' must be a formula with a response", call. = FALSE)
  lhs <- formula[[2L]]
  response_name <- if (is.symbol(lhs)) as.character(lhs) else deparse1(lhs)
  if (is.symbol(lhs) && !response_name %in% names(data))
    stop("The response variable '", response_name, "' is not a column in 'data'", call. = FALSE)
  y <- eval(lhs, data, environment(formula))
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) != nrow(data) ||
      any(is.infinite(y) | y < 0 | y != floor(y), na.rm = TRUE)) {
    stop("The response must be a finite non-negative integer vector", call. = FALSE)
  }

  if (is.null(degpd.args$m)) {
    degpd.args$m <- 1
  }
  if (!identical(as.integer(degpd.args$m), 1L)) {
    stop("fit_spliced_degpd() currently supports only degpd.args = list(m = 1)", call. = FALSE)
  }

  if (!is.logical(fit_tail) || length(fit_tail) != 1L || is.na(fit_tail))
    stop("'fit_tail' must be TRUE or FALSE.", call. = FALSE)
  bulk_probs <- .normalize_spliced_bulk_probs(
    bulk_probs = bulk_probs,
    threshold_prob = threshold_prob,
    fit_tail = isTRUE(fit_tail)
  )
  if (!fit_tail) threshold_prob <- NA_real_
  threshold_index <- match(threshold_prob, bulk_probs)

  original <- if (isTRUE(keep_data)) data
  bulk_data <- data
  bulk_formula <- formula
  bulk_lhs <- response_name
  if (!is.symbol(lhs) || make.names(response_name) != response_name) {
    ## qgam needs a syntactic response name. Bind the evaluated counts once,
    ## without overwriting a predictor or exposing the binding to the tail.
    bulk_lhs <- .aux_response_name(".egpd_response", c(names(data), all.vars(formula)))
    bulk_formula[[2L]] <- as.name(bulk_lhs)
    bulk_data[[bulk_lhs]] <- y
  }
  cc <- .gam_complete_rows(list(bulk_formula), bulk_data, required = bulk_lhs)
  if (!all(cc)) {
    if (!any(cc)) stop("No complete bulk observations remain.", call. = FALSE)
    warning(sprintf("Omitted %d incomplete bulk observation(s).", sum(!cc)), call. = FALSE)
    data <- data[cc, , drop = FALSE]
    bulk_data <- bulk_data[cc, , drop = FALSE]
    y <- y[cc]
  }
  bulk_formula <- .gam_offset_formula(bulk_formula, bulk_data)
  backend <- .spliced_bulk_terms(bulk_formula, .gam_offset_data(bulk_formula, bulk_data))
  qgam_args <- modifyList(
    list(
      form = backend$formula,
      data = backend$data,
      qu = bulk_probs,
      control = modifyList(list(verbose = FALSE, progress = FALSE), qgam_control)
    ),
    list()
  )
  bulk_fit <- do.call(qgam::mqgam, qgam_args)
  bulk_fit$egpd.bulk.terms <- backend$pred
  if (length(attr(bulk_formula, "egpd.offset", exact = TRUE)))
    bulk_fit$egpd.offset.formula <- bulk_formula

  tail_fit <- NULL
  n_exceed <- 0L
  tail_formula_norm <- NULL

  if (isTRUE(fit_tail)) {
    bulk_train <- .spliced_degpd_sorted_matrix(
      .spliced_degpd_qdo_matrix(bulk_fit, bulk_probs, backend$data)
    )
    threshold_shift <- .spliced_degpd_shift(bulk_train[, threshold_index])
    exceedance <- y - threshold_shift
    tail_data <- data[exceedance >= 0, , drop = FALSE]
    rhs_vars <- if (is.list(tail_formula)) unlist(lapply(tail_formula, all.vars)) else
      all.vars(tail_formula)
    lhs <- .aux_response_name("exceedance", c(names(tail_data), rhs_vars))
    tail_data[[lhs]] <- exceedance[exceedance >= 0]
    tail_formula_norm <- .normalize_spliced_tail_formula(tail_formula, response = lhs)
    tail_fit <- .fit_spliced_degpd_tail(
      tail_formula = tail_formula_norm,
      tail_data = tail_data,
      fit_control = egpd_control,
      training_data = data
    )
    n_exceed <- tail_fit$nobs
  }

  structure(
    list(
      bulk_fit = bulk_fit,
      tail_fit = tail_fit,
      formula = formula,
      tail_formula = tail_formula_norm,
      bulk_probs = bulk_probs,
      threshold_prob = threshold_prob,
      threshold_index = threshold_index,
      degpd.args = degpd.args,
      response_name = response_name,
      n = nrow(data),
      n_exceed = n_exceed,
      call = match.call(),
      data = original
    ),
    class = "spliced_degpd"
  )
}

#' @rdname fit_spliced_degpd
#' @param x an object of class \code{spliced_degpd} or
#'   \code{summary.spliced_degpd}.
#' @param ... unused additional arguments.
#' @export
print.spliced_degpd <- function(x, ...) {
  cat("Spliced discrete EGPD fit\n")
  cat("Bulk probabilities:", paste(signif(x$bulk_probs, 3), collapse = ", "), "\n")
  if (!is.null(x$tail_fit)) cat("Threshold probability:", signif(x$threshold_prob, 3), "\n")
  cat("Tail fitted:", !is.null(x$tail_fit), "\n")
  if (!is.null(x$tail_fit)) {
    cat("Tail exceedances:", x$n_exceed, "\n")
  }
  invisible(x)
}

#' @rdname fit_spliced_degpd
#' @param object an object of class \code{spliced_degpd}.
#' @export
summary.spliced_degpd <- function(object, ...) {
  bulk_edf <- tryCatch(
    sum(vapply(object$bulk_fit$fit, function(x) sum(x$edf), numeric(1))),
    error = function(e) NA_real_
  )
  structure(
    list(
      bulk_probs = object$bulk_probs,
      threshold_prob = object$threshold_prob,
      n = object$n,
      n_exceed = object$n_exceed,
      bulk_edf = bulk_edf,
      has_tail = !is.null(object$tail_fit),
      tail_summary = if (!is.null(object$tail_fit)) {
        data.frame(parameter = c(names(object$tail_fit$coefficients), "xi", "kappa"),
                   estimate = c(unname(object$tail_fit$coefficients),
                                object$tail_fit$xi,
                                object$tail_fit$kappa))
      } else {
        NULL
      },
      call = object$call
    ),
    class = "summary.spliced_degpd"
  )
}

#' @rdname fit_spliced_degpd
#' @export
print.summary.spliced_degpd <- function(x, ...) {
  cat("Spliced discrete EGPD fit\n")
  cat("Bulk probabilities:", paste(signif(x$bulk_probs, 3), collapse = ", "), "\n")
  if (isTRUE(x$has_tail)) cat("Threshold probability:", signif(x$threshold_prob, 3), "\n")
  cat("Observations:", x$n, "\n")
  if (isTRUE(x$has_tail)) {
    cat("Tail exceedances:", x$n_exceed, "\n")
  }
  if (is.finite(x$bulk_edf)) {
    cat("Bulk effective degrees of freedom:", round(x$bulk_edf, 2), "\n")
  }
  if (!is.null(x$tail_summary)) {
    cat("\nTail model summary:\n")
    print(x$tail_summary, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname fit_spliced_degpd
#' @param newdata a data frame of covariates. If omitted, the original training
#'   data are used only when \code{keep_data = TRUE}. An unseen level in a
#'   tail-only factor is an error when predicting the tail.
#' @param type prediction type: \code{"quantile"}, \code{"cdf"},
#'   \code{"pmf"}, \code{"parameter"}, or \code{"bulk"}.
#' @param prob probabilities used for \code{type = "quantile"}.
#' @param at non-negative integer support values used for
#'   \code{type = "cdf"} or \code{"pmf"}.
#' @export
predict.spliced_degpd <- function(object,
                                  newdata,
                                  type = c("quantile", "cdf", "pmf", "parameter", "bulk"),
                                  prob = NULL,
                                  at = NULL,
                                  ...) {
  type <- match.arg(type)
  newdata <- .spliced_degpd_prepare_newdata(object, newdata)
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame", call. = FALSE)
  }

  ## All bulk quantiles share the same design. Filter once, before qdo()
  ## can recycle known rows, then restore their positions for every type.
  bulk_data <- .spliced_bulk_data(object$bulk_fit, newdata)
  ok <- if (nrow(bulk_data)) .gam_estimable_rows(object$bulk_fit$fit[[1L]], bulk_data) else logical()
  rows <- newdata[ok, , drop = FALSE]
  bulk <- .spliced_degpd_sorted_matrix(.spliced_degpd_qdo_matrix(object$bulk_fit, object$bulk_probs,
                                                                 bulk_data[ok, , drop = FALSE]))
  have <- stats::complete.cases(bulk)
  ok[ok] <- have
  pred <- .spliced_predict_rows(object, rows[have, , drop = FALSE], type, prob, at,
                               bulk[have, , drop = FALSE])
  if (all(ok)) return(pred)
  out <- pred[match(seq_along(ok), which(ok)), , drop = FALSE]
  rownames(out) <- rownames(newdata)
  out
}

.spliced_predict_rows <- function(object, newdata, type, prob, at, bulk_mat) {
  row_names <- rownames(newdata)
  threshold_values <- unname(bulk_mat[, object$threshold_index])
  threshold_shift <- .spliced_degpd_shift(threshold_values)

  if (type == "bulk") {
    out <- as.data.frame(bulk_mat, row.names = row_names)
    names(out) <- .egpd_qnames(object$bulk_probs)
    return(out)
  }

  has_tail <- !is.null(object$tail_fit)
  upper_index <- if (!has_tail) length(object$bulk_probs) else object$threshold_index
  upper_prob <- if (!has_tail) 1 else object$threshold_prob
  bulk_probs <- object$bulk_probs[seq_len(upper_index)]
  bulk_sub <- bulk_mat[, seq_len(upper_index), drop = FALSE]
  alpha_tilde <- if (has_tail) {
    .spliced_degpd_alpha_tilde(bulk_sub, bulk_probs, threshold_values, upper_prob)
  } else {
    rep(1, nrow(bulk_sub))
  }

  tail_pars <- if (has_tail) .spliced_degpd_tail_parameters(object, newdata) else NULL

  if (type == "parameter") {
    out <- data.frame(
      threshold = threshold_shift,
      threshold_quantile = threshold_values,
      alpha_tilde = alpha_tilde
    )
    if (!is.null(tail_pars)) {
      out <- cbind(out, tail_pars)
    }
    rownames(out) <- row_names
    return(out)
  }

  if (type == "quantile") {
    if (is.null(prob)) {
      prob <- object$bulk_probs
    }
    prob <- as.numeric(prob)
    if (any(!is.finite(prob)) || any(prob < 0 | prob > 1)) {
      stop("'prob' must contain probabilities in [0, 1]", call. = FALSE)
    }

    out <- matrix(NA_real_, nrow = nrow(bulk_sub), ncol = length(prob))
    for (i in seq_len(nrow(out))) {
      out[i, ] <- .spliced_degpd_row_quantile(
        probs = prob,
        bulk_row = bulk_sub[i, ],
        bulk_probs = bulk_probs,
        upper_prob = upper_prob,
        has_tail = has_tail,
        shift = threshold_shift[i],
        alpha_tilde = alpha_tilde[i],
        scale = if (has_tail) tail_pars$scale[i] else 1,
        shape = if (has_tail) tail_pars$shape[i] else 0,
        kappa = if (has_tail) tail_pars$kappa[i] else 1
      )
    }

    out <- as.data.frame(out, row.names = row_names)
    names(out) <- .egpd_qnames(prob)
    return(out)
  }

  if (is.null(at)) {
    stop("'at' must be supplied for type = 'cdf' or 'pmf'", call. = FALSE)
  }
  at <- as.numeric(at)
  if (!length(at) || any(!is.finite(at)) || any(at != floor(at)) || any(at < 0)) {
    stop("'at' must contain non-negative integers", call. = FALSE)
  }

  ## evaluate on the sorted unique support points, then return one column per
  ## requested value, in the order supplied (duplicates included)
  at_u <- sort(unique(at))
  cols <- match(at, at_u)
  at_lab <- format(at, scientific = FALSE, trim = TRUE)
  if (type == "cdf") {
    cdf <- .spliced_degpd_cdf_matrix(object, newdata, at_u, bulk_mat, tail_pars)
    out <- cdf[, cols, drop = FALSE]
    dimnames(out) <- list(row_names, paste0("F:", at_lab))
    return(as.data.frame(out))
  }

  pmf <- matrix(0, nrow = nrow(bulk_sub), ncol = length(at_u))
  for (i in seq_len(nrow(pmf))) {
    above <- has_tail & at_u >= threshold_shift[i]
    if (any(above)) {
      pmf[i, above] <- (1 - alpha_tilde[i]) * ddiscegpd(
        at_u[above] - threshold_shift[i], sigma = tail_pars$scale[i],
        xi = tail_pars$shape[i], kappa = tail_pars$kappa[i], type = 1
      )
    }
    if (any(!above)) {
      bulk_cdf <- function(k) pmin(alpha_tilde[i], .spliced_degpd_bulk_cdf(
        bulk_sub[i, ], bulk_probs, k, upper_prob
      ))
      k <- at_u[!above]
      pmf[i, !above] <- pmax(bulk_cdf(k) - bulk_cdf(k - 1), 0)
    }
  }
  out <- pmf[, cols, drop = FALSE]
  dimnames(out) <- list(row_names, paste0("P:", at_lab))
  as.data.frame(out)
}
