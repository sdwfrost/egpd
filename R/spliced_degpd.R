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

.spliced_degpd_terms <- function(formula) {
  terms_obj <- terms.formula(formula, specials = c("s", "te", "ti"))
  labels <- attr(terms_obj, "term.labels")
  intercept <- attr(terms_obj, "intercept") == 1
  if (length(labels) == 0L && intercept) {
    labels <- "1"
  }
  labels
}

.spliced_degpd_rhs <- function(formula) {
  reformulate(.spliced_degpd_terms(formula))
}

.spliced_degpd_with_response <- function(formula, response) {
  reformulate(.spliced_degpd_terms(formula), response = response)
}

.spliced_degpd_has_predictors <- function(formula) {
  length(attr(terms.formula(formula), "term.labels")) > 0L
}

.normalize_spliced_bulk_probs <- function(bulk_probs, threshold_prob, fit_tail) {
  if (!is.numeric(bulk_probs) || length(bulk_probs) == 0L || any(!is.finite(bulk_probs))) {
    stop("'bulk_probs' must be a non-empty numeric vector", call. = FALSE)
  }
  if (!is.numeric(threshold_prob) || length(threshold_prob) != 1L || !is.finite(threshold_prob) ||
      threshold_prob <= 0 || threshold_prob >= 1) {
    stop("'threshold_prob' must be a single probability in (0, 1)", call. = FALSE)
  }

  bulk_probs <- sort(unique(as.numeric(bulk_probs)))
  if (any(bulk_probs <= 0 | bulk_probs >= 1)) {
    stop("'bulk_probs' must lie strictly between 0 and 1", call. = FALSE)
  }

  if (fit_tail) {
    bulk_probs <- bulk_probs[bulk_probs <= threshold_prob]
  }
  sort(unique(c(bulk_probs, threshold_prob)))
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
        "'tail_formula' must have names from lsigma, lxi, lkappa, or length 1/3 if unnamed",
        call. = FALSE
      )
    }
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

  tail_formula$lsigma <- .spliced_degpd_with_response(tail_formula$lsigma, response)
  tail_formula$lxi <- .spliced_degpd_rhs(tail_formula$lxi)
  tail_formula$lkappa <- .spliced_degpd_rhs(tail_formula$lkappa)
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

.fit_spliced_degpd_tail <- function(tail_formula, tail_data, fit_control = list()) {
  if (.spliced_degpd_has_specials(tail_formula$lsigma)) {
    stop("tail_formula$lsigma currently supports parametric terms only", call. = FALSE)
  }
  if (.spliced_degpd_has_predictors(tail_formula$lxi)) {
    stop("tail_formula$lxi currently supports only ~1", call. = FALSE)
  }
  if (.spliced_degpd_has_predictors(tail_formula$lkappa)) {
    stop("tail_formula$lkappa currently supports only ~1", call. = FALSE)
  }

  scale_formula <- .spliced_degpd_rhs(tail_formula$lsigma)
  X <- model.matrix(scale_formula, data = tail_data)
  y <- tail_data$exceedance

  fix.arg <- fit_control$fix.arg
  if (is.null(fix.arg)) {
    fix.arg <- list(kappa = 1)
  }
  if (!is.list(fix.arg)) {
    stop("tail optimizer 'fix.arg' must be a list", call. = FALSE)
  }
  bad_fix <- setdiff(names(fix.arg), c("xi", "kappa"))
  if (length(bad_fix) > 0L) {
    stop("Unknown tail fixed parameter(s): ", paste(bad_fix, collapse = ", "), call. = FALSE)
  }

  init_fit <- tryCatch(
    fitegpd(y, type = 1, family = "degpd", fix.arg = fix.arg, hessian = FALSE, optim.method = "BFGS"),
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
  beta0[1] <- log(sigma0)

  start <- fit_control$start
  if (!is.null(start)) {
    if (!is.list(start)) {
      stop("tail optimizer 'start' must be a list", call. = FALSE)
    }
    if (!is.null(start$beta)) {
      if (length(start$beta) != ncol(X)) {
        stop("tail optimizer 'start$beta' must have length equal to the scale design rank", call. = FALSE)
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

  theta0 <- beta0
  if (is.null(fix.arg$xi)) {
    theta0 <- c(theta0, log(pmax(xi0, 1e-8)))
  }
  if (is.null(fix.arg$kappa)) {
    theta0 <- c(theta0, log(pmax(kappa0, 1e-8)))
  }

  unpack <- function(theta) {
    beta <- theta[seq_len(ncol(X))]
    pos <- ncol(X)
    xi <- if (is.null(fix.arg$xi)) {
      pos <- pos + 1L
      exp(theta[pos])
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
    sigma <- drop(exp(X %*% pars$beta))
    dens <- ddiscegpd(y, sigma = sigma, xi = pars$xi, kappa = pars$kappa, type = 1)
    if (any(!is.finite(dens)) || any(dens <= 0)) {
      return(1e20)
    }
    -sum(log(dens))
  }

  optim_method <- fit_control$optim.method
  if (is.null(optim_method)) {
    optim_method <- "BFGS"
  }
  optim_control <- fit_control$control
  if (is.null(optim_control)) {
    optim_control <- list()
  }

  opt <- optim(
    par = theta0,
    fn = nll,
    method = optim_method,
    control = modifyList(list(maxit = 500L), optim_control),
    hessian = FALSE
  )
  pars <- unpack(opt$par)

  structure(
    list(
      formula = scale_formula,
      coefficients = setNames(pars$beta, colnames(X)),
      xi = pars$xi,
      kappa = pars$kappa,
      fix.arg = fix.arg,
      convergence = opt$convergence,
      loglik = -opt$value,
      optim = opt
    ),
    class = "spliced_tail_degpd"
  )
}

.spliced_degpd_response_name <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula", call. = FALSE)
  }
  terms_obj <- terms.formula(formula)
  if (attr(terms_obj, "response") != 1L) {
    stop("'formula' must include a response", call. = FALSE)
  }
  response <- as.character(formula)[2]
  if (!(response %in% names(data))) {
    stop("The response variable '", response, "' is not a column in 'data'", call. = FALSE)
  }
  response
}

.spliced_degpd_qdo_matrix <- function(bulk_fit, probs, newdata) {
  pred <- qgam::qdo(obj = bulk_fit, qu = probs, fun = predict, newdata = newdata)
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
  out <- t(apply(mat, 1, sort))
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = 1L)
  }
  colnames(out) <- colnames(mat)
  rownames(out) <- rownames(mat)
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
  values <- as.numeric(values)
  xout <- as.numeric(xout)

  if (length(probs) == 1L) {
    out <- rep(values[1], length(xout))
    out[xout < probs[1]] <- 0
    return(out)
  }

  stats::approx(
    x = probs,
    y = values,
    xout = xout,
    yleft = 0,
    yright = values[length(values)],
    ties = "ordered"
  )$y
}

.spliced_degpd_bulk_cdf <- function(values, probs, at, upper_prob) {
  values <- as.numeric(values)
  probs <- as.numeric(probs)
  at <- as.numeric(at)

  if (length(values) == 1L) {
    return(ifelse(at < values[1], 0, upper_prob))
  }

  stats::approx(
    x = values,
    y = probs,
    xout = at,
    yleft = 0,
    yright = upper_prob,
    ties = min
  )$y
}

.spliced_degpd_tail_parameters <- function(object, newdata) {
  X <- model.matrix(object$tail_fit$formula, data = newdata)
  data.frame(
    scale = drop(exp(X %*% object$tail_fit$coefficients)),
    shape = rep(object$tail_fit$xi, nrow(newdata)),
    kappa = rep(object$tail_fit$kappa, nrow(newdata))
  )
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

.spliced_degpd_cdf_matrix <- function(object, newdata, at, bulk_mat = NULL, tail_pars = NULL) {
  at <- sort(unique(as.integer(at)))
  if (any(!is.finite(at))) {
    stop("'at' must contain finite integer values", call. = FALSE)
  }

  if (is.null(bulk_mat)) {
    bulk_mat <- .spliced_degpd_sorted_matrix(
      .spliced_degpd_qdo_matrix(object$bulk_fit, object$bulk_probs, newdata)
    )
  }

  upper_index <- if (is.null(object$tail_fit)) length(object$bulk_probs) else object$threshold_index
  upper_prob <- if (is.null(object$tail_fit)) max(object$bulk_probs) else object$threshold_prob
  threshold_idx <- object$threshold_index
  threshold_values <- bulk_mat[, threshold_idx]
  threshold_shift <- pmax(ceiling(threshold_values), 0)
  bulk_probs <- object$bulk_probs[seq_len(upper_index)]
  bulk_mat <- bulk_mat[, seq_len(upper_index), drop = FALSE]

  cdf <- matrix(0, nrow = nrow(bulk_mat), ncol = length(at))
  for (i in seq_len(nrow(bulk_mat))) {
    cdf[i, ] <- .spliced_degpd_bulk_cdf(
      values = bulk_mat[i, ],
      probs = bulk_probs,
      at = at,
      upper_prob = upper_prob
    )
  }

  if (!is.null(object$tail_fit)) {
    if (is.null(tail_pars)) {
      tail_pars <- .spliced_degpd_tail_parameters(object, newdata)
    }
    for (i in seq_len(nrow(cdf))) {
      tail_cdf <- .spliced_degpd_tail_cdf(
        x = at,
        shift = threshold_shift[i],
        scale = tail_pars$scale[i],
        shape = tail_pars$shape[i],
        kappa = tail_pars$kappa[i]
      )
      idx <- at >= threshold_shift[i]
      cdf[i, idx] <- object$threshold_prob +
        (1 - object$threshold_prob) * tail_cdf[idx]
    }
  }

  cdf[cdf < 0] <- 0
  cdf[cdf > 1] <- 1
  cdf <- t(apply(cdf, 1, cummax))
  if (is.null(dim(cdf))) {
    cdf <- matrix(cdf, nrow = 1L)
  }
  colnames(cdf) <- paste0("F:", at)
  rownames(cdf) <- rownames(bulk_mat)
  cdf
}

#' Fit a spliced discrete EGPD forecasting model
#'
#' Fits a hybrid count-distribution model inspired by the X-flexForecast
#' framework of Maia, Castro-Camilo, and Browell (2026): multiple additive
#' quantile regressions describe the bulk of the count distribution, while
#' threshold exceedances are modeled with a type-1 DEGPD tail.
#'
#' @param formula a formula for the bulk additive quantile regression.
#' @param data a data frame.
#' @param bulk_probs probabilities used for the bulk quantile regression fit.
#' @param threshold_prob probability threshold where the tail model begins.
#' @param fit_tail logical; should a DEGPD tail be fitted above the threshold?
#' @param tail_formula NULL, a formula, or a named list of formulae for the
#'   tail model. If a single formula is supplied, it is used for \code{lsigma}
#'   and the shape and kappa submodels default to \code{~1}. If a list is
#'   supplied, valid names are \code{lsigma}, \code{lxi}, and \code{lkappa}.
#' @param degpd.args list of DEGPD family arguments for the tail fit. Currently
#'   only \code{m = 1} is supported.
#' @param qgam_control optional control list passed to \code{qgam::mqgam()}.
#' @param egpd_control optional list of tail optimizer options. Supported
#'   entries are \code{start}, \code{fix.arg}, \code{optim.method}, and
#'   \code{control}.
#' @param keep_data logical; should the training data be stored in the fitted
#'   object?
#'
#' @return An object of class \code{"spliced_degpd"}.
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
  response_name <- .spliced_degpd_response_name(formula, data)
  y <- data[[response_name]]
  if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) {
    stop("The response must be a finite non-negative integer vector", call. = FALSE)
  }

  if (is.null(degpd.args$m)) {
    degpd.args$m <- 1
  }
  if (!identical(as.integer(degpd.args$m), 1L)) {
    stop("fit_spliced_degpd() currently supports only degpd.args = list(m = 1)", call. = FALSE)
  }

  bulk_probs <- .normalize_spliced_bulk_probs(
    bulk_probs = bulk_probs,
    threshold_prob = threshold_prob,
    fit_tail = isTRUE(fit_tail)
  )
  threshold_index <- match(threshold_prob, bulk_probs)

  qgam_args <- modifyList(
    list(
      form = formula,
      data = data,
      qu = bulk_probs,
      control = modifyList(list(verbose = FALSE, progress = FALSE), qgam_control)
    ),
    list()
  )
  bulk_fit <- do.call(qgam::mqgam, qgam_args)

  bulk_train <- .spliced_degpd_sorted_matrix(
    .spliced_degpd_qdo_matrix(bulk_fit, bulk_probs, data)
  )

  tail_fit <- NULL
  n_exceed <- 0L
  tail_formula_norm <- NULL

  if (isTRUE(fit_tail)) {
    threshold_shift <- pmax(ceiling(bulk_train[, threshold_index]), 0)
    exceedance <- y - threshold_shift
    tail_data <- data[exceedance >= 0, , drop = FALSE]
    tail_data$exceedance <- exceedance[exceedance >= 0]
    n_exceed <- nrow(tail_data)
    if (n_exceed < 10L) {
      stop("Not enough exceedances were available to fit the tail model", call. = FALSE)
    }

    tail_formula_norm <- .normalize_spliced_tail_formula(tail_formula, response = "exceedance")
    tail_fit <- .fit_spliced_degpd_tail(
      tail_formula = tail_formula_norm,
      tail_data = tail_data,
      fit_control = egpd_control
    )
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
      data = if (isTRUE(keep_data)) data else NULL
    ),
    class = "spliced_degpd"
  )
}

#' @export
print.spliced_degpd <- function(x, ...) {
  cat("Spliced discrete EGPD fit\n")
  cat("Bulk probabilities:", paste(signif(x$bulk_probs, 3), collapse = ", "), "\n")
  cat("Threshold probability:", signif(x$threshold_prob, 3), "\n")
  cat("Tail fitted:", !is.null(x$tail_fit), "\n")
  if (!is.null(x$tail_fit)) {
    cat("Tail exceedances:", x$n_exceed, "\n")
  }
  invisible(x)
}

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

#' @export
print.summary.spliced_degpd <- function(x, ...) {
  cat("Spliced discrete EGPD fit\n")
  cat("Bulk probabilities:", paste(signif(x$bulk_probs, 3), collapse = ", "), "\n")
  cat("Threshold probability:", signif(x$threshold_prob, 3), "\n")
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

  bulk_mat <- .spliced_degpd_sorted_matrix(
    .spliced_degpd_qdo_matrix(object$bulk_fit, object$bulk_probs, newdata)
  )
  threshold_values <- bulk_mat[, object$threshold_index]
  threshold_shift <- pmax(ceiling(threshold_values), 0)

  if (type == "bulk") {
    out <- as.data.frame(bulk_mat)
    names(out) <- paste0("q:", format(object$bulk_probs, trim = TRUE))
    return(out)
  }

  tail_pars <- if (!is.null(object$tail_fit)) .spliced_degpd_tail_parameters(object, newdata) else NULL

  if (type == "parameter") {
    out <- data.frame(
      threshold = threshold_shift,
      threshold_quantile = as.numeric(threshold_values)
    )
    if (!is.null(tail_pars)) {
      out <- cbind(out, tail_pars)
    }
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

    upper_index <- if (is.null(object$tail_fit)) length(object$bulk_probs) else object$threshold_index
    bulk_probs <- object$bulk_probs[seq_len(upper_index)]
    out <- matrix(NA_real_, nrow = nrow(bulk_mat), ncol = length(prob))

    for (i in seq_len(nrow(out))) {
      bulk_q <- .spliced_degpd_bulk_quantile(
        probs = bulk_probs,
        values = bulk_mat[i, seq_len(upper_index)],
        xout = prob
      )
      if (!is.null(tail_pars)) {
        tail_idx <- prob > object$threshold_prob
        if (any(tail_idx)) {
          p_tail <- (prob[tail_idx] - object$threshold_prob) / (1 - object$threshold_prob)
          bulk_q[tail_idx] <- threshold_shift[i] + qdiscegpd(
            p = p_tail,
            sigma = tail_pars$scale[i],
            xi = tail_pars$shape[i],
            kappa = tail_pars$kappa[i],
            type = 1
          )
        }
      }
      out[i, ] <- bulk_q
    }

    out <- as.data.frame(out)
    names(out) <- paste0("q:", format(prob, trim = TRUE))
    return(out)
  }

  if (is.null(at)) {
    stop("'at' must be supplied for type = 'cdf' or 'pmf'", call. = FALSE)
  }
  if (any(!is.finite(at)) || any(at != floor(at)) || any(at < 0)) {
    stop("'at' must contain non-negative integers", call. = FALSE)
  }

  cdf <- .spliced_degpd_cdf_matrix(
    object = object,
    newdata = newdata,
    at = at,
    bulk_mat = bulk_mat,
    tail_pars = tail_pars
  )

  if (type == "cdf") {
    return(as.data.frame(cdf))
  }

  lower_at <- pmax(sort(unique(as.integer(at))) - 1L, -1L)
  cdf_lower <- .spliced_degpd_cdf_matrix(
    object = object,
    newdata = newdata,
    at = lower_at,
    bulk_mat = bulk_mat,
    tail_pars = tail_pars
  )
  pmf <- cdf - cdf_lower
  pmf[pmf < 0] <- 0
  colnames(pmf) <- paste0("P:", sort(unique(as.integer(at))))
  as.data.frame(pmf)
}
