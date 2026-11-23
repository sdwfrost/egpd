## Setup and fitting helper functions
## Adapted from evgam by Ben Youngman

############ .setup.formulae ##########################

## Parameter-name table for one family. Names given in a named formula list or in
## `fixed` are matched exactly against each parameter's natural name (e.g. "sigma",
## "xi", "pi", "lambda", "sdlog"), its internal and long link-prefixed names
## ("lsigma", "logscale", "logitpi") and the link-prefixed forms of the natural name
## for the parameter's actual link ("l"/"log" for log links, "logit" for logit
## links). Prefixes are never stripped from user input: stripping a link prefix and
## then a leading "l" turned gpois's "logitlambda" into "ambda", so `lambda = ~ 1` was
## rejected. EGPD-type families (those with a shape xi) also accept the scale
## aliases "sigma", "psi" and "scale" and the shape alias "shape". The historical
## spellings "lpi" and "lxi" for parameters that are not log-linked are accepted with
## a deprecation warning.
.param_table <- function(nms, nms2 = NULL, nat = NULL) {
  if (is.null(nms2)) nms2 <- nms
  ## custom families: no natural names are known, so use the link-prefix convention
  if (is.null(nat)) nat <- .strip_link(nms)
  links <- vapply(seq_along(nms), function(i) {
    pre <- if (endsWith(nms[i], nat[i]))
      substr(nms[i], 1L, nchar(nms[i]) - nchar(nat[i])) else ""
    switch(pre, l = , log = "log", logit = "logit", probit = "probit",
           cloglog = "cloglog", "identity")
  }, character(1))
  syn <- as.list(nat)
  if ("xi" %in% nat) {
    i <- which(nat %in% c("sigma", "psi"))
    if (length(i) == 1L) syn[[i]] <- unique(c(nat[i], "sigma", "psi", "scale"))
    j <- which(nat == "xi")
    syn[[j]] <- c("xi", "shape")
  }
  accepted <- lapply(seq_along(nms), function(i) {
    pre <- switch(links[i], log = c("l", "log"), logit = "logit", probit = "probit",
                  cloglog = "cloglog", character(0))
    unique(c(nms[i], nms2[i], syn[[i]], as.vector(outer(pre, syn[[i]], paste0))))
  })
  legacy <- lapply(seq_along(nms), function(i)
    if (links[i] != "log" && nat[i] %in% c("pi", "xi")) paste0("l", nat[i]) else character(0))
  list(nms = nms, nms2 = nms2, nat = nat, links = links,
       accepted = accepted, legacy = legacy)
}

## Index of the parameter that a user-supplied name refers to (see .param_table).
.match_param_name <- function(nm, ptab, what = "parameter") {
  hit <- which(vapply(ptab$accepted, function(a) nm %in% a, logical(1)))
  if (length(hit) > 1L)
    stop(sprintf("'%s' is ambiguous for this family; use one of: %s", nm,
                 paste(ptab$nms[hit], collapse = ", ")), call. = FALSE)
  if (length(hit) == 1L) return(hit)
  hit <- which(vapply(ptab$legacy, function(a) nm %in% a, logical(1)))
  if (length(hit) == 1L) {
    alt <- if (ptab$nms[hit] != ptab$nat[hit]) sprintf(" (or '%s')", ptab$nms[hit]) else ""
    warning(structure(class = c("deprecatedWarning", "warning", "condition"), list(
      message = sprintf("'%s' is a deprecated name for the %s '%s'; use '%s'%s.",
                        nm, what, ptab$nat[hit], ptab$nat[hit], alt),
      call = NULL)))
    return(hit)
  }
  valid <- ifelse(ptab$nat == ptab$nms, ptab$nat,
                  paste0(ptab$nat, " (or ", ptab$nms, ")"))
  stop(sprintf("'%s' is not a parameter of this family; valid names: %s", nm,
               paste(valid, collapse = ", ")), call. = FALSE)
}

## Rebuild a formula from term labels, keeping the environment of the formula it
## replaces: reformulate() would otherwise attach the calling frame, and smooth
## arguments such as s(x, k = kk) that refer to the caller's local variables then
## fail with "object not found".
.reformulate_env <- function(termlabels, response = NULL, intercept = TRUE, like) {
  f <- reformulate(termlabels, response = response, intercept = intercept)
  e <- environment(like)
  if (!is.null(e)) environment(f) <- e
  f
}

.offset_labels <- function(tt) {
  vars <- attr(tt, "variables")
  vapply(attr(tt, "offset"), function(k) deparse1(vars[[k + 1L]]), "")
}

## mgcv's formula interpreter only keeps the first offset. Give it one
## column containing model.offset()'s sum, evaluated in the original
## formula environment on both training data and newdata.
.gam_offset_formula <- function(f, data) {
  tt <- terms.formula(f)
  offs <- .offset_labels(tt)
  if (length(offs) < 2L) return(f)
  nm <- .aux_response_name(".egpd_offset", c(names(data), all.vars(f)))
  backend <- .reformulate_env(c(attr(tt, "term.labels"), paste0("offset(", nm, ")")),
                              if (length(f) == 3L) f[[2L]],
                              intercept = attr(tt, "intercept") == 1L, like = f)
  attr(backend, "egpd.offset") <- setNames(list(.reformulate_env(offs, like = f)), nm)
  backend
}

.gam_offset_data <- function(f, data) {
  spec <- attr(f, "egpd.offset", exact = TRUE)
  if (length(spec)) {
    mf <- stats::model.frame(spec[[1L]], data, na.action = stats::na.pass)
    if (nrow(mf) != nrow(data))
      stop("Evaluated offsets must have one value per data row; supply row-specific offsets as columns.",
           call. = FALSE)
    data[[names(spec)]] <- stats::model.offset(mf)
  }
  data
}

## Only an NA-intolerant term (e.g. poly) needs its inputs screened
## before model.frame. Other terms can legitimately resolve a raw NA.
.gam_complete_rows <- function(formulae, data, required = character()) {
  rhs <- lapply(formulae, function(f)
    stats::as.formula(call("~", f[[length(f)]]), env = environment(f)))
  ok <- stats::complete.cases(data[, intersect(required, names(data)), drop = FALSE])
  frames <- lapply(rhs, function(f) {
    backend <- .gam_offset_formula(f, data)
    ff <- mgcv::interpret.gam(backend)$fake.formula
    environment(ff) <- environment(f)
    list(backend = backend, frame = ff)
  })
  for (fr in frames) {
    spec <- attr(fr$backend, "egpd.offset", exact = TRUE)
    for (f in c(list(fr$frame), unname(spec))) {
      for (v in as.list(attr(stats::terms(f), "variables"))[-1L]) {
        cols <- intersect(.call_vars(v), names(data))
        raw <- stats::complete.cases(data[, cols, drop = FALSE])
        if (!any(ok & !raw)) next
        if (is.symbol(v)) {
          ok <- ok & raw
        } else {
          ## poly() rejects missing inputs before model.frame can omit them.
          ## Inspect the failing basis argument, not every variable in v.
          repeat {
            d <- data[ok, , drop = FALSE]
            keep <- tryCatch({ eval(v, d, environment(f)); NULL }, error = function(err) {
              cl <- conditionCall(err)
              if (!is.call(cl) || !(identical(cl[[1L]], as.name("poly")) ||
                                    identical(cl[[1L]], quote(stats::poly)))) stop(err)
              arg <- match.call(stats::poly, cl, expand.dots = FALSE)$x
              cc <- stats::complete.cases(eval(arg, d, environment(f)))
              if (length(cc) != nrow(d) || all(cc)) stop(err)
              cc
            })
            if (is.null(keep)) break
            ok[ok] <- keep
            if (!any(ok)) break
          }
        }
      }
    }
  }
  if (any(ok)) {
    d <- data[ok, , drop = FALSE]
    keep <- rep(TRUE, nrow(d))
    for (fr in frames) {
      mf <- stats::model.frame(fr$frame, data = .gam_offset_data(fr$backend, d),
                               na.action = stats::na.pass)
      if (nrow(mf) != nrow(d))
        stop("Evaluated predictor must have one row per data row; supply row-specific covariates as columns.",
             call. = FALSE)
      keep <- keep & stats::complete.cases(mf)
    }
    ok[ok] <- keep
  }
  ok
}

.reorder_named_formulae <- function(formula, ptab) {
  nms <- ptab$nms
  if (!is.list(formula) || is.null(nms) || length(formula) != length(nms))
    return(formula)
  fnms <- names(formula)
  if (is.null(fnms) || !any(nzchar(fnms)))
    return(formula)
  pos <- integer(length(formula))
  for (j in seq_along(formula)) {
    nm <- fnms[j]
    if (is.na(nm) || !nzchar(nm)) next
    idx <- .match_param_name(nm, ptab)
    if (idx %in% pos[seq_len(j - 1L)])
      stop(sprintf("Duplicate formula specification for parameter '%s'", nms[idx]),
           call. = FALSE)
    pos[j] <- idx
  }
  unnamed <- which(pos == 0L)
  if (length(unnamed))
    pos[unnamed] <- setdiff(seq_len(length(nms)), pos)
  ## The list is used as it is when it is already in canonical order and its first
  ## formula has the response (or none has). It used to be returned as it was whenever
  ## the order was canonical, so list(lsigma = ~ x, lxi = y ~ 1, lkappa = ~ 1) found no
  ## response, although a named list may carry it on any element.
  two <- vapply(formula, function(f) inherits(f, "formula") && length(f) == 3L, logical(1))
  if (identical(pos, seq_len(length(nms))) && (two[1L] || !any(two))) {
    names(formula) <- nms
    return(formula)
  }
  has_resp <- vapply(formula, function(f) attr(terms.formula(f), "response") == 1L, logical(1))
  ## the response as written (a name, or the backquoted column of a response computed
  ## from an expression, see .materialise_response), not its text, which would be parsed
  ## back into the expression
  resp <- if (any(has_resp)) formula[[which(has_resp)[1L]]][[2L]] else NULL
  out <- vector("list", length(nms))
  for (j in seq_along(formula)) out[[pos[j]]] <- formula[[j]]
  if (!is.null(resp) && attr(terms.formula(out[[1L]]), "response") != 1L) {
    tt <- terms.formula(out[[1L]], specials = c("s", "te", "ti"))
    tl <- attr(tt, "term.labels")
    offs <- .offset_labels(tt)
    ic <- attr(tt, "intercept") == 1L
    if (!length(tl) && ic) tl <- "1"
    out[[1L]] <- .reformulate_env(c(tl, offs), response = resp, intercept = ic,
                                  like = out[[1L]])
  }
  names(out) <- nms
  out
}

## A response written as an expression, such as I(y + 1) or log(y), is evaluated once,
## as by model.frame(): in `data` and the environment of the formula that carries it --
## the first, or for a named list (`named`) the first that has a response, which
## .reorder_named_formulae() then moves to the front. A call is evaluated even when
## `data` has a column named like it (only a name refers to a column), and its value is
## stored under its text unless that would shadow another binding. Every formula
## with that response refers to the resulting column by name, so the expression
## is never evaluated again in a different formula environment.
.materialise_response <- function(formula, data, named = TRUE) {
  single <- inherits(formula, "formula")
  fl <- if (single) list(formula) else formula
  if (!is.list(fl) || !length(fl))
    return(list(formula = formula, data = data))
  two <- vapply(fl, function(f) inherits(f, "formula") && length(f) == 3L, logical(1))
  k <- if (named && !single && any(nzchar(names(fl)))) which(two)[1L] else 1L
  if (is.na(k) || !two[k])
    return(list(formula = formula, data = data))
  lhs <- fl[[k]][[2L]]
  ## (interval-censored responses are rejected by .setup.formulae)
  if (!is.call(lhs) || identical(lhs[[1L]], as.name("cens")))
    return(list(formula = formula, data = data))
  nm <- as.character(fl[[k]])[2L]
  used <- c(names(data), unlist(lapply(fl, all.vars)))
  if (nm %in% used) nm <- .aux_response_name(".egpd_response", used)
  data[[nm]] <- eval(lhs, data, environment(fl[[k]]))
  for (j in which(two))
    if (identical(fl[[j]][[2L]], lhs)) fl[[j]][[2L]] <- as.name(nm)
  list(formula = if (single) fl[[1L]] else fl, data = data)
}

## The data columns a model term reads: the column it names, if it is a name (backquoted
## or not), otherwise the symbols of the expression that are columns of `data` or that
## `env` does not provide (and so must be missing columns). As in base R and mgcv, an
## expression such as log(x) is computed from its variables even when `data` also has a
## column named "log(x)": that column took its place, and predict() then required it in
## `newdata`. Function positions, including pkg::fun, are skipped.
.term_vars <- function(v, data, env) {
  e <- tryCatch(str2lang(v), error = function(err) NULL)
  if (is.symbol(e)) return(as.character(e))
  if (!is.call(e)) return(v)
  vars <- .call_vars(e)
  if (!is.environment(env)) env <- emptyenv()
  vars[vars %in% names(data) | !vapply(vars, exists, logical(1), envir = env)]
}

.call_vars <- function(e) {
  if (is.symbol(e)) {
    v <- as.character(e)
    return(if (nzchar(v)) v else character(0))
  }
  if (!is.call(e)) return(character(0))
  ## The name after $ or @ is not looked up in data or in the formula
  ## environment. Likewise the operands of :: and ::: are literal names.
  head <- e[[1L]]
  if (is.symbol(head)) {
    op <- as.character(head)
    if (op %in% c("$", "@")) return(.call_vars(e[[2L]]))
    if (op %in% c("::", ":::", "quote")) return(character(0))
  }
  unique(unlist(lapply(seq_along(e)[-1L], function(i) .call_vars(e[[i]])),
                use.names = FALSE))
}

.setup.formulae <- function(formula, npar, npar2, data, trace, nms) {
if (inherits(formula, "formula"))
  formula <- lapply(seq_len(npar), function(i) formula)
## evgam's interval-censored responses, cens(left, right) ~ ..., were carried over
## but no family implements a censored likelihood (every d0 stops on it), so reject
## them here rather than after the model has been set up.
lhs <- vapply(formula, function(f)
  if (inherits(f, "formula") && length(f) == 3L) paste(deparse(f[[2L]]), collapse = "") else "",
  character(1))
if (any(startsWith(gsub(" ", "", lhs), "cens(")))
  stop("censored responses, cens(left, right), are not supported: no egpd() family ",
       "implements a censored likelihood", call. = FALSE)
if (npar == 1) {
  if (!(length(formula) %in% c(npar, 1)))
    stop("length(formula) for this family should be 1")
} else {
  if (!(length(formula) %in% c(npar, 1)))
    stop(paste("length(formula) for this family should be", npar, "(or 1 if all parameters are to have the same formula)"))
}
## Variables the formulae read from `data`. mgcv's fake.names lists a parametric term
## as written -- log(x), factor(g), I(x^2), poly(x, 2) -- and every such expression used
## to be rejected as a variable "not supplied to `data'". What the data must supply is
## the variables an expression is computed from (as evgam, which checks all.vars());
## any other symbol in it that the formula's environment provides, such as a degree or
## threshold held in a local variable, is a constant. Smooth terms contribute their
## covariates (and `by` variables) only, so s(x, k = kk) can use a local kk.
terms.list <- lapply(formula, terms.formula, specials=c("s", "te", "ti"))
offterms.list <- lapply(terms.list, .offset_labels)
fake <- lapply(formula, function(f) {
  gp <- mgcv::interpret.gam(.gam_offset_formula(f, data))
  ## mgcv clones the first term for each id, dropping sp declarations on
  ## later terms. Don't silently lose an explicit fixed value that way.
  ss <- gp$smooth.spec
  ids <- vapply(ss, function(s) if (is.null(s$id)) NA_character_ else s$id, "")
  for (j in which(!is.na(ids) & duplicated(ids))) {
    sp <- ss[[j]]$sp
    fixed <- which(sp >= 0)
    if (length(fixed)) {
      first <- ss[[match(ids[j], ids)]]$sp
      if (!identical(as.numeric(first[fixed]), as.numeric(sp[fixed])))
        stop("Specify fixed `sp` on the first smooth with id = ", ids[j], ".", call. = FALSE)
    }
  }
  gp$fake.names
})
pred.vars <- unique(unlist(lapply(seq_along(formula), function(i) {
  fk <- c(fake[[i]][!grepl("^offset\\(", fake[[i]])], offterms.list[[i]])
  unlist(lapply(fk, .term_vars, data = data, env = environment(formula[[i]])))
})))
if (!all(pred.vars %in% names(data))) {
  missing.vars <- pred.vars[!(pred.vars %in% names(data))]
  stop(paste("Variable(s) '", paste(missing.vars, collapse=", "), "' not supplied to `data'.", sep=""))
}
got.specials <- sapply(lapply(terms.list, function(x) unlist(attr(x, "specials"))), any)
termlabels.list <- lapply(terms.list, attr, "term.labels")
got.intercept <- sapply(terms.list, attr, "intercept") == 1
for (i in seq_along(termlabels.list)) {
  if (length(termlabels.list[[i]]) == 0) {
    if (got.intercept[i]) {
      termlabels.list[[i]] <- "1"
    } else if (length(offterms.list[[i]]) > 0) {
      ## offset-only predictor, e.g. ~ offset(z) - 1: a fixed/known linear
      ## predictor with no free coefficients (0-column design). The parameter
      ## is pinned to its offset. Leave term labels empty; the offset is
      ## reattached and the dropped intercept is respected below.
    } else {
      stop(paste("formula element", i, "incorrectly specified"))
    }
  }
}
got.response <- sapply(terms.list, attr, "response") == 1
if (!got.response[1]) {
  stop("formula has no response")
} else {
  ## a name is the column it names (as.character() of the formula would backquote a
  ## non-syntactic one, such as the column of a computed response)
  lhs <- formula[[1]][[2L]]
  response.name <- if (is.symbol(lhs)) as.character(lhs) else as.character(formula[[1]])[2]
}
## Rebuilt formulae keep their original environments (see .reformulate_env), and get
## the response as written in the first formula: a name, or the backquoted column of a
## response computed from an expression (.materialise_response). Its text would be
## parsed back into the expression, to be evaluated again in their environments.
if (any(!got.response)) {
  for (i in which(!got.response)) {
    formula[[i]] <- .reformulate_env(c(termlabels.list[[i]], offterms.list[[i]]),
                                     response=formula[[1]][[2L]], intercept=got.intercept[i],
                                     like=formula[[i]])
  }
}
stripped.formula <- lapply(seq_along(termlabels.list), function(i)
  .reformulate_env(c(termlabels.list[[i]], offterms.list[[i]]),
                   intercept=got.intercept[i], like=formula[[i]]))
attr(formula, "response.name") <- response.name
pred.vars <- pred.vars[!(pred.vars %in% response.name)]
attr(formula, "predictor.names") <- pred.vars
attr(formula, "stripped") <- stripped.formula
attr(formula, "smooths") <- got.specials
attr(formula, "npar") <- npar
for (i in seq_along(formula)) {
  attr(formula[[i]], "intercept") <- got.intercept[i]
  attr(formula[[i]], "smooth") <- got.specials[i]
}
names(formula) <- nms
formula
}

############ .predictable.gam ##########################

## mgcv: log(full.sp) = L %*% log(sp) + lsp0. For fixed and linked
## smooths each row selects one free parameter, or none (a fixed multiplier).
## Keep that index local to this predictor; the id strings are not global.
.sp_index <- function(G) {
  nr <- length(G$lsp0)
  nc <- length(G$sp)
  if (nr != sum(vapply(G$smooth, function(s) length(s$S), 0L)) ||
      any(vapply(G$smooth, function(s) !is.null(s$updateS), logical(1))) ||
      !is.null(G[["H"]]) || any(!is.finite(G$lsp0)))
    stop("Unsupported mgcv penalty setup: expected finite multipliers of smooth penalty matrices.",
         call. = FALSE)
  L <- G[["L"]]
  if (is.null(L) && nr == nc) return(seq_len(nr))
  if (!is.matrix(L) || !all(dim(L) == c(nr, nc)) || !all(L %in% 0:1) ||
      any(rowSums(L) > 1) || any(colSums(L) == 0))
    stop("Unsupported mgcv log-sp map: each penalty must be fixed or select one free parameter.",
         call. = FALSE)
  as.integer(L %*% seq_len(nc))
}

.predictable.gam <- function(G, formula) {
## xlevels and contrasts let predict.gam() rebuild factor-valued parametric terms, such
## as factor(g), on newdata holding only some of their levels
keep <- c("dev.extra", "pterms", "nsdf", "X", "terms", "mf", "smooth", "sp", "term.names", "offset",
          "xlevels", "contrasts", "lsp0", "sp.index")
G$sp.index <- .sp_index(G)
## mgcv sets up some smooths (t2() terms) in a "fit" basis that differs from the
## prediction basis of PredictMat() and predict.gam(), and returns G$P, which maps the
## coefficients of the first to those of the second: G$X = Xpred %*% P, with Xpred the
## prediction-basis design. egpd fits and keeps the coefficients in the fit basis, so
## P is kept to map every design rebuilt by predict.gam() into it (.X.egpd), and
## coefficients handed to mgcv's own code into the prediction basis (plot.egpd).
P <- G[["P"]]
G <- G[keep]
if (!is.null(P)) G$P <- P
## Offset-only predictor (e.g. ~ offset(z) - 1): mgcv returns a 0-column design,
## which would create zero-width Hessian blocks (Armadillo submat out-of-bounds in
## gradHess). Inject a single all-zero column so the parameter is fixed to its
## offset (X %*% beta == 0) without any C++ change; its coefficient is inert (zero
## gradient/Hessian) and is absorbed by the existing pinv/perturb path.
if (ncol(G$X) == 0) {
  G$X <- matrix(0, nrow = nrow(G$X), ncol = 1)
  colnames(G$X) <- "(fixed)"
  G$term.names <- "(fixed)"
}
G$nb <- ncol(G$X)
G$coefficients <- numeric(G$nb)
old <- c("mf", "pP", "cl")
new <- c("model", "paraPen", "call")
is.in <- !is.na(match(old, names(G)))
if (any(is.in)) names(G)[match(old[is.in], names(G))] <- new[is.in]
G$formula <- formula
class(G) <- "gamlist"
G
}

############ .X.egpd ##########################

## predict.gam can convert an unseen fixed-effect level to NA *after*
## its NA handling, then recycle the remaining rows of the model matrix.
## Screen the actual values (not the factor's declared levels). Use the
## fitted parametric terms; new random-effect levels keep mgcv's semantics.
.gam_estimable_rows <- function(o, newdata) {
  ok <- rep(TRUE, nrow(newdata))
  tt <- o$pterms
  pv <- attr(tt, "predvars")
  for (nm in names(o$xlevels)) {
    j <- match(nm, names(attr(tt, "dataClasses")))
    v <- if (identical(attr(newdata, "terms", exact = TRUE), o$terms) && nm %in% names(newdata)) newdata[[nm]] else
      eval(pv[[j + 1L]], newdata, environment(tt))
    if (NROW(v) != nrow(newdata))
      stop("Evaluated factor must have one value per newdata row.", call. = FALSE)
    unseen <- !is.na(v) & !v %in% o$xlevels[[nm]]
    if (any(unseen))
      warning(sprintf("Factor '%s': level(s) %s not in original fit; returning NA for affected rows.",
                      nm, paste(unique(v[unseen]), collapse = ", ")), call. = FALSE)
    ok <- ok & !unseen & !is.na(v)
  }
  ok
}

.gam_predict_rows <- function(o, newdata, type) {
  ## predict.gam restores omitted rows with napredict(), which can drop
  ## the lpmatrix's model.offset attribute for *every* row. Restore the
  ## positions ourselves, using only this predictor's inputs.
  tt <- stats::delete.response(o$terms)
  spec <- attr(o$formula, "egpd.offset", exact = TRUE)
  is_frame <- identical(attr(newdata, "terms", exact = TRUE), o$terms)
  ## A frame made for another predictor is just covariate data. Its terms
  ## must not tell model.matrix/predict.gam the fitted design is ready.
  if (!is_frame) attr(newdata, "terms") <- NULL
  vars <- if (is_frame) names(attr(o$terms, "dataClasses"))[-attr(o$terms, "response")] else
    c(setdiff(.call_vars(tt), names(spec)), if (length(spec)) .call_vars(spec[[1L]]))
  ok <- if (is_frame) stats::complete.cases(newdata[, intersect(vars, names(newdata)), drop = FALSE])
        else rep(TRUE, nrow(newdata))
  if (any(ok)) {
    d <- newdata[ok, , drop = FALSE]
    if (!is_frame) d <- .gam_offset_data(o$formula, d)
    ## Use the fitted terms (including the training polynomial/spline
    ## bases) to catch missingness introduced by a transformation.
    mf <- if (is_frame) d[, intersect(vars, names(d)), drop = FALSE] else
      stats::model.frame(tt, d, na.action = stats::na.pass)
    if (nrow(mf) != nrow(d))
      stop("Evaluated predictor must have one row per newdata row; supply row-specific covariates and offsets as columns.",
           call. = FALSE)
    keep <- .gam_estimable_rows(o, d) & stats::complete.cases(mf)
    ok[ok] <- keep
    if (any(ok)) {
      ## predict.gam's preliminary reformulate(all.vars(...)) reparses
      ## backquoted names as expressions and $ fields as data variables.
      ## The row checks above have already validated the fitted terms.
      nd <- d[keep, , drop = FALSE]
      if (!is_frame) nd <- nd[, intersect(.call_vars(tt), names(nd)), drop = FALSE]
      if (!ncol(nd)) nd[[".row"]] <- seq_len(nrow(nd))
      pred <- mgcv::predict.gam(o, newdata = nd, type = type, newdata.guaranteed = TRUE)
      if (NROW(pred) != nrow(nd) || !length(attr(pred, "model.offset")) %in% c(0L, 1L, nrow(nd)))
        stop("Predicted design and offsets must align with the newdata rows.", call. = FALSE)
      if (all(ok)) return(pred)
    }
  }
  if (type == "lpmatrix") {
    out <- matrix(NA_real_, nrow(newdata), length(o$coefficients),
                  dimnames = list(rownames(newdata), names(o$coefficients)))
    if (any(ok)) {
      out[ok, ] <- if (ncol(pred)) pred else 0  # injected offset-only column
      off <- attr(pred, "model.offset")
      if (length(off)) {
        offsets <- rep(NA_real_, length(ok))
        offsets[ok] <- off
        attr(out, "model.offset") <- offsets
      }
    }
  } else {
    out <- rep(NA_real_, length(ok))
    if (any(ok)) out[ok] <- pred
  }
  out
}

.X.egpd <- function(object, newdata) {
object <- object[sapply(object, inherits, what="gamlist")]
if (missing(newdata)) {
  X <- lapply(object, function(x) x$X)
  offsets <- lapply(object, function(x) as.numeric(x$offset))
} else {
  for (i in seq_along(object))
    class(object[[i]]) <- "gam"
  ## predict.gam() builds the design in the prediction basis; a predictor whose smooths
  ## were set up in another basis (t2(), see .predictable.gam) is mapped into the basis
  ## of its coefficients, X = Xpred %*% P, as on the data it was set up on. Without P
  ## the design of a t2() term on newdata, or on all rows of a compact or maxspline
  ## fit, did not match its coefficients or its penalties.
  X <- lapply(object, function(o) {
    x <- .gam_predict_rows(o, newdata, "lpmatrix")
    P <- o[["P"]]
    if (is.null(P)) return(x)
    off <- attr(x, "model.offset")
    x <- x %*% P
    attr(x, "model.offset") <- off
    x
  })
  ## mirror the offset-only zero-column injection from .predictable.gam so that
  ## newdata prediction matches the fitted design (preserve the model offset).
  X <- lapply(X, function(x) {
    if (ncol(x) == 0) {
      off <- attr(x, "model.offset")
      z <- matrix(0, nrow = nrow(x), ncol = 1)
      if (!is.null(off)) attr(z, "model.offset") <- off
      z
    } else x
  })
  ## mgcv reports a predictor without an offset as the scalar offset 0. Expand it to
  ## one value per row: the likelihood code adds any non-empty offset elementwise,
  ## so a scalar made every maxspline < nrow(data) fit fail with "incompatible matrix
  ## dimensions".
  offsets <- lapply(X, function(x) {
    off <- as.numeric(attr(x, "model.offset"))
    if (length(off) == 1L) rep(off, nrow(x)) else off
  })
}
for (i in seq_along(X))
    colnames(X[[i]]) <- object[[i]]$term.names
names(X) <- names(object)
attr(X, "offsets") <- offsets
X
}

############ .setup.data ##########################

## mgcv::gam(fit = FALSE) for one formula. gam() looks the variables of the model's
## terms that are not in `data` up in its caller's frame, so it is called from the
## formula's environment: a constant in a term, e.g. deg in poly(x, deg), is then found
## where the formula was written, as when calling mgcv::gam() there. (Variables in
## `data` are unaffected.)
.gam_setup <- function(formula, data, knots) {
  env <- environment(formula)
  if (!is.environment(env) || !exists("::", envir = env)) env <- parent.frame()
  ## A response that egpd() computed from an expression is the column named by the
  ## expression's text, which the formula refers to by that backquoted name
  ## (.materialise_response). mgcv rebuilds the formula from text, which turns the name
  ## back into the expression, and evaluates it again, here in this formula's
  ## environment, where a local constant in it may be unknown. The response plays no
  ## part in the set-up, so mgcv reads it from a column with a syntactic name.
  lhs <- if (length(formula) == 3L) formula[[2L]] else NULL
  if (is.symbol(lhs) && !identical(make.names(as.character(lhs)), as.character(lhs)) &&
      as.character(lhs) %in% names(data)) {
    nm <- .aux_response_name(".egpd_response", c(names(data), all.vars(formula)))
    data[[nm]] <- data[[as.character(lhs)]]
    formula[[2L]] <- as.name(nm)
  }
  eval(quote(mgcv::gam(.egpd_formula, data = .egpd_data, fit = FALSE,
                       knots = .egpd_knots, method = "REML")),
       list(.egpd_formula = formula, .egpd_data = data, .egpd_knots = knots), env)
}

## Row indices `id` (a basis subsample of `data`, see `maxspline`) extended by the first
## row of each level, of each factor-like variable of `formula`, that `id` misses: the
## data columns the formula refers to, and the factor-valued expressions among its terms
## and the covariates and `by` variables of its smooths, such as factor(g) in
## s(x, by = factor(g)) (mgcv's fake.formula; its parametric part, `pf`, missed those,
## so a level that the subsample missed lost its smooth without a warning). Unchanged
## when the subsample already contains every level.
.cover_levels <- function(id, formula, data) {
  rhs <- formula[[length(formula)]]
  cols <- lapply(intersect(.call_vars(rhs), names(data)), function(v) data[[v]])
  ff <- tryCatch(mgcv::interpret.gam(formula)$fake.formula, error = function(e) NULL)
  if (!is.null(ff)) {
    ## covariates only: the response is not evaluated again
    if (length(ff) == 3L) ff <- ff[-2L]
    environment(ff) <- environment(formula)
    mf <- tryCatch(stats::model.frame(ff, data = data, na.action = stats::na.pass),
                   error = function(e) NULL)
    if (!is.null(mf)) cols <- c(cols, as.list(mf))
  }
  for (v in cols) {
    if (!(is.factor(v) || is.character(v) || is.logical(v)) || NCOL(v) != 1L)
      next
    miss <- setdiff(unique(v), v[id])
    if (length(miss))
      id <- c(id, match(miss, v))
  }
  id
}

## Families with a probability mass function: their log-likelihood cannot exceed 0.
.egpd_count_families <- c("degpd", "zidegpd", "ztdegpd", "pig", "zipig", "gpig",
                          "gpignat", "zigpig", "zigpignat", "bell", "bellnat",
                          "zibell", "zibellnat", "cmp", "prl", "gpois", "gwaring",
                          "ztgwaring", "plnorm")

.setup.data <- function(data, responsename, formula, family, nms, removeData,
knots, maxdata, maxspline, compact, sargs, outer, trace, gamma) {

## All distributional predictors must use the same observations.
needed <- unique(c(responsename, sargs$id))
complete <- .gam_complete_rows(formula, data, required = needed)
if (any(!complete)) {
  warning(sprintf("Omitting %d rows with missing response or predictor values.",
                  sum(!complete)), call.=FALSE)
  data <- data[complete, , drop=FALSE]
}
if (!nrow(data)) stop("No complete observations remain.", call.=FALSE)
if (family %in% setdiff(.egpd_count_families, c("ztdegpd", "ztgwaring"))) {
  y <- as.matrix(data[, responsename, drop=FALSE])
  if (any(!is.finite(y) | y < 0 | y != floor(y)))
    stop("Discrete families require finite non-negative integer responses.",
         call.=FALSE)
}
if (family %in% c("ztdegpd", "ztgwaring")) {
  y <- as.matrix(data[, responsename, drop=FALSE])
  if (any(!is.finite(y) | y < 1 | y != floor(y)))
    stop("Zero-truncated discrete families require finite positive integer responses (y >= 1).",
         call.=FALSE)
}
## The composite Pareto density is supported on y > 0 only. At y <= 0 its objective
## and analytic gradient disagree, and zeros used to be fitted silently.
if (family == "comppareto") {
  y <- as.matrix(data[, responsename, drop=FALSE])
  if (any(!is.finite(y) | y <= 0))
    stop("family = \"comppareto\" requires finite, strictly positive responses (y > 0).",
         call.=FALSE)
}

if (nrow(data) > maxdata) {
    id <- sort(sample(nrow(data), maxdata))
    data <- data[id, , drop=FALSE]
    if (trace >= 0)
      message("`data' truncated to `maxdata' rows. Re-supply `data' to, e.g., `predict.egpd'")
}

subsampling <- FALSE

## gams
gams <- list()

for (i in seq_along(formula)) {
  f <- .gam_offset_formula(formula[[i]], data)
  d <- .gam_offset_data(f, data)
  if (nrow(data) > maxspline) {
    ## The model is set up on the subsample and its design then evaluated on all rows,
    ## so the subsample must contain every level of every factor in the formula: a
    ## level it missed was dropped from the set-up, and the full-data design no longer
    ## matched it ("number of items to replace is not a multiple of replacement length").
    id <- .cover_levels(sample(nrow(data), maxspline), f, d)
    gams[[i]] <- .gam_setup(f, d[id, , drop=FALSE], knots)
  } else {
    gams[[i]] <- .gam_setup(f, d, knots)
  }
  gams[[i]] <- .predictable.gam(gams[[i]], f)
}

gc()

## likelihood
lik.data <- list()
lik.data$control <- list()
lik.data$outer <- outer
lik.data$control$outer <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-2, stepmax=3)
lik.data$control$inner <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-4, stepmax=1e2)
lik.data$y <- as.matrix(data[,responsename, drop=FALSE])
lik.data$nobs <- nrow(lik.data$y)
## Censored responses are rejected in .setup.formulae; the flag remains because the
## family likelihoods check it.
lik.data$censored <- FALSE
lik.data$sandwich <- !is.null(sargs$id)
if (lik.data$sandwich)
  lik.data$sandwich.split <- data[,sargs$id]
if (lik.data$sandwich && family %in% .egpd_bhhh_families)
  lik.data$observed.step <- .Machine$double.eps^(1/3)

## Evaluate on the fitting rows, never on a deduplicated batch. A term
## or offset such as x - mean(x) depends on the multiplicities.
lik.data$X <- if (nrow(data) > maxspline) .X.egpd(gams, data) else .X.egpd(gams)
lik.data$offsets <- attr(lik.data$X, "offsets")
if (compact) {
  dup.id <- .compact_design_id(lik.data$X, lik.data$offsets, nrow(data))
  unq.id <- match(seq_len(max(dup.id)), dup.id)
  lik.data$X <- lapply(lik.data$X, function(X) X[unq.id, , drop = FALSE])
  attr(lik.data$X, "offsets") <- lapply(lik.data$offsets, function(off)
    if (length(off) > 1L) off[unq.id] else off)
  lik.data$dupid <- dup.id - 1
  lik.data$duplicate <- 1
} else {
  lik.data$dupid <- 0
  lik.data$duplicate <- 0
}
for (i in seq_along(gams)) {
  if (removeData)
    gams[[i]]$y <- NULL
}
if (length(lik.data$X) == 1 & length(nms) > 1) {
  for (i in 2:length(nms)) {
    lik.data$X[[i]] <- lik.data$X[[1]]
    gams[[i]] <- gams[[1]]
    lik.data$offsets[[i]] <- lik.data$offsets[[1]]
  }
}
nbk <- sapply(lik.data$X, ncol)
lik.data$nb <- sum(nbk)
lik.data$idpars <- rep(seq_along(lik.data$X), nbk)
lik.data$LAid <- lik.data$idpars > 0
lik.data$subsampling <- subsampling
gotsmooth <- which(vapply(gams, function(x) length(x$smooth) > 0L, logical(1)))
## gamma tempers the likelihood (k = 1/gamma) for the fit and the smoothing-parameter
## selection; the covariance matrices are then put back on the scale of the unscaled
## likelihood (.VpVc), and the smoothing parameters reported for it (.swap, .finalise)
lik.data$k <- 1 / gamma
lik.data$gamma <- gamma
if (is.null(sargs$id)) {
  lik.data$adjust <- 0
} else {
  if (gamma != 1)
    stop("Can't have gamma != 1 and sandwich adjustment.")
  if (is.null(sargs$method))
    sargs$method <- "magnitude"
  if (sargs$method == "curvature") {
    if (trace > 0)
      message(paste("Sandwich adjustment method: curvature"))
    lik.data$adjust <- 2
  } else {
    if (trace > 0)
      message(paste("Sandwich adjustment method: magnitude"))
    lik.data$adjust <- 1
  }
}
if (is.null(sargs$force))
  sargs$force <- FALSE
lik.data$force <- sargs$force
lik.data$npar <- attr(formula, "npar")
list(lik.data=lik.data, gotsmooth=gotsmooth, data=data, gams=gams, sandwich=lik.data$adjust > 0)
}

## Refine groups by exact equality in the joint evaluated design. Only
## integer labels are serialized, and column/row names are not data.
.compact_design_id <- function(X, offsets, n) {
  id <- rep.int(1L, n)
  first <- rep.int(1L, n)
  for (block in c(X, offsets)) {
    if (!length(block)) next
    block <- as.matrix(block)
    if (nrow(block) == 1L) next
    for (j in seq_len(ncol(block))) {
      v <- block[, j]
      if (all(v == v[first])) next
      key <- paste(id, match(v, unique(v)), sep = ":")
      id <- match(key, unique(key))
      if (max(id) == n) return(id)
      first <- match(id, id)
    }
  }
  id
}

############ .setup.inner.inits ##########################

## A DEGPD shape start on the scale of the shape predictor: log xi under the default
## link, xi itself under the identity link, and under the bounded link
## xi = xi.max (1 - exp(-exp(eta) / xi.max)) its inverse, after moving xi to at most
## `cap` * xi.max. Under the bounded link a start written as log xi lands near the cap
## (the default log xi = 0.05, i.e. xi = 1.05, became xi = 0.199 for xi.max = 0.2),
## on the plateau where the likelihood no longer depends on eta: Newton's steps there
## are erratic and the gradient criterion holds wherever they stop.
.shape_eta <- function(xi, identity = FALSE, xi.max = Inf, cap = 0.5) {
  if (identity)
    return(xi)
  if (is.null(xi.max) || !is.finite(xi.max))
    return(log(xi))
  xi <- pmin(xi, cap * xi.max)
  log(-xi.max * log1p(-xi / xi.max))
}

## Why a fit's bounded DEGPD shape sits on one of its link's flat asymptotes -- xi -> 0
## (eta < -bound) or xi -> xi.max (exp(eta) / xi.max > bound) -- or NULL. There the
## gradient with respect to eta vanishes whether or not the likelihood could still
## increase in xi, so convergence says nothing about optimality.
.shape_plateau <- function(fit, xi.max, bound = 9) {
  if (is.null(xi.max) || !is.finite(xi.max) || isTRUE(fit$likfns$identity.xi))
    return(NULL)
  ## a shape held by `fixed` is recorded under its internal name ("lxi"); the fit's
  ## own element names are the long ones ("logshape"), which never matched
  if (any(c("lxi", "xi") %in% names(fit$fixed)))
    return(NULL)
  eta <- fit[[2L]]$fitted
  if (!length(eta) || anyNA(eta))
    return(NULL)
  if (any(eta < -bound))
    return("has its shape on the xi -> 0 asymptote of the bounded link")
  if (any(exp(eta) / xi.max > bound))
    return(sprintf("has its shape on the xi -> xi.max = %g asymptote of the bounded link", xi.max))
  NULL
}

## Which linear predictors of a fit sit on a plateau of their link, where the likelihood
## no longer depends on them and the parameter is not identified: by the criterion of
## .reml_regular() (R/egpd_fit.R), a predictor other than the scale's (the first) and
## other than one held by `fixed`, with a fitted |eta| > `bound`, or a bounded DEGPD shape
## on either asymptote of its link (.shape_plateau). The Laplace REML criterion is then
## ill-defined: it depends on how far the inner fit drifted along the plateau (continuous
## carrier 5 with kappa -> 0: negREML 447.15 at log kappa -63, 463.22 when refitted at the
## same sp, reaching log kappa -18, with the same log-likelihood). A named logical, one
## per linear predictor, recorded as fit$plateau by .finalise(). The criterion suits the
## EGPD-type families and the generalized Waring, whose parameters other than the scale
## (or mean) are shapes, weights and probabilities on log or logit links; for the other
## families (where, e.g., a CompPareto threshold on the log scale, or a custom family's
## location, can legitimately exceed |eta| = 9) the flags are NA.
.plateau_flags <- function(fit, family, bound = 9) {
  X <- fit$likdata$X
  out <- rep(NA, length(X))
  names(out) <- names(fit)[seq_along(X)]
  if (!(as.vector(family) %in% c("egpd", "degpd", "zidegpd", "ztdegpd", "gwaring", "ztgwaring")))
    return(out)
  free <- vapply(X, function(Xi) ncol(Xi) > 0L && !.is_pinned_block(Xi), logical(1))
  for (i in seq_along(X))
    out[i] <- i > 1L && free[i] && isTRUE(any(abs(fit[[i]]$fitted) > bound))
  if (length(X) > 1L && free[2L] && !is.null(.shape_plateau(fit, fit$likdata$xi.max, bound)))
    out[2L] <- TRUE
  out
}

## The response-scale parameters of the count families whose likelihoods clamp them to
## a box: each linear predictor's name, inverse link and bounds, read from the C++
## (cf_bounds_cpp(), prl_bounds_cpp(), cmp_bounds_cpp()) so that they are the numbers
## the likelihood enforces. Beyond a bound the likelihood is evaluated AT the bound and
## the score in that direction is zero. `wall` marks bounds that are the edge of the
## parameter space itself rather than a numerical limit.
.count_box <- function(family) {
  cf <- function(nm) cf_bounds_cpp()[[nm]]
  box <- function(par, inv, b, wall = NULL)
    list(par = par, inv = inv, lo = unname(b[paste0(par, "_lo")]),
         hi = unname(b[paste0(par, "_hi")]), wall = wall)
  switch(as.vector(family),
    gwaring = , ztgwaring = box(c("mu", "k", "rho"), list(exp, exp, exp), cf("gwaring"),
                                wall = c(rho = "rho > 1")),
    gpois = box(c("mu", "lambda"), list(exp, stats::plogis), cf("gpois")),
    plnorm = box(c("mu", "sigma"), list(exp, exp), cf("plnorm")),
    prl = box("mu", list(exp), prl_bounds_cpp(), wall = c(mu = "mu > 4")),
    cmp = box(c("mu", "nu"), list(exp, exp), cmp_bounds_cpp()),
    NULL)
}

## Each box parameter's value at every observation for coefficients `beta`.
.count_box_values <- function(beta, likdata, bx) {
  b <- split(beta, factor(likdata$idpars, levels = seq_along(likdata$X)))
  lapply(seq_along(bx$par), function(j) {
    X <- likdata$X[[j]]
    eta <- if (ncol(X)) as.vector(X %*% b[[j]]) else numeric(nrow(likdata$y))
    if (ncol(X) && isTRUE(likdata$duplicate == 1L)) eta <- eta[likdata$dupid + 1L]
    if (length(likdata$offsets[[j]])) eta <- eta + likdata$offsets[[j]]
    bx$inv[[j]](eta)
  })
}

## The generalized Waring's rho cannot leave its lower bound once there (the score is
## zero beyond it). A fit that ends there either overshot into the clamp or follows
## data that call for rho <= 1, an infinite mean, which the mean parameterisation
## (a = mu (rho - 1) / k) can only approach as mu grows without bound; in neither case
## did it stop at an optimum. Returns the report, or NULL. (PRL's floor mu > 4 is not
## treated so: a constrained optimum can lie on it.)
.count_wall_report <- function(beta, likdata, family) {
  if (!(as.vector(family) %in% c("gwaring", "ztgwaring")))
    return(NULL)
  bx <- .count_box(family)
  j <- match("rho", bx$par)
  rho <- .count_box_values(beta, likdata, bx)[[j]]
  if (!any(rho <= bx$lo[j]))
    return(NULL)
  sprintf(paste0("rho reached its lower bound %s (rho > 1), where the likelihood is ",
                 "clamped and cannot be improved on in rho: the fit stopped there, not ",
                 "at an optimum"), format(bx$lo[j], digits = 10))
}

## Stop if coefficients `beta` put a count family's parameter outside its box (.count_box)
## for any observation: the optimiser could never move it back, as its score is zero
## there, and it used to stop at once and report convergence (gwaring from
## inits = c(0, 0, 0), i.e. rho = 1: logLik -764.0 against -719.9). A parameter held by
## `fixed` is checked at its fixed value, and the error names `fixed`.
.check_count_box <- function(beta, likdata, family, supplied) {
  bx <- .count_box(family)
  if (is.null(bx))
    return(invisible(NULL))
  vals <- .count_box_values(beta, likdata, bx)
  for (j in seq_along(bx$par)) {
    X <- likdata$X[[j]]
    v <- vals[[j]]
    nm <- bx$par[j]
    bad <- !is.finite(v) | v <= bx$lo[j] | v >= bx$hi[j]
    if (!any(bad))
      next
    pinned <- !ncol(X) || .is_pinned_block(X)
    who <- if (pinned) "`fixed`" else if (supplied) "the supplied `inits`" else
      sprintf("the default starting values for family '%s'", family)
    v1 <- v[bad][1]
    if (nm %in% names(bx$wall) && isTRUE(v1 <= bx$lo[j]))
      stop(sprintf(paste0("%s put %s at %s, at or below its lower bound %s: outside the ",
                          "parameter space (%s)%s"),
                   who, nm, format(v1, digits = 15), format(bx$lo[j], digits = 10),
                   bx$wall[[nm]], if (pinned) "" else "; supply other starting values"),
           call. = FALSE)
    if (pinned)
      stop(sprintf(paste0("`fixed` put %s at %s, outside [%s, %s], the range in which ",
                          "the likelihood is evaluated; fix it inside that range"),
                   nm, format(v1, digits = 10), format(bx$lo[j], digits = 10),
                   format(bx$hi[j], digits = 10)),
           call. = FALSE)
    stop(sprintf(paste0("%s put %s at %s, outside [%s, %s], the range in which the ",
                        "likelihood is evaluated (beyond it the fit cannot move %s); ",
                        "supply other starting values"),
                 who, nm, format(v1, digits = 15), format(bx$lo[j], digits = 10),
                 format(bx$hi[j], digits = 10), nm),
         call. = FALSE)
  }
  invisible(NULL)
}

## Is a design block the single all-zero "(fixed)" column of a parameter pinned by
## `fixed` (an offset-only predictor, see .predictable.gam)? Its coefficient is inert.
.is_pinned_block <- function(X) {
  is.matrix(X) && ncol(X) == 1L && identical(colnames(X), "(fixed)") && all(X == 0)
}

## Every parameter held by `fixed` leaves nothing to estimate; say so instead of
## failing inside the optimiser.
.check_not_all_fixed <- function(fixed, ptab) {
  if (!length(fixed) || is.null(names(fixed)) || any(!nzchar(names(fixed))) ||
      is.null(ptab))
    return(invisible(NULL))
  hit <- suppressWarnings(vapply(names(fixed), .match_param_name, integer(1), ptab = ptab))
  if (all(seq_along(ptab$nms) %in% hit))
    stop("every parameter is held by `fixed`: there is nothing to estimate", call. = FALSE)
  invisible(NULL)
}

## The small preliminary model uses one coefficient per block. Keep offsets,
## fixed parameters and, for a formula without an intercept, its actual first
## column rather than fitting an intercept and transferring it to a slope.
.start_model <- function(likdata) {
  n <- nrow(likdata$y)
  likdata$X <- lapply(likdata$X, function(X) {
    if (!ncol(X)) return(matrix(0, n, 1L, dimnames = list(NULL, "(fixed)")))
    z <- X[, 1L, drop = FALSE]
    if (isTRUE(likdata$duplicate == 1L)) z <- z[likdata$dupid + 1L, , drop = FALSE]
    z
  })
  np <- length(likdata$X)
  likdata$nb <- np
  likdata$dupid <- likdata$duplicate <- 0
  likdata$idpars <- seq_len(np)
  likdata$S <- diag(0, np)
  likdata$CH <- diag(np)
  likdata$compmode <- numeric(np)
  likdata
}

## Project a proposed *total* linear predictor onto the actual start model.
## For an intercept this subtracts the mean offset over the observations.
.start_coefficients <- function(eta, ld, robust = FALSE) {
  vapply(seq_along(ld$X), function(i) {
    z <- ld$X[[i]][, 1L]
    if (all(z == 0)) return(0)
    e <- eta[[i]]
    if (length(ld$offsets[[i]])) e <- e - ld$offsets[[i]]
    if (all(z == 1)) {
      if (robust) stats::median(e) else mean(e)
    } else sum(z * e) / sum(z^2)
  }, 0)
}

## Package-generated restarts are proposals for predictor levels, not user
## coefficients. Translate them only after row selection and offset evaluation.
.eta_start <- function(x, idx = seq_along(x)) {
  if (!is.null(x)) attr(x, "egpd.eta.start") <- idx
  x
}

## `trace` controls the force = TRUE messages of the sandwich adjustment below. It
## used to be read without being an argument, so it resolved to base::trace() and
## every force = TRUE fit that needed the perturbation failed with "comparison (>=)
## is possible only for atomic and list types" (inherited from evgam).
.setup.inner.inits <- function(inits, likdata, likfns, npar, family, trace = 0) {

## Keep the BHHH machinery for obtaining the unadjusted MLE. From that point
## onwards bread, adjusted likelihood curvature, REML and covariance must agree.
obs.step <- likdata$observed.step
likdata$observed.step <- NULL
## Parameters held by `fixed` have one inert, all-zero "(fixed)" design column. The
## intercept-only start fit used to give them a free intercept of ones on top of their
## fixed offset, so it was fitted with those parameters free, and their inert
## coefficient kept that intercept: coef() and summary() then showed an arbitrary
## number (with z = -Inf, "p = 0") for a parameter that was not estimated. Keep them
## pinned in the start fit too, and start their inert coefficient at 0.
pinned <- vapply(likdata$X, function(X) !ncol(X) || .is_pinned_block(X), logical(1))
likdata0 <- .start_model(likdata)

if (is.null(inits)) {
  if (!is.null(likfns$inits)) {
    inits <- likfns$inits(likdata)
  } else {
  if (family %in% c("egpd", "degpd")) {
    inits <- numeric(npar)
    ## The second linear predictor is lxi under the default log link but xi itself
    ## under the identity link (.degpd1idfns, flagged by identity.xi). Seeding both
    ## with 0.05 therefore starts the identity-link fit from xi = 0.05 rather than
    ## xi = exp(0.05) = 1.051 -- a tail 20x lighter. On series whose largest counts
    ## sit far out relative to sigma = mean(y), the DEGPD log-density then underflows
    ## at those counts, log p = -Inf, and the fit aborts with "Some gradient
    ## non-finite" before the optimiser can take a step. Seed the same xi under both,
    ## and under the bounded link (xi.max) map it inside the cap (.shape_eta).
    xm <- if (is.null(likdata$xi.max)) Inf else likdata$xi.max
    xi0 <- if (isTRUE(likfns$identity.xi)) exp(.05) else
      if (is.finite(xm)) .shape_eta(exp(.05), xi.max = xm) else .05
    inits[1:2] <- c(log(mean(likdata$y[,1])), xi0)
    if (attr(family, "type") == 6)
      inits <- c(inits[1:2], 0, 0, 0)
  } else {
    if (family == "zidegpd") {
      inits <- numeric(npar)
      inits[1:2] <- c(log(mean(likdata$y[,1])), .05)
      if (attr(family, "type") == 6)
        inits <- c(inits[1:2], 0, 0, 0, 0)
    } else if (family %in% c("pig", "zipig")) {
      ## PIG/ZIPIG: seed log(mu) at the log mean and log(sigma) at a moment
      ## estimate of the dispersion (var = mu + sigma*mu^2). zipig's logitpi
      ## is seeded near the observed zero proportion.
      y1 <- likdata$y[, 1]
      inits <- numeric(npar)
      inits[1] <- log(max(mean(y1), 1e-3))
      s0 <- (stats::var(y1) - mean(y1)) / max(mean(y1)^2, 1e-8)
      inits[2] <- log(min(max(s0, 0.1), 10))
      if (family == "zipig")
        inits[3] <- stats::qlogis(min(max(mean(y1 == 0), 1e-3), 0.5))
    } else if (family %in% c("gpig", "gpignat", "zigpig", "zigpignat")) {
      ## GPIG/ZIGPIG (Zhu & Joe 2009). Seed a = 1/2 (the PIG special case),
      ## derive c from the empirical dispersion index D via D-1=(1-a)c/(1-c),
      ## and mu (or the implied b) from the sample mean.
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3)
      vv <- stats::var(y1)
      D <- if (is.finite(vv) && ybar > 0) max(vv / ybar, 1.05) else 2
      a0 <- 0.5
      cc <- min(max((D - 1) / ((1 - a0) + (D - 1)), 0.05), 0.95)
      inits <- numeric(npar)
      if (family %in% c("gpig", "zigpig")) {        # mean: logmu, logita, logitc
        inits[1] <- log(ybar)
        inits[2] <- stats::qlogis(a0)
        inits[3] <- stats::qlogis(cc)
      } else {                                      # native: logita, logb, logitc
        b0 <- ybar * (1 - cc)^(1 - a0) / (a0 * cc)
        inits[1] <- stats::qlogis(a0)
        inits[2] <- log(max(b0, 1e-3))
        inits[3] <- stats::qlogis(cc)
      }
      if (family %in% c("zigpig", "zigpignat"))
        inits[4] <- stats::qlogis(min(max(mean(y1 == 0), 1e-3), 0.5))
    } else if (family == "cmp") {
      ## CMP: Var ~ mu/nu asymptotically, so seed nu at 1/DI (clamped into the valid range).
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3); vv <- stats::var(y1)
      DI <- if (is.finite(vv) && vv > 0) max(vv / ybar, 1e-3) else 1
      inits <- c(log(ybar), log(min(max(1 / DI, 1e-5), 50)))
    } else if (family == "prl") {
      ## PRL cannot represent a mean below 4 (mu = tau^2/(tau-1) is minimised at tau = 2).
      y1 <- likdata$y[, 1]
      inits <- log(max(mean(y1), 4 + 1e-3))
    } else if (family == "gpois") {
      ## Generalized Poisson: Var/mean = 1/(1-lambda)^2, so lambda = 1 - sqrt(mean/var).
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3); vv <- stats::var(y1)
      lam0 <- if (is.finite(vv) && vv > ybar) 1 - sqrt(ybar / vv) else 0.1
      lam0 <- min(max(lam0, 0.01), 0.98)
      inits <- c(log(ybar), stats::qlogis(lam0))
    } else if (family %in% c("gwaring", "ztgwaring")) {
      ## Generalized Waring: seed the mean at ybar, k = 1 (the Waring special case), and
      ## rho = 2.5, just inside the finite-variance region rho > 2. The zero-truncated
      ## form has the same parameters but its truncation term rejects rho <= 1, so it
      ## must not fall through to the all-zero start (rho = exp(0) = 1): there the
      ## objective is infinite, its gradient is zeroed, and the optimiser used to stop
      ## at once and report convergence.
      y1 <- likdata$y[, 1]
      inits <- c(log(max(mean(y1), 1e-3)), 0, log(2.5))
    } else if (family == "plnorm") {
      ## Poisson-lognormal: matching the first two moments gives
      ## sigma^2 = log(1 + (var - mean)/mean^2).
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3); vv <- stats::var(y1)
      s2 <- if (is.finite(vv) && vv > ybar) log(1 + (vv - ybar) / ybar^2) else 0.25
      inits <- c(log(ybar), 0.5 * log(min(max(s2, 1e-4), 25)))
    } else if (family %in% c("bell", "bellnat", "zibell", "zibellnat")) {
      ## Bell/ZIBell (Castellares et al. 2018). Native theta is seeded at the
      ## MLE theta.hat = W0(ybar); the mean parameterisation at log(ybar).
      ## ZI variants seed logitpi near the observed zero proportion.
      y1 <- likdata$y[, 1]
      ybar <- max(mean(y1), 1e-3)
      inits <- numeric(npar)
      if (family %in% c("bell", "zibell")) {          # mean: log mu
        inits[1] <- log(ybar)
      } else {                                        # native: log theta
        inits[1] <- log(max(bell_W0_cpp(ybar), 1e-3))
      }
      if (family %in% c("zibell", "zibellnat"))
        inits[2] <- stats::qlogis(min(max(mean(y1 == 0), 1e-3), 0.5))
    } else {
      inits <- numeric(npar)
      ## ZT-DEGPD starts from zeros, i.e. xi = 1 under the log link; under the
      ## bounded link eta = 0 would sit on the plateau near xi.max
      if (family == "ztdegpd" && !isTRUE(likfns$identity.xi) &&
          is.finite(if (is.null(likdata$xi.max)) Inf else likdata$xi.max))
        inits[2] <- .shape_eta(1, xi.max = likdata$xi.max)
    }
  }
  }
  ## A custom initialiser may already return coefficients; likewise the
  ## specialised ZT-Waring start has already fitted this offset model.
  project <- is.null(likfns$inits) || family == "comppareto"
  if (project) {
    eta0 <- inits
    inits <- .start_coefficients(eta0, likdata0)
  }
  inits[pinned] <- 0
  .check_count_box(inits, likdata0, family, supplied = FALSE)
  if (.egpd_invalid_objective(.nllh.nopen(inits, likdata0, likfns))) {
    ## One extreme offset can pull the mean so far that other rows overflow.
    ## A median-centred intercept still translates with the exposure units.
    if (project && any(lengths(likdata0$offsets) > 0L))
      inits <- .start_coefficients(eta0, likdata0, robust = TRUE)
    .check_count_box(inits, likdata0, family, supplied = FALSE)
    if (.egpd_invalid_objective(.nllh.nopen(inits, likdata0, likfns)))
      stop(sprintf(paste0("the default starting values for family '%s' give a non-finite ",
                          "log-likelihood on these data; supply `inits`"), family),
           call. = FALSE)
  }
  beta0 <- .newton_step_inner(inits, .nllh.nopen, .search.nopen, likdata=likdata0, likfns=likfns, control=likdata$control$inner)$par
  inits.supplied <- FALSE
} else {
  inits.supplied <- TRUE
  if (is.list(inits)) {
    betamat <- expand.grid(inits)
    betanllh <- numeric(nrow(betamat))
    for (i in seq_len(nrow(betamat))) {
      beta0 <- unlist(betamat[i,])
      betanllh[i] <- likfns$d0(beta0, likdata0)
    }
    beta0 <- betamat[which.min(betanllh),]
  } else {
    beta0 <- inits
    eta <- attr(inits, "egpd.eta.start")
    if (length(eta)) beta0[eta] <- .start_coefficients(inits, likdata0)[eta]
  }
}
beta0 <- unlist(lapply(seq_len(npar), function(i) {
  nci <- ncol(likdata$X[[i]])
  ## a 0-column (offset-only) block contributes no coefficients; a pinned block's
  ## coefficient is inert and starts (and so stays) at 0
  if (nci == 0) numeric(0) else if (pinned[i]) 0 else c(beta0[i], rep(0, nci - 1))
}))
compmode <- 0 * beta0
CH <- diag(compmode + 1)
k <- likdata$k
likdata[c("k", "CH", "compmode")] <- list(k, CH, compmode)
## A start outside a count family's clamp box cannot be improved on in that parameter
## (its score is zero there), nor can one at which the likelihood cannot be evaluated:
## every objective maps such points to 1e20 with a zero gradient. Either way the
## optimiser would stop at once, so fail as optim() does instead of returning an
## unmoved fit. (After the default start's intercept-only fit a parameter may sit
## beyond a bound legitimately, e.g. gpois lambda at its Poisson limit, so only the
## supplied values are checked against the box here.)
if (inits.supplied)
  .check_count_box(beta0, likdata, family, supplied = TRUE)
if (.egpd_invalid_objective(.nllh.nopen(beta0, likdata, likfns)))
  stop(if (inits.supplied)
         "the supplied `inits` give a non-finite log-likelihood on these data; supply other starting values"
       else
         sprintf(paste0("the default starting values for family '%s' give a non-finite ",
                        "log-likelihood on these data; supply `inits`"), family),
       call. = FALSE)
diagH <- diag(.gH.nopen(beta0, likdata=likdata, likfns=likfns)[[2]])
if (likdata$sandwich) {
  mle <- .fit_inner(beta0, .nllh.nopen, .search.nopen, likdata, likfns)
  if (mle$convergence != 0L)
    stop("Cannot compute sandwich adjustment: unpenalised fit did not converge (",
         mle$report, ").", call. = FALSE)
  beta0 <- mle$par
  likdata$observed.step <- obs.step
  H <- .gH.nopen(beta0, likdata=likdata, likfns=likfns, sandwich=TRUE)
  scores <- rowsum(t(H[[1]]), likdata$sandwich.split, reorder = FALSE)
  if (nrow(scores) < 2L)
    stop("Sandwich adjustment needs at least two independent clusters.", call. = FALSE)
  J <- crossprod(scores)
  ## Compare also to the individual scores: when all cluster sums vanish at
  ## the MLE (cluster fixed effects), optimizer residuals are not information.
  hcJ <- .hessian_check(J, J + tcrossprod(H[[1]]), likdata)
  hcJ$rank <- min(hcJ$rank, nrow(scores) - 1L)  # all cluster scores sum to zero
  H <- H[[2]]
  diagH <- diag(H)
  curvature <- likdata$adjust == 2L
  if (curvature || !is.null(obs.step))
    hcH <- .hessian_check(H, H, likdata)
  if (!is.null(obs.step))
    .check_observed_bread(beta0, H, H, likdata, likfns, hcH)
  force.curv <- curvature && is.null(obs.step) && isTRUE(likdata$force)
  if (is.na(hcJ$rank) || (curvature && !force.curv))
    .require_sandwich_rank(hcJ, "cluster score covariance")
  if (hcJ$rank < hcJ$npar)
    warning(warningCondition(sprintf("The cluster score covariance is rank deficient (rank %d/%d): %s.",
      hcJ$rank, hcJ$npar, if (curvature)
        "forced pseudo-inverse inference is unsupported in the missing directions"
      else "the magnitude adjustment uses its trace, not a full-rank sandwich covariance"),
      class = "egpd_cluster_rank", call = NULL))
  ## The inert coefficients of parameters held by `fixed` have exactly zero rows and
  ## columns in H and J. Left in, they made H singular, so the adjustment failed (and
  ## with force = TRUE the magnitude k counted them, which inflated it by nb/(nb - q)
  ## and shrank the standard errors). Adjust the free coefficients only: k is their
  ## number over tr(H^-1 J), and CH leaves the inert coefficients alone.
  free <- !(likdata$idpars %in% which(pinned))
  H <- H[free, free, drop=FALSE]
  J <- J[free, free, drop=FALSE]
  if (curvature) {
    if (!force.curv || is.na(hcH$rank))
      .require_sandwich_rank(hcH, "observed information")
    if (hcH$rank < hcH$npar || hcJ$rank < hcJ$npar) {
      if (hcH$rank < hcH$npar)
        warning("The unpenalised Hessian does not support full-rank curvature inference.",
                call. = FALSE)
      if (trace >= 0)
        message("Using a pseudo-inverse for the curvature adjustment (force = TRUE).")
      CH[free, free] <- .forced_curvature(H, J, hcH$rank, hcJ$rank)
    } else {
      ## CH' H CH = H J^-1 H without forming/inverting the sandwich.
      A <- backsolve(chol(J), H, transpose = TRUE)
      if (!is.null(obs.step)) {
        CH[free, free] <- backsolve(chol(H), A)
      } else {
        ## Keep the existing orientation of the exact-Hessian adjusted likelihood.
        sH <- svd(H)
        sHA <- svd(crossprod(A))
        CH[free, free] <- solve(sqrt(sH$d) * t(sH$v), sqrt(sHA$d) * t(sHA$v))
      }
    }
    compmode <- beta0
  } else {
    cholH <- try(chol(H), silent=TRUE)
    if (inherits(cholH, "try-error")) {
      if (!likdata$force)
        stop("Hessian of unpenalised MLE not positive definite.\n  Supply `force=TRUE' to `sandwich.args' for a pseudo-inverse.")
      if (trace >= 0)
        message("Hessian of the unpenalised MLE is not positive definite; using its ",
                "pseudo-inverse for the sandwich adjustment (force = TRUE).")
      iH <- pinv(H)
    } else iH <- chol2inv(cholH)
    k <- 1 / mean(diag(crossprod(iH, J)))
    ## with a pseudo-inverse (force = TRUE) of an indefinite Hessian this need not be
    ## a valid scaling, and a non-positive k would turn the fit into a maximisation
    if (!is.finite(k) || k <= 0)
      stop(sprintf(paste0("the sandwich magnitude adjustment is not positive (k = %g); ",
                          "the unpenalised fit is not at a maximum of the likelihood"), k),
           call. = FALSE)
  }
  if (!is.null(obs.step)) {
    C <- CH[free, free, drop = FALSE]
    diagH[free] <- k * diag(crossprod(C, H %*% C))
  }
}
attr(beta0, "k") <- k
attr(beta0, "CH") <- CH
attr(beta0, "compmode") <- compmode
attr(beta0, "diagH") <- diagH
beta0
}

############ .guess ##########################

.guess <- function(x, d, s) {
okay <- s != 0
val <- d / (d + exp(x) * s)
mean(val[okay]) - .4
}

############ .sandwich ##########################

.sandwich <- function(likdata, beta) {
likdata$k <- attr(beta, "k")
likdata$CH <- attr(beta, "CH")
likdata$compmode <- attr(beta, "compmode")
bigX <- do.call(cbind, likdata$X)
CHX <- bigX %*% likdata$CH
CHX <- lapply(unique(likdata$idpars), function(i) CHX[,likdata$idpars == i])
likdata$CHX <- CHX
likdata
}

############ .outer ##########################

.outer <- function(rho0, beta, likfns, likdata, Sdata, control, correctV, outer, trace) {

tryCatch(
  .outer.fit(rho0, beta, likfns, likdata, Sdata, control, correctV, outer, trace),
  egpd_d34_unavailable = function(e) {
    warning(structure(class = c("egpd_d34_fallback", "warning", "condition"),
                      list(message = paste0(conditionMessage(e),
                                            "; using finite-difference outer optimization."),
                           call = NULL)))
    likfns$d340 <- NULL
    .outer.fit(rho0, beta, likfns, likdata, Sdata, control, correctV, "fd", trace)
  }
)

}

.outer.fit <- function(rho0, beta, likfns, likdata, Sdata, control, correctV, outer, trace) {

attr(rho0, "beta") <- beta

if (outer == "fixed") {

  fit.reml <- .reml0_fixed(rho0, likfns=likfns, likdata=likdata, Sdata=Sdata)
  ## .reml0_fixed minimises the penalised negative log-likelihood, which used to be
  ## reported as negREML. Report the REML criterion at the supplied smoothing
  ## parameters instead, as for estimated ones.
  laml <- .reml0(rho0, likfns, likdata, Sdata, beta=fit.reml$beta)
  fit.reml$negREML <- if (.egpd_invalid_objective(laml)) NA_real_ else as.numeric(laml)

} else {

if (is.null(likfns$d340) & outer != "fd")
  outer <- "fd"

if (outer == "newton") {
  fit.reml <- .newton_step_inner(rho0, .reml0, .search.reml, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
} else {
  if (outer == "fd") {
    ## .reml1.fd tightens the inner tolerance for its own differences. The line
    ## search compares objective values at the default tolerance, which is cheaper
    ## (tightening every evaluation doubled the cost of these fits); if it then runs
    ## out of resolution before the gradient criterion is met -- typically on a flat
    ## REML plateau as a smoothing parameter grows -- finish from there with
    ## objective values as accurate as the differences.
    fit.reml <- .BFGS(rho0, .reml0, .reml1.fd, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
    likdata_fd <- .tighten_inner(likdata, 1e-8)
    if (fit.reml$convergence == 0L &&
        !identical(fit.reml$report, "gradient tolerance reached") &&
        !identical(likdata_fd$control$inner$gradtol, likdata$control$inner$gradtol)) {
      rho.mid <- fit.reml$par
      attr(rho.mid, "beta") <- fit.reml$beta
      fit.reml <- .BFGS(rho.mid, .reml0, .reml1.fd, likfns=likfns, likdata=likdata_fd, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
    }
  } else {
    fit.reml <- .BFGS(rho0, .reml0, .reml1, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
  }
  rho1 <- fit.reml$par
  attr(rho1, "beta") <- fit.reml$beta
  if (correctV) {
    ## finite differences of the analytic REML gradient where third derivatives
    ## exist, otherwise Richardson-extrapolated second differences of the criterion
    fit.reml$Hessian <- try(.reml12(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)[[2]], silent=TRUE)
    if (inherits(fit.reml$Hessian, "try-error"))
      fit.reml$Hessian <- try(.reml2.fdfd(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata), silent=TRUE)
  }
}

if (correctV && outer != "fixed") {
  ## A failure here only affects the smoothing-parameter correction to Vc, so it
  ## must not discard an otherwise complete fit.
  fit.reml$invHessian <- tryCatch(
    if (inherits(fit.reml$Hessian, "try-error"))
      stop(attr(fit.reml$Hessian, "condition")) else .reml_vcov(fit.reml$Hessian),
    error = function(e) {
      warning("the REML Hessian could not be computed (", conditionMessage(e),
              "); Vc omits smoothing-parameter uncertainty and equals Vp", call.=FALSE)
      NULL
    })
}

fit.reml$trace <- trace

if (trace == 1) {
  report <- "\n Final max(|grad|))"
  likdata$S <- .makeS(Sdata, exp(fit.reml$par))
  report <- c(report, paste("   Inner:", signif(max(abs(.gH.pen(fit.reml$beta, likdata, likfns)[[1]])), 3)))
  report <- c(report, paste("   Outer:", signif(max(abs(fit.reml$gradient)), 3)))
  report <- c(report, "", "")
  cat(paste(report, collapse="\n"))
}

}

fit.reml

}

############ .outer.nosmooth ##########################

.outer.nosmooth <- function(beta, likfns, likdata, control, trace) {

fit.inner <- .fit_inner(beta, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns)

fit.inner$beta <- fit.inner$par
fit.inner

}

############ .VpVc ##########################

.VpVc <- function(fitreml, likfns, likdata, Sdata, correctV, sandwich, smooths, trace) {
lsp <- fitreml$par
H0 <- .gH.nopen(fitreml$beta, likdata, likfns)[[2]]
if (smooths) {
  sp <- exp(lsp)
  H <- H0 + likdata$S
} else {
  H <- H0
}
## Coefficients pinned by an offset-only predictor (`fixed=`, see .predictable.gam)
## carry an all-zero design column, so their row and column of H are exactly zero and
## H is singular. That direction is a nuisance -- orthogonal to every estimated
## coefficient -- but left alone it corrupts the whole covariance matrix through
## pinv() and the smoothing-parameter correction below; on one test fit it produced a
## negative variance for the shape parameter, so confint() returned NA. Giving those
## directions a unit diagonal makes H block-diagonal and invertible without altering
## the inverse of the free block; their rows and columns are zeroed out afterwards,
## which also leaves them exactly zero effective degrees of freedom in .edf().
inert <- which(rowSums(abs(H)) == 0 & colSums(abs(H)) == 0)
## Whether H is singular, or not positive definite, on the coefficients that are
## estimated: .finalise() records the rank and warns if the standard errors cannot show
## it (.warn_singular_hessian). This replaces a message() that only caught an H failing
## chol(), which an exactly singular H can pass through rounding. A diagnostic that fails
## (e.g. an eigen() error) leaves the rank NA, "not assessed", rather than failing the fit.
hessian <- tryCatch(.hessian_check(H, H0, likdata),
                    error = function(e) list(rank = NA_integer_, npar = NA_integer_))
if (!is.null(likdata$observed.step))
  .check_observed_bread(fitreml$beta, H, H0, likdata, likfns, hessian)
else if (isTRUE(likdata$adjust == 2L) && !likdata$force)
  .require_sandwich_rank(hessian, "adjusted information at the fitted coefficients")
if (length(inert)) diag(H)[inert] <- 1

Vc <- Vp <- pinv(H)
## gamma (as in mgcv) only makes the smoothing-parameter selection smoother. The fit
## minimises the log-likelihood divided by gamma plus the penalty, i.e. the unscaled
## negative log-likelihood plus the penalty with smoothing parameters gamma * sp, and
## the covariance of the coefficients is that of this unscaled penalised likelihood,
## (H_raw + gamma S)^-1 = (H_raw / gamma + S)^-1 / gamma. Leaving it at
## (H_raw / gamma + S)^-1 inflated every standard error and interval by sqrt(gamma).
## The derivatives of the coefficients with respect to the log smoothing parameters
## and the REML covariance of those (Vc below) do not depend on the scaling.
g <- if (is.null(likdata$gamma)) 1 else likdata$gamma
if (g != 1)
  Vc <- Vp <- Vp / g
if (smooths) {
if (correctV) {
## Vc corrects Vp for uncertainty in the smoothing parameters. That correction only
## exists if they were estimated: when all are fixed (in the terms or by `sp`)
## they carry no uncertainty, fitreml then has no invHessian to propagate,
## and the correction term is exactly zero, so Vc = Vp. Previously this branch ran
## regardless and failed on the missing invHessian ("requires numeric/complex matrix
## arguments"), which made any `sp = ` fit unusable at the default correctV = TRUE.
Vrho <- fitreml$invHessian
if (is.null(Vrho) || !is.matrix(Vrho) || nrow(Vrho) != length(lsp)) {
  Vrho <- 0
} else {
  cholVp <- try(chol(Vp), silent=TRUE)
  if (inherits(cholVp, "try-error")) {
      cholVp <- attr(.perturb(Vp), "chol")
  }
  attr(lsp, "beta") <- fitreml$beta
  spSl <- Map("*", attr(Sdata, "Sl"), exp(lsp))
  dbeta <- .d1beta(lsp, fitreml$beta, spSl, .Hdata(H))$d1
  Vbetarho <- tcrossprod(dbeta %*% Vrho, dbeta)
  VR <- matrix(0, nrow=likdata$nb, ncol=likdata$nb)
  Vc <- .perturb(Vp + Vbetarho + VR)
}
} else {
  Vrho <- 0
}
} else {
  Vrho <- 0
}
## a pinned coefficient has no sampling variability and no covariance with anything
if (length(inert)) {
  Vp[inert, ] <- 0; Vp[, inert] <- 0
  Vc[inert, ] <- 0; Vc[, inert] <- 0
}
## Hessians of the unscaled likelihood, to match Vp in .edf()
if (g != 1) {
  H0 <- g * H0
  H <- g * H
}
list(Vp=Vp, Vc=Vc, Vlsp=Vrho, H0=H0, H=H, hessian=hessian)
}

############ .hessian_check ##########################

## Families whose Hessian is the outer product of the per-observation scores (BHHH, a
## documented decision): those of src/countfams.cpp (gpois, gwaring, plnorm), GPIG and
## ZIGPIG (src/gpig.cpp), CMP, PRL, and the zero-truncated generalized Waring
## (R/hurdle.R). The other families use the exact observed information.
.egpd_bhhh_families <- c("gpois", "gwaring", "ztgwaring", "plnorm", "gpig", "gpignat",
                         "zigpig", "zigpignat", "cmp", "prl")

.require_sandwich_rank <- function(hc, what) {
  if (!isTRUE(hc$rank == hc$npar))
    stop(sprintf("Cannot compute sandwich adjustment: %s is non-finite, numerically singular or not positive definite (rank %s/%s).",
                 what, hc$rank, hc$npar), call. = FALSE)
}

## Legacy exact-Hessian force option. Build the map on the resolved positive
## subspaces, so a chol() that happens to pass on roundoff cannot invert an
## unsupported direction. There C'H C = (H^+ J H^+)^+.
.forced_curvature <- function(H, J, rankH, rankJ) {
  if (min(rankH, rankJ) < 1L)
    stop("Cannot force a curvature adjustment of zero rank.", call. = FALSE)
  eh <- eigen(H, symmetric = TRUE)
  ej <- eigen(J, symmetric = TRUE)
  ih <- seq_len(rankH); ij <- seq_len(rankJ)
  U <- sweep(eh$vectors[, ih, drop = FALSE], 2L, sqrt(eh$values[ih]), "/")
  Z <- sweep(ej$vectors[, ij, drop = FALSE], 2L, sqrt(ej$values[ij]), "*")
  v <- svd(tcrossprod(U) %*% Z, nv = 0L)
  r <- seq_len(min(rankH, sum(v$d > max(v$d) * sqrt(.Machine$double.eps))))
  if (!length(r)) stop("Cannot force a curvature adjustment of zero rank.", call. = FALSE)
  U[, r, drop = FALSE] %*% (t(v$u[, r, drop = FALSE]) / v$d[r])
}

## Test resolution in the metric of the requested covariance: compare two
## score differences in every curvature direction, not just the largest one.
.check_observed_bread <- function(beta, H, H0, ld, lf, hc = .hessian_check(H, H0, ld)) {
  .require_sandwich_rank(hc, "observed information")
  R <- chol(H[hc$free, hc$free, drop = FALSE])
  ld$observed.step <- 2 * ld$observed.step
  E <- (.hess.nopen(beta, ld, lf) - H0)[hc$free, hc$free, drop = FALSE]
  E <- backsolve(R, t(backsolve(R, E, transpose = TRUE)), transpose = TRUE)
  if (!all(is.finite(E)) || norm(E, "2") > 0.01)
    stop("Cannot compute sandwich adjustment: observed information is unstable across score-difference steps.",
         call. = FALSE)
}

## Each coefficient's scale: the root mean square of its design column over the
## observations (1 for a column of zeros). Measured in these units every coefficient moves
## its linear predictor by about one unit, so the rank test of .hessian_check() does not
## depend on the units of the covariates (an unscaled test called `y ~ year`, with year
## around 2010, singular), while the intercepts of intercept-only predictors -- the
## coefficients of the generalized Waring a = k case -- are left as they are.
.coef_scale <- function(likdata) {
  unlist(lapply(likdata$X, function(X) {
    if (!ncol(X)) return(numeric(0))
    cnt <- if (isTRUE(likdata$duplicate == 1)) tabulate(likdata$dupid + 1L, nbins = nrow(X))
    ms <- vapply(seq_len(ncol(X)), function(j)
      if (is.null(cnt)) mean(X[, j]^2) else sum(cnt * X[, j]^2) / sum(cnt), 0)
    s <- sqrt(ms)
    s[!is.finite(s) | s <= 0] <- 1
    s
  }))
}

## Is the Hessian H whose pseudo-inverse is Vp (.VpVc: H0, the Hessian of the negative
## log-likelihood, plus the penalty) singular, or not positive definite, on the
## coefficients that are estimated -- all but the inert "(fixed)" ones of parameters held
## by `fixed`? Returns, for .finalise() and .warn_singular_hessian():
## - rank: the numerical rank, the number of eigenvalues above tol times the largest
##   eigenvalue of H0, with each coefficient in the units of .coef_scale(), and less any
##   direction to which Vp gives zero or negative variance (below); NA if H is not
##   finite. tol = sqrt(eps) ~ 1.5e-8 is the relative tolerance at which MASS::ginv() and
##   ginv.egpd() treat a singular value as zero, and that of .edf_clamped()'s
##   semi-definiteness test; it is about the optimisers' relative function tolerance
##   (fntol 1e-8), below which a fit cannot resolve the curvature, and a direction this
##   flat has a variance at least 1/tol ~ 7e7 times the smallest. At a = k the generalized
##   Waring's outer-product Hessian is exactly singular; at the returned estimates (a and k
##   equal to 1e-6 to 3e-5 relative) the ratio was 1e-16 to 3e-12, while 95% of the
##   converged fits of the test suite with no predictor beyond |eta| = 9 are above 2e-5
##   (median 4e-3; below 1e-6 only fits approaching a plateau and other singular ones,
##   see the round-6 notes). The reference is the largest
##   eigenvalue of H0, not of H, so that a large smoothing parameter (a large penalty
##   eigenvalue) does not make the coefficients that the data determine look singular;
##   the penalty can only raise H's eigenvalues above H0's.
## - npar: the number of estimated coefficients.
## - free, Hf, Hs, cut, thr: the estimated coefficients, H on them unscaled and in the
##   units of .coef_scale(), pinv()'s cut-off (Armadillo's: the matrix size times its
##   largest |eigenvalue| times eps, for H with a unit diagonal on its zero rows, as
##   .VpVc() inverts it; eigenvalues of smaller magnitude get zero variance, negative ones
##   of larger magnitude negative variance), and the rank threshold.
.hessian_check <- function(H, H0, likdata, tol = sqrt(.Machine$double.eps)) {
  pinned <- vapply(likdata$X, .is_pinned_block, logical(1))
  free <- which(!(likdata$idpars %in% which(pinned)))
  out <- list(rank = NA_integer_, npar = length(free), free = free)
  if (!length(free))
    return(out)
  Hf <- H[free, free, drop = FALSE]
  H0f <- H0[free, free, drop = FALSE]
  if (!all(is.finite(Hf)) || !all(is.finite(H0f)))
    return(out)
  ev <- eigen(Hf, symmetric = TRUE, only.values = TRUE)$values
  cut <- nrow(H) * max(abs(ev), if (any(rowSums(abs(H)) == 0)) 1) * .Machine$double.eps
  d <- .coef_scale(likdata)[free]
  Hs <- Hf / outer(d, d)
  evs <- eigen(Hs, symmetric = TRUE, only.values = TRUE)$values
  ref <- if (identical(Hf, H0f)) evs[1L] else
    eigen(H0f / outer(d, d), symmetric = TRUE, only.values = TRUE)$values[1L]
  if (!is.finite(ref) || ref <= 0)
    ref <- max(abs(evs))
  thr <- tol * ref
  out$rank <- as.integer(min(sum(evs > thr), sum(ev >= cut)))
  out[c("Hf", "Hs", "cut", "thr")] <- list(Hf, Hs, cut, thr)
  out
}

## Warn, once, with condition class "egpd_singular_hessian", when the Hessian behind the
## standard errors (.hessian_check) is singular, or not positive definite, in a way that
## neither the standard errors nor `plateau` show. Among the estimated coefficients of the
## predictors that are not on a plateau of their link (`plateau`, .plateau_flags: those
## are reported there, and by print(), and are not warned about, by the round-5 decision):
## (a) Vp gives some direction zero or negative variance -- an eigenvalue of H at
##     pinv()'s cut-off or below, as at an exact singularity (e.g. collinear terms, or a
##     parameter whose score is zero, beyond the range in which a count family's
##     likelihood is evaluated) or where H is not positive definite: the standard errors
##     there are understated, or not defined; or
## (b) the family's Hessian is the outer product of the scores (BHHH), and it is
##     numerically singular (rank threshold of .hessian_check). The outer product is
##     singular where the scores are linearly dependent, which need not make the
##     likelihood flat: at the generalized Waring's a = k it is exactly singular while the
##     observed information is not, so its very large (or, once pinv() cuts it off, zero)
##     variance does not describe the likelihood.
## A merely near-singular exact (observed-information) Hessian is not warned about: its very
## large variances are the likelihood's own flatness, e.g. a parameter approaching a
## plateau. Such fits, and those with a predictor on a plateau, still get `hessian.rank`
## below the number of coefficients and a note in print() and summary(). Only converged
## fits are assessed: the non-convergence warning already says that the estimates are
## unreliable.
.warn_singular_hessian <- function(hc, fit, family) {
  if (is.na(hc$rank) || hc$rank >= hc$npar)
    return(invisible(FALSE))
  family <- as.vector(family)
  bhhh <- family %in% .egpd_bhhh_families && is.null(fit$likdata$observed.step)
  plat <- fit$plateau
  K <- which(!(fit$likdata$idpars[hc$free] %in% which(!is.na(plat) & plat)))
  if (!length(K))
    return(invisible(FALSE))
  e <- eigen(hc$Hf[K, K, drop = FALSE], symmetric = TRUE)
  zero <- abs(e$values) < hc$cut
  neg <- e$values <= -hc$cut
  ## each coefficient's share of the offending directions, to name them
  w <- rowSums(e$vectors[, zero | neg, drop = FALSE]^2)
  near <- FALSE
  if (bhhh) {
    es <- eigen(hc$Hs[K, K, drop = FALSE], symmetric = TRUE)
    near <- es$values <= hc$thr
    w <- pmax(w, rowSums(es$vectors[, near, drop = FALSE]^2))
  }
  if (!any(zero, neg, near))
    return(invisible(FALSE))
  o <- order(w, decreasing = TRUE)
  o <- o[w[o] >= 0.1 * w[o[1L]]]
  who <- names(fit$coefficients)[hc$free[K]][o[seq_len(min(4L, length(o)))]]
  pen <- paste(c(if (isTRUE(fit$likdata$adjust > 0)) "sandwich-adjusted ",
                  if (!is.null(fit$Sdata)) "penalised "), collapse = "")
  msg <- if (any(neg))
    sprintf(paste0("egpd(): the Hessian of the %snegative log-likelihood is not positive ",
                   "definite at the estimates (%d negative eigenvalue(s); numerical rank %d of ",
                   "%d estimated coefficients), so they may not be at a maximum."),
            pen, sum(neg), hc$rank, hc$npar)
  else
    sprintf(paste0("egpd(): the Hessian of the %snegative log-likelihood is singular at the ",
                   "estimates: numerical rank %d of %d estimated coefficients."),
            pen, hc$rank, hc$npar)
  if (bhhh)
    msg <- paste0(msg, sprintf(paste0(
      " For family '%s' it is the outer product of the per-observation scores (BHHH), ",
      "which is singular when the scores are linearly dependent or vanish -- e.g. the ",
      "generalized Waring at a = k (a = mu (rho - 1) / k), whose likelihood is ",
      "symmetric in a and k, or a parameter at a limit of its range."), family))
  msg <- paste0(msg, sprintf(paste0(
    " Standard errors, vcov(), confint(), predict(se.fit = TRUE) and simulate() do not ",
    "reflect the uncertainty along the unidentified direction(s)%s, mainly in %s. ",
    "Consider fixing a parameter (`fixed`), simplifying the model, or profile ",
    "likelihood intervals (confint(method = \"profile\"))."),
    if (any(neg)) " (Vp gives them negative variance)"
    else if (any(zero)) " (Vp gives them zero variance)" else "",
    paste(who, collapse = ", ")))
  warning(structure(class = c("egpd_singular_hessian", "warning", "condition"),
                    list(message = msg, call = NULL)))
  invisible(TRUE)
}

############ .edf ##########################

## Effective degrees of freedom per coefficient, diag((H0 + S)^-1 H0), with H0 the
## Hessian of the negative log-likelihood and S the penalty (both unscaled by gamma, as
## Vp). The formula assumes H0 is positive semi-definite. At a valid penalised optimum
## H0 + S is positive definite, but H0 need not be: with smooths on weakly identified
## parameters it had negative eigenvalues at all twelve of twelve simulated zero-truncated
## DEGPD optima, and the trace, unbounded below, gave edf -51.5 (AIC below -2 logLik).
## There the edf follow mgcv's rule (.edf_clamped). Vp and Vc are not touched: they stay
## the inverse of the exact penalised Hessian.
.edf <- function(beta, likfns, likdata, VpVc, sandwich) {
edf <- diag(crossprod(VpVc$Vp, VpVc$H0))
clamped <- .edf_clamped(VpVc$H0, likdata)
if (is.null(clamped)) edf else clamped
}

## The edf by mgcv's rule (the clamp loop of mgcv:::gam.fit5.post.proc) for a penalised
## fit whose H0 is not positive semi-definite, or NULL otherwise. The eigenvalues of H0
## below tol * (its largest eigenvalue) are raised to that level, with tol = 0 first and
## then 1e-7, 1.1e-6, 1.11e-5, ... until the clamped H0 plus S has full rank (pivoted
## Cholesky), and the edf are diag((Hc + S)^-1 Hc) for the clamped Hc. (mgcv clamps the
## Hessian after a diagonal rescaling; this clamps H0 itself.) NULL -- the plain formula,
## bit for bit -- when the fit has no penalty, when H0 passes chol(), or when its smallest
## eigenvalue is at least -tau times its largest (`tau` = sqrt(eps), the tolerance at which
## MASS::ginv() and ginv.egpd() count an eigenvalue as zero): an H0 that is positive
## semi-definite up to rounding, e.g. singular in an unidentified direction that the plain
## formula (through pinv() in .VpVc) gives no degree of freedom. Without a penalty the plain
## formula already gives each coefficient one degree of freedom. The inert coefficients of
## parameters held by `fixed` (zero rows of H0 and S, see .VpVc) are left out of the
## eigenproblem, where clamping their zero eigenvalues up would give them edf 1, and keep
## edf 0.
.edf_clamped <- function(H0, likdata, tau = sqrt(.Machine$double.eps)) {
  if (is.null(likdata$S) || !any(likdata$S != 0) || !all(is.finite(H0)))
    return(NULL)
  g <- if (is.null(likdata$gamma)) 1 else likdata$gamma
  S <- g * likdata$S
  free <- which(rowSums(abs(H0)) > 0 | rowSums(abs(S)) > 0)
  if (!length(free))
    return(NULL)
  H0f <- H0[free, free, drop = FALSE]
  if (!inherits(try(chol(H0f), silent = TRUE), "try-error"))
    return(NULL)
  eh <- eigen(H0f, symmetric = TRUE)
  mev <- max(eh$values)
  if (!is.finite(mev) || mev <= 0 || min(eh$values) >= -tau * mev)
    return(NULL)
  Sf <- S[free, free, drop = FALSE]
  tol <- 0
  dtol <- 1e-7
  repeat {
    Hc <- eh$vectors %*% (pmax(eh$values, tol * mev) * t(eh$vectors))
    L <- suppressWarnings(chol(Hc + Sf, pivot = TRUE))
    if (attr(L, "rank") == length(free) || tol > 1)
      break
    tol <- tol + dtol
    dtol <- 10 * dtol
  }
  if (attr(L, "rank") < length(free))
    return(NULL)
  ip <- order(attr(L, "pivot"))
  V <- chol2inv(L)[ip, ip, drop = FALSE]
  edf <- numeric(ncol(H0))
  edf[free] <- rowSums(V * Hc)
  edf
}

############ .swap ##########################

## The smoothing parameters of the penalty on the unscaled log-likelihood, as mgcv
## reports them: the fit with log smoothing parameters `lsp` on the likelihood divided
## by gamma has the penalty gamma * exp(lsp) on the likelihood itself. Supplying these
## as `sp` reproduces the fit whatever gamma is.
.reported_sp <- function(lsp, likdata) {
  g <- if (is.null(likdata$gamma)) 1 else likdata$gamma
  exp(as.vector(lsp)) * g
}

.swap <- function(fitreml, gams, likdata, VpVc, edf, smooths) {
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (smooths) sp <- .reported_sp(fitreml$par, likdata)
start <- 0L
for (i in seq_along(gams)) {
  idi <- likdata$idpars == i
  gams[[i]]$coefficients <- fitreml$beta[idi]
  names(gams[[i]]$coefficients) <- gams[[i]]$term.names
  gams[[i]]$Vp <- Vp[idi, idi, drop = FALSE]
  gams[[i]]$Vc <- Vc[idi, idi, drop = FALSE]
  if (smooths) {
    g <- gams[[i]]
    g$sp[] <- sp[start + seq_along(g$sp)]
    start <- start + length(g$sp)
    g$full.sp <- exp(g$lsp0) * c(1, g$sp)[1L + g$sp.index]
    gams[[i]] <- g
  }
  gams[[i]]$edf <- edf[idi]
}
gams
}

############ .finalise ##########################

## Non-convergence must never pass silently: the returned estimates are used by
## predict(), simulate() and predictive_coverage() without further checks.
.warn_nonconvergence <- function(code, report, family, outer = FALSE) {
  reason <- if (nzchar(report)) report else paste("convergence code", code)
  hint <- paste0("Try other starting values via `inits`",
                 if (family %in% c("egpd", "degpd", "zidegpd")) " (or `restarts = TRUE`)" else "",
                 if (identical(as.integer(code), 1L))
                   sprintf(", or raise control$%s$itlim", if (outer) "outer" else "inner")
                 else "", ".")
  msg <- sprintf(paste0("egpd() did not converge (%s; convergence = %d): logLik, AIC ",
                        "and BIC are not reported and the estimates are unreliable. %s"),
                 reason, as.integer(code), hint)
  warning(structure(class = c("egpd_nonconvergence", "warning", "condition"),
                    list(message = msg, call = NULL)))
}

## The distinct values of a smooth's covariates (its `term`s) in `data`: the columns
## they name, or for an expression such as log(x) or factor(g), its values, computed in
## `data` and the environment `env` of the smooth's formula, as for the set-up.
.smooth_data <- function(term, env, data) {
  e <- lapply(term, function(v) tryCatch(str2lang(v), error = function(err) NULL))
  col <- vapply(seq_along(term), function(j)
    if (is.symbol(e[[j]])) as.character(e[[j]]) else term[j], character(1))
  if (all(vapply(e, is.symbol, logical(1))) && all(col %in% names(data)))
    return(unique(data[, col, drop = FALSE]))
  out <- data[, integer(0), drop = FALSE]
  for (j in seq_along(term))
    out[[col[j]]] <- if (is.call(e[[j]])) eval(e[[j]], data, env) else data[[col[j]]]
  unique(out)
}

.finalise <- function(gams, data, likfns, likdata, Sdata, fitreml, VpVc, family, gotsmooth,
formula, responsenm, removeData, edf, linkNames) {
names(gams) <- linkNames
has.penalty <- !is.null(Sdata)
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (has.penalty) {
  gams$full.sp <- unlist(lapply(gams, "[[", "full.sp"))
  gams$sp <- .reported_sp(fitreml$par, likdata)
}
gams$nobs <- likdata$nobs
gams$convergence <- if (is.null(fitreml$convergence)) NA_integer_ else as.integer(fitreml$convergence)
report <- if (is.null(fitreml$report)) "" else paste(fitreml$report, collapse = " ")
## The log-likelihood of the data at the fitted coefficients is likfns$d0 itself, not
## .nllh.nopen: the latter is the objective being minimised, which a sandwich
## "magnitude" adjustment or gamma != 1 multiplies by k = 1/gamma (and a "curvature"
## adjustment evaluates at transformed coefficients). Reporting that scaled value
## made AIC/BIC favour adjusted fits for no reason.
ll <- -as.numeric(likfns$d0(as.vector(fitreml$beta), likdata))
## Objectives map invalid parameter values to 1e20, so a fit that never left an
## invalid start "converges" with logLik -1e20. That is a failure, not a fit, and
## the most basic one, so it is reported in preference to the optimiser's own code.
if (length(ll) != 1L || !is.finite(ll) || ll <= -1e19) {
  gams$convergence <- 2L
  report <- "the log-likelihood is not finite at the returned coefficients"
} else if (family %in% .egpd_count_families &&
           ll > sqrt(.Machine$double.eps) * max(1, likdata$nobs)) {
  ## A probability mass function is at most one, so a positive log-likelihood for a
  ## count family means its numerics broke down where the optimiser stopped (e.g. a
  ## truncation term computed by cancellation as f(0) -> 1, or a parameter driven
  ## far beyond the range the family's code can evaluate).
  gams$convergence <- 2L
  report <- sprintf(paste0("the log-likelihood (%.6g) is positive, which is impossible ",
                           "for a discrete family: the likelihood is numerically ",
                           "unreliable at the returned coefficients"), ll)
} else if (isTRUE(gams$convergence == 0L) &&
           !is.null(wall <- .count_wall_report(as.vector(fitreml$beta), likdata, family))) {
  gams$convergence <- 3L
  report <- paste0(report, "; but ", wall)
}
gams$convergence.report <- report
gams$logLik <- NA_real_
## NB: optimizers return a *double* 0 on success, so the old
## identical(gams$convergence, 0L) test was always FALSE and logLik/AIC/BIC
## were never populated. Compare numerically instead.
converged <- isTRUE(gams$convergence == 0L)
gams$df <- sum(edf)
if (converged) {
## The unpenalised log-likelihood at the fitted coefficients, and nothing else.
## This used to have likdata$const = .5 * Mp * log(2 * pi) subtracted from it,
## which is the REML criterion's normalising constant and has no place in a
## log-likelihood: it is not a Laplace correction (that would *add* the term and
## carry .5*log|S|+ - .5*log|H| with it), and .nllh.nopen is already used bare
## as the log-likelihood by the profile intervals in confint(). Because Mp counts
## one per linear predictor plus each smooth's null-space dimension, the offset
## grew with model complexity, so it biased AIC/BIC against the more flexible
## model rather than cancelling in a comparison.
gams$logLik <- ll
gams$AIC <- -2 * gams$logLik + 2 * gams$df
gams$BIC <- -2 * gams$logLik + log(likdata$nobs) * gams$df
if (!isTRUE(likdata$k == 1) || isTRUE(likdata$adjust > 0))
  gams$logLik.adjusted <- -.nllh.nopen(fitreml$beta, likdata, likfns)
} else {
gams$AIC <- gams$BIC <- 1e20
.warn_nonconvergence(gams$convergence, report, family,
                     outer = has.penalty && !identical(likdata$outer, "fixed"))
}
attr(gams, "df") <- sum(edf)
gams$simulate <- list(mu=fitreml$beta, Sigma=Vp)
gams$family <- family
gams$idpars <- likdata$idpars
nms <- names(gams)[seq_along(formula)]
if (family != "custom") {
  ## strip one leading link prefix only: gsub() also removed the trailing "log" of
  ## CompPareto's "logsdlog", naming it "sd"
  nms <- sub("^(logit|probit|log)", "", nms)
}
gams$predictor.names <- attr(formula, "predictor.names")
formula <- attr(formula, "stripped")
names(formula) <- nms
gams$call <- formula
gams$response.name <- responsenm
gams$gotsmooth <- gotsmooth
if (!removeData) {
    gams$data <- data
}
gams$Vc <- Vc
gams$Vp <- Vp
gams$Vlsp <- VpVc$Vlsp
## Retain the joint, unpenalised information only for fits needing a boundary
## test. Inverting Vp and subtracting a large penalty can lose this information.
if (any(vapply(Sdata, .boundary_smooth, logical(1))) && !likdata$sandwich)
  gams$re.info <- VpVc$H0
## The REML criterion at the returned smoothing parameters; without smooths, the
## minimised negative log-likelihood (scaled by k = 1/gamma or the sandwich factor).
gams$negREML <- if (!is.null(fitreml$negREML)) fitreml$negREML else
  if (is.null(fitreml$objective)) NA_real_ else as.numeric(fitreml$objective)
gams$coefficients <- as.vector(fitreml$beta)
for (i in seq_along(likdata$X)) {
  P <- gams[[i]][["P"]]
  if (!is.null(P) && rcond(P) <= sqrt(.Machine$double.eps) &&
      nrow(gams[[i]]$model) != nrow(data)) {
    ## A basis subsample does not contain all the rows needed for smooth
    ## tests. A singular P cannot reconstruct the prediction design from X.
    g <- gams[[i]]
    class(g) <- "gam"
    gams[[i]]$summary.X <- .gam_predict_rows(g, data, "lpmatrix")
  }
  gams[[i]]$X <- likdata$X[[i]]
  if (likdata$duplicate == 1) {
    gams[[i]]$X <- gams[[i]]$X[likdata$dupid + 1, , drop = FALSE]
    if (!is.null(rownames(gams[[i]]$X))) rownames(gams[[i]]$X) <- rownames(data)
  }
  gams[[i]]$fitted <- as.vector(gams[[i]]$X %*% gams[[i]]$coefficients)
  if (length(likdata$offsets[[i]]) > 0)
    gams[[i]]$fitted <- gams[[i]]$fitted + likdata$offsets[[i]]
  names(gams[[i]]$coefficients) <- colnames(gams[[i]]$X)
}
gams$likdata <- likdata
gams$likfns <- likfns
gams$plateau <- .plateau_flags(gams, family)
if (has.penalty) gams$Sdata <- Sdata
gams$formula <- formula
gams$compacted <- likdata$duplicate == 1
if (gams$compacted) gams$compactid <- likdata$dupid + 1
## The covariate values of each smooth, once per distinct set of covariates. A smooth
## of an expression, such as s(log(x)), te(log(x), z) or the factor(g) of an "fs"
## smooth, has the expression's values, evaluated as by the set-up: selecting columns
## by the terms' text failed for those with "undefined columns selected", after the
## model had been fitted.
smooth.terms <- unlist(lapply(gotsmooth, function(i) lapply(gams[[i]]$smooth, function(z) z$term)),
                       recursive = FALSE)
smooth.envs <- unlist(lapply(gotsmooth, function(i)
  rep(list(environment(formula[[i]])), length(gams[[i]]$smooth))), recursive = FALSE)
first <- which(!duplicated(smooth.terms))
gams$plotdata <- lapply(first, function(k) .smooth_data(smooth.terms[[k]], smooth.envs[[k]], data))
names(gams$coefficients) <- unlist(lapply(seq_along(likdata$X), function(i) paste(names(gams)[i], names(gams[[i]]$coefficients), sep = "_")))
gams$ngam <- length(formula)
## The numerical rank of the Hessian behind Vp and Vc on the estimated coefficients
## (.hessian_check), and a warning where the standard errors cannot show that it is
## singular (.warn_singular_hessian). It needs `plateau` and the coefficient names.
gams$hessian.rank <- c(rank = VpVc$hessian$rank, npar = VpVc$hessian$npar)
if (converged)
  .warn_singular_hessian(VpVc$hessian, gams, family)
class(gams) <- "egpd"
return(gams)
}

############ .setup.fixed.egpd ##########################

## Map a response pin back through a custom unlink. A supplied forward link
## also lets the caller select the valid branch/domain. Never accept a trial
## eta without checking it against the very function prediction will use.
.custom_fixed_eta <- function(fn, v, nm) {
  scalar <- function(f, x) {
    if (!is.function(f)) return(NA_real_)
    y <- tryCatch(suppressWarnings(f(x)), error = function(e) NA_real_)
    if (is.numeric(y) && length(y) == 1L && is.finite(y)) as.numeric(y) else NA_real_
  }
  link <- attr(fn, "link", exact = TRUE)
  if (is.null(link)) {
    if (identical(fn, exp)) link <- log
    if (identical(fn, identity)) link <- identity
    if (identical(fn, plogis)) link <- qlogis
    if (identical(fn, pnorm)) link <- qnorm
  }
  eta <- NA_real_
  if (!is.null(link)) {
    eta <- scalar(link, v)
  } else {
    ## Bracket an interior value in either orientation, without assuming
    ## that zero (or either infinity) is in the domain. Keep non-finite
    ## probes as breaks. A flat tail rounded to v is not a strict bracket.
    steps <- 2^(-20:20)
    x <- c(-rev(steps), 0, steps)
    residual <- function(e) scalar(fn, e) - v
    y <- vapply(x, residual, numeric(1))
    keep <- is.na(y) | y != 0
    x <- x[keep]; y <- y[keep]
    i <- which(abs(diff(sign(y))) == 2)
    if (length(i) == 1L)
      eta <- tryCatch(uniroot(residual, x[i + 0:1], f.lower = y[i], f.upper = y[i+1L],
                              tol = .Machine$double.eps, check.conv = TRUE)$root,
                      error = function(e) NA_real_)
  }
  if (!is.finite(eta) || !isTRUE(abs(scalar(fn, eta) - v) <= 1e-12 * abs(v)))
    stop("Cannot map fixed$", nm, " through its custom unlink to a finite link-scale value. ",
         "Check the range or supply the forward link as the unlink function's 'link' attribute.",
         call. = FALSE)
  eta
}

## Sugar for holding distribution parameters at known values.
##
## `fixed` is a named list of RESPONSE-scale values, e.g. list(kappa = 1). Each is
## translated into the offset-only predictor `~ offset(<col>) - 1`, which pins that
## parameter without any free coefficient (see .predictable.gam). The translation
## exists because the offset has to be supplied on the LINK scale, which is the part
## users get wrong: log-linked parameters (lsigma, lkappa, lxi, ...) take log(value),
## logit-linked ones (logitpi, ...) qlogis(value), and identity-linked ones (e.g. `xi`
## under degpd.args link = "identity") are pinned on their natural scale.
##
## Names are matched by .match_param_name() exactly as in a named formula list: the
## natural name (kappa, xi, pi, lambda, sdlog, ...), the internal and long
## link-prefixed names (lkappa, logkappa, logitpi, ...) and link-prefixed forms of
## the natural name all refer to the same parameter, and the value is always on the
## RESPONSE scale whichever spelling is used, so fixed = list(kappa = 1) and
## list(lkappa = 1) both mean "kappa = 1".
.setup.fixed.egpd <- function(fixed, formula, data, nms, xi.max=Inf, ptab=NULL, unlink=NULL) {
  if (is.null(names(fixed)) || any(!nzchar(names(fixed))))
    stop("`fixed` must be a named list, e.g. fixed = list(kappa = 1)", call. = FALSE)
  if (is.null(nms))
    stop("this family does not expose named parameters, so `fixed` cannot be used",
         call. = FALSE)
  if (!is.list(formula))
    stop("`fixed` requires the multi-parameter list form of `formula`", call. = FALSE)
  if (is.null(ptab)) {
    ## Direct callers may supply a built-in family's short internal names.
    nat <- unname(.natural_param_names[nms])
    nat[is.na(nat)] <- .strip_link(nms[is.na(nat)])
    ptab <- .param_table(nms, nat = nat)
  }

  n    <- nrow(data)
  offsets <- list(); values <- list()
  done <- integer(0)
  vars <- unlist(lapply(formula, all.vars))

  for (nm in names(fixed)) {
    v <- fixed[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
      stop(sprintf("fixed$%s must be a single finite number", nm), call. = FALSE)

    i <- .match_param_name(nm, ptab)
    if (i %in% done)
      stop(sprintf("`fixed` gives parameter '%s' more than once", ptab$nat[i]),
           call. = FALSE)
    done <- c(done, i)

    if (i <= length(unlink) && !is.null(unlink[[i]])) {
      eta <- .custom_fixed_eta(unlink[[i]], v, nm)
    } else {
      link <- ptab$links[i]
      if (link %in% c("logit", "probit", "cloglog") && (v <= 0 || v >= 1))
        stop(sprintf("fixed$%s must be strictly between zero and one", nm), call. = FALSE)
      if (link == "log" && v <= 0)
        stop(sprintf("fixed$%s must be positive (it is fitted on the log scale)", nm),
             call. = FALSE)
      eta <- switch(link, log = log(v), logit = qlogis(v), probit = qnorm(v),
                    cloglog = log(-log1p(-v)), identity = v)
      if (nms[i] == "lxi" && is.finite(xi.max)) {
        if (v >= xi.max)
          stop("Fixed xi must be strictly below xi.max.", call.=FALSE)
        eta <- log(-xi.max * log1p(-v / xi.max))
      }
    }

    col <- .aux_response_name(paste0(".fixedpar_", nms[i]), c(names(data), vars))
    data[[col]] <- rep(eta, n)
    ## the response as written, not its text (see .reorder_named_formulae)
    response <- if (i == 1L) formula[[1]][[2L]] else NULL
    formula[[i]] <- .reformulate_env(paste0("offset(", col, ")"), response=response,
                                     intercept=FALSE, like=formula[[i]])
    offsets[[col]] <- eta            # link scale, for predict() on newdata
    values[[nms[i]]] <- v            # response scale, for the user
  }

  list(formula = formula, data = data, offsets = offsets, values = values)
}
