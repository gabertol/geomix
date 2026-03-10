#' Bootstrap Uncertainty Estimation for Source Unmixing
#'
#' Performs parametric bootstrap to estimate uncertainty in mixing
#' proportions (A matrix) and source signatures (B matrices).
#' For distributional (KDE) blocks, resamples individual grains per sample
#' and rebuilds KDE. For compositional blocks, resamples counts via
#' multinomial. This approach correctly propagates grain-level sampling
#' uncertainty and avoids the NMF label-switching problem via Hungarian
#' alignment at each iteration.
#'
#' @param data_list Named list of processed data matrices (output of
#'   build_kde_block / build_compositional_block — the $data_mat matrices)
#' @param data_types Character vector specifying data types per block
#'   ("continuous" or "simplex")
#' @param K Number of sources
#' @param raw_data_list Named list of raw data.frames (long format), one per
#'   KDE block. Each must have a "sample" column. Used for grain resampling.
#'   Blocks not present here fall back to row-resampling.
#' @param kde_params Named list of build_kde_block parameters per block.
#'   Each element: list(vars, n_points, bw_method, bw_multiplier).
#' @param n_boot Number of bootstrap iterations (default: 100)
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param parallel Logical, use parallel processing? (default: FALSE)
#' @param n_cores Number of cores if parallel=TRUE
#' @param verbose Logical, print progress? (default: TRUE)
#' @param fit_original Optional: existing unmix() fit to use as alignment
#'   reference and avoid re-fitting. If NULL, fit is computed internally.
#' @param ... Additional arguments passed to unmix()
#'
#' @return A list of class "bootstrap_result"
#' @export
bootstrap_unmix <- function(data_list,
                            data_types   = NULL,
                            K,
                            raw_data_list = NULL,
                            kde_params    = NULL,
                            n_boot        = 100,
                            conf_level    = 0.95,
                            parallel      = FALSE,
                            n_cores       = NULL,
                            verbose       = TRUE,
                            fit_original  = NULL,
                            ...) {

  if (verbose) message(sprintf("Starting parametric bootstrap with %d iterations...", n_boot))

  # ── Original fit ───────────────────────────────────────────────────────────
  if (is.null(fit_original)) {
    fit_orig <- unmix(
      data_list  = data_list,
      data_types = data_types,
      K          = K,
      seed       = 123,
      verbose    = FALSE,
      ...
    )
  } else {
    fit_orig <- fit_original
  }

  # Derive K from fit (avoids mismatch when fit_original is provided)
  K        <- ncol(fit_orig$A)
  n_blocks <- length(data_list)

  # Resolve sample names BEFORE N (N derived from samp_names length)
  samp_names <- rownames(data_list[[1]])
  if (is.null(samp_names) || length(samp_names) == 0)
    samp_names <- rownames(fit_orig$A)
  if ((is.null(samp_names) || length(samp_names) == 0) && !is.null(raw_data_list)) {
    raw1 <- as.data.frame(raw_data_list[[1]])
    if ("sample" %in% names(raw1))
      samp_names <- sort(unique(raw1[["sample"]]))
  }
  N <- if (!is.null(samp_names) && length(samp_names) > 0) length(samp_names) else fit_orig$N
  if (is.null(samp_names) || length(samp_names) == 0)
    samp_names <- paste0("S", seq_len(N))

  # ── Storage ────────────────────────────────────────────────────────────────
  A_boot <- array(NA, dim = c(n_boot, N, K),
                  dimnames = list(NULL, samp_names,
                                  colnames(fit_orig$A)))
  B_boot_list <- lapply(seq_along(data_list), function(i) {
    array(NA, dim = c(n_boot, K, ncol(fit_orig$B_list[[i]])))
  })
  names(B_boot_list) <- names(data_list)

  # ── Hungarian alignment helper ─────────────────────────────────────────────
  align_sources <- function(A_ref, A_new) {
    K_ <- ncol(A_ref)
    # Correlation-based cost: more robust than SSE for label switching
    cost <- matrix(0, K_, K_)
    for (i in seq_len(K_))
      for (j in seq_len(K_)) {
        r <- tryCatch(stats::cor(A_ref[, i], A_new[, j]),
                      error = function(e) 0)
        cost[i, j] <- 1 - ifelse(is.finite(r), r, 0)
      }
    perm <- if (requireNamespace("clue", quietly = TRUE)) {
      as.integer(clue::solve_LSAP(cost))
    } else {
      p <- integer(K_); rem <- seq_len(K_)
      for (i in seq_len(K_)) {
        best <- rem[which.min(cost[i, rem])]
        p[i] <- best; rem <- rem[rem != best]
      }
      p
    }
    perm
  }

  # ── Resample one block ─────────────────────────────────────────────────────
  resample_block <- function(nm, samp_names_) {
    dtype  <- if (!is.null(data_types)) data_types[[nm]] else "simplex"
    X_orig <- data_list[[nm]]

    if (dtype == "continuous" && !is.null(raw_data_list[[nm]])) {
      # ── Parametric bootstrap: resample grains, rebuild KDE ────────────────
      raw   <- as.data.frame(raw_data_list[[nm]])  # force data.frame, drop tibble
      kpar  <- kde_params[[nm]]
      vars_ <- if (!is.null(kpar$vars))          kpar$vars          else setdiff(names(raw), "sample")
      npts  <- if (!is.null(kpar$n_points))       kpar$n_points      else 129
      bwm   <- if (!is.null(kpar$bw_method))      kpar$bw_method     else "SJ"
      bwmul <- if (!is.null(kpar$bw_multiplier))  kpar$bw_multiplier else 1

      pieces <- lapply(samp_names_, function(s) {
        grains    <- raw[raw[["sample"]] == s, , drop = FALSE]
        resampled <- grains[sample(nrow(grains), nrow(grains), replace = TRUE), ]
        resampled[["sample"]] <- s
        resampled
      })
      boot_raw <- do.call(rbind.data.frame, pieces)
      rownames(boot_raw) <- NULL
      if (!"sample" %in% names(boot_raw))
        stop("Internal: sample column lost after rbind — pieces[[1]] names: ",
             paste(names(pieces[[1]]), collapse=","))
      blk <- build_kde_block(boot_raw, vars = vars_, n_points = npts,
                             bw_method = bwm, bw_multiplier = bwmul)
      return(blk$data_mat)

    } else if (dtype == "simplex") {
      # ── Parametric bootstrap: multinomial resample of counts ──────────────
      # X_orig rows are proportions; back-convert to counts using row sums
      # stored in the original block, or assume n=100 if unavailable
      n_counts <- round(rowSums(X_orig) * 100)
      n_counts[n_counts < 1] <- 100
      Xb <- t(sapply(seq_len(nrow(X_orig)), function(i) {
        p   <- X_orig[i, ]
        p   <- pmax(p, 0); p <- p / sum(p)
        cnt <- as.vector(stats::rmultinom(1, n_counts[i], p))
        cnt / sum(cnt)
      }))
      rownames(Xb) <- rownames(X_orig)
      colnames(Xb) <- colnames(X_orig)
      return(Xb)

    } else {
      # fallback: row-resample
      idx <- sample(N, N, replace = TRUE)
      Xb  <- X_orig[idx, , drop = FALSE]
      rownames(Xb) <- rownames(X_orig)
      return(Xb)
    }
  }

  # ── Bootstrap iterator ─────────────────────────────────────────────────────
  bootstrap_iter <- function(b) {
    data_boot <- lapply(names(data_list), resample_block, samp_names_ = samp_names)
    names(data_boot) <- names(data_list)

    # Multiple restarts — pick best alignment to reference
    n_restarts <- 3L
    best_fit  <- NULL
    best_cost <- Inf
    for (restart in seq_len(n_restarts)) {
      fb <- tryCatch({
        unmix(data_list  = data_boot,
              data_types = data_types,
              K          = K,
              seed       = b * 100L + restart,
              verbose    = FALSE, ...)
      }, error = function(e) NULL)
      if (is.null(fb)) next
      perm_r <- align_sources(fit_orig$A, fb$A)
      # Cost = total SSE after alignment
      A_aligned <- fb$A[, perm_r, drop = FALSE]
      cost_r    <- sum((fit_orig$A - A_aligned)^2)
      if (cost_r < best_cost) {
        best_cost <- cost_r
        best_fit  <- fb
        best_perm <- perm_r
      }
    }

    if (is.null(best_fit)) {
      if (verbose) message(sprintf("  iter %d: all restarts failed", b))
      return(NULL)
    }
    fit_boot <- best_fit

    # Align sources
    perm <- best_perm
    if (verbose && b <= 3)
      message(sprintf("  iter %d: perm = [%s]", b, paste(perm, collapse = ",")))

    list(
      A      = fit_boot$A[, perm, drop = FALSE],
      B_list = lapply(fit_boot$B_list, function(B) B[perm, , drop = FALSE])
    )
  }

  # ── Run ────────────────────────────────────────────────────────────────────
  if (parallel) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    if (verbose) message(sprintf("Using %d cores", n_cores))
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl,
                            c("data_list", "data_types", "K", "N", "samp_names",
                              "fit_orig", "raw_data_list", "kde_params", "unmix",
                              "build_kde_block", "align_sources", "resample_block"),
                            envir = environment())
    boot_results <- parallel::parLapply(cl, seq_len(n_boot), bootstrap_iter)
  } else {
    boot_results <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      if (verbose && b %% 10 == 0)
        message(sprintf("Bootstrap iteration %d/%d", b, n_boot))
      boot_results[[b]] <- bootstrap_iter(b)
    }
  }

  # ── Store results ──────────────────────────────────────────────────────────
  for (b in seq_along(boot_results)) {
    if (!is.null(boot_results[[b]])) {
      A_boot[b, , ] <- boot_results[[b]]$A
      for (i in seq_along(B_boot_list)) {
        tryCatch({
          B_boot_list[[i]][b, , ] <- boot_results[[b]]$B_list[[i]]
        }, error = function(e) NULL)  # skip if B dims differ (e.g. KDE grid mismatch)
      }
    }
  }

  valid <- !apply(A_boot, 1, function(x) all(is.na(x)))
  A_boot <- A_boot[valid, , , drop = FALSE]
  for (i in seq_along(B_boot_list))
    B_boot_list[[i]] <- B_boot_list[[i]][valid, , , drop = FALSE]
  n_valid <- sum(valid)

  if (verbose) message(sprintf("Completed %d/%d iterations", n_valid, n_boot))

  # ── Statistics ─────────────────────────────────────────────────────────────
  alpha   <- 1 - conf_level
  A_mean  <- apply(A_boot, c(2, 3), mean,     na.rm = TRUE)
  A_se    <- apply(A_boot, c(2, 3), sd,       na.rm = TRUE)
  A_lower <- apply(A_boot, c(2, 3), quantile, probs = alpha / 2,     na.rm = TRUE)
  A_upper <- apply(A_boot, c(2, 3), quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  rownames(A_mean) <- rownames(A_se) <-
    rownames(A_lower) <- rownames(A_upper) <- samp_names
  colnames(A_mean) <- colnames(A_se) <-
    colnames(A_lower) <- colnames(A_upper) <- colnames(fit_orig$A)

  result <- list(
    fit_original = fit_orig,
    A_boot       = A_boot,
    B_boot_list  = B_boot_list,
    A_mean       = A_mean,
    A_lower      = A_lower,
    A_upper      = A_upper,
    A_se         = A_se,
    conf_level   = conf_level,
    n_boot       = n_valid,
    call         = match.call()
  )
  class(result) <- "bootstrap_result"
  return(result)
}

#' @export
print.bootstrap_result <- function(x, ...) {
  cat("Bootstrap Unmixing Results\n==========================\n\n")
  cat("Iterations:", x$n_boot, "\n")
  cat("Confidence level:", x$conf_level * 100, "%\n")
  cat("Mean SE:", round(mean(x$A_se), 4), "\n")
  invisible(x)
}

#' @export
summary.bootstrap_result <- function(object, ...) {
  print(object)
  cat("\nA with CI (first 3 samples):\n")
  for (i in seq_len(min(3, nrow(object$A_mean)))) {
    cat(sprintf("\nSample %s:\n", rownames(object$A_mean)[i]))
    for (k in seq_len(ncol(object$A_mean))) {
      cat(sprintf("  Source%d: %.3f (CI: %.3f-%.3f, SE: %.3f)\n",
                  k, object$fit_original$A[i,k],
                  object$A_lower[i,k], object$A_upper[i,k], object$A_se[i,k]))
    }
  }
  invisible(object)
}
