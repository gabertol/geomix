#' Grid search for optimal K in unmixing models
#'
#' Runs `unmix()` for a sequence of K values and multiple random seeds,
#' summarizing model fit statistics (loss, AIC, BIC). Intended for both
#' interactive exploration (e.g., in Shiny) and scripted benchmarking.
#'
#' @param data_list Named list of matrices (samples × features) passed to `unmix()`.
#' @param data_types Named character vector indicating type of each dataset:
#'   `"continuous"`, `"compositional"`, or `"clr"`.
#' @param K_seq Integer vector of K values to test.
#' @param reps Number of replicate runs per K (different seeds).
#' @param max_iter Maximum iterations passed to `unmix()`.
#' @param tol Convergence tolerance passed to `unmix()`.
#' @param seeds Optional integer vector of length = `reps`. Defaults to `1:reps`.
#' @param verbose Logical; print progress messages.
#'
#' @return A tibble with one row per run containing:
#'   K, replicate number, seed, loss, iterations, convergence flag,
#'   number of samples, and dimensions of each block (as JSON).
#' @export
k_optimize <- function(data_list, data_types, K_seq, reps = 8,
                       max_iter = 1000, tol = 1e-8, seeds = NULL, verbose = FALSE) {
  stopifnot(length(data_list) > 0, length(data_types) > 0)
  if (is.null(seeds)) seeds <- seq_len(reps)
  if (length(seeds) != reps) stop("length(seeds) must equal reps")

  dims <- .geomix_dims(data_list)

  out <- lapply(K_seq, function(K) {
    lapply(seq_len(reps), function(r) {
      set.seed(seeds[r])
      fit <- try(unmix(data_list, data_types, K = K, max_iter = max_iter, tol = tol, verbose = FALSE),
                 silent = TRUE)
      if (inherits(fit, "try-error") || !is.list(fit)) {
        return(tibble::tibble(
          K = K, rep = r, seed = seeds[r],
          loss = NA_real_, iters = NA_integer_, converged = FALSE,
          n_samples = dims$n_samples, B_dims_json = dims$B_dims_json
        ))
      }
      tibble::tibble(
        K = K, rep = r, seed = seeds[r],
        loss = as.numeric(fit$loss),
        iters = length(fit$loss_history),
        converged = isTRUE(fit$converged),
        n_samples = dims$n_samples,
        B_dims_json = dims$B_dims_json
      )
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()

  out
}

#' Summarize K optimization results
#'
#' Computes mean and standard deviation of the loss (and AIC/BIC if
#' `data_list` is provided). Returns both the summary table and the
#' best K according to the selected metric.
#'
#' @param kgrid Output from [k_optimize()].
#' @param metric Criterion used for ranking: `"loss"`, `"aic"`, or `"bic"`.
#' @param data_list Optional; needed for AIC/BIC computation.
#'
#' @return A list with:
#' \item{summary}{A tibble summarized by K.}
#' \item{k_best}{The optimal K according to the chosen metric.}
#' @export
k_summary <- function(kgrid, metric = c("loss", "aic", "bic"), data_list = NULL) {
  metric <- match.arg(metric)
  if (nrow(kgrid) == 0) stop("Empty kgrid input.")

  sum_loss <- kgrid |>
    dplyr::group_by(K) |>
    dplyr::summarise(
      mean_loss = mean(loss, na.rm = TRUE),
      sd_loss   = stats::sd(loss, na.rm = TRUE),
      n_ok      = sum(is.finite(loss)),
      .groups = "drop"
    )

  if (metric %in% c("aic", "bic")) {
    if (is.null(data_list)) stop("data_list must be provided for AIC/BIC.")
    dims <- .geomix_dims(data_list)
    k_params <- dims$n_samples * sum_loss$K + sum_loss$K * dims$total_features
    n_obs    <- dims$n_samples * dims$total_features
    sum_loss$aic <- 2 * k_params + 2 * sum_loss$mean_loss
    sum_loss$bic <- k_params * log(n_obs) + 2 * sum_loss$mean_loss
  }

  order_idx <- switch(metric,
                      loss = order(sum_loss$mean_loss),
                      aic  = order(sum_loss$aic),
                      bic  = order(sum_loss$bic))
  k_best <- sum_loss$K[order_idx][1]
  list(summary = sum_loss, k_best = k_best)
}

#' Detect an "elbow" in the K–loss curve
#'
#' Simple geometric heuristic (distance to diagonal method).
#'
#' @param x Numeric vector of K values (increasing order).
#' @param y Numeric vector of loss values.
#' @return Integer K corresponding to the elbow (or `NA`).
#' @export
k_find_elbow <- function(x, y) {
  if (length(x) < 3) return(NA_integer_)
  xn <- (x - min(x)) / (max(x) - min(x) + 1e-12)
  yn <- (y - min(y, na.rm = TRUE)) / (max(y, na.rm = TRUE) - min(y, na.rm = TRUE) + 1e-12)
  d  <- abs(yn - xn) / sqrt(2)
  idx <- which.max(d)
  as.integer(x[idx])
}

#' Plot the K-curve using ggplot2
#'
#' Creates a base ggplot object (optionally shaded by ±1 SD) suitable for
#' wrapping with `plotly::ggplotly()` in Shiny apps.
#'
#' @param ksum Summary tibble returned by [k_summary()].
#' @param show_sd Logical; include shaded standard deviation area.
#' @param elbow_k Optional K value to highlight with a dashed line.
#' @return A `ggplot` object.
#' @export
k_plot <- function(ksum, show_sd = TRUE, elbow_k = NULL) {
  stopifnot(all(c("K", "mean_loss", "sd_loss") %in% names(ksum)))
  library(ggplot2)
  p <- ggplot(ksum, aes(K, mean_loss)) +
    geom_line() + geom_point() +
    labs(x = "K", y = "Mean loss", title = "K optimization curve")
  if (show_sd) {
    p <- p + geom_ribbon(aes(ymin = pmax(mean_loss - sd_loss, 0),
                             ymax = mean_loss + sd_loss), alpha = 0.2)
  }
  if (!is.null(elbow_k)) {
    p <- p + geom_vline(xintercept = elbow_k, linetype = 2)
  }
  p
}

# -------------------- internal utilities --------------------

# Compute dimensional metadata for AIC/BIC
.geomix_dims <- function(data_list) {
  stopifnot(is.list(data_list), length(data_list) >= 1)
  n_samples <- nrow(data_list[[1]])
  feats <- vapply(data_list, ncol, 1L)
  total_features <- sum(feats)
  B_dims <- lapply(data_list, dim)
  B_dims_json <- tryCatch(jsonlite::toJSON(B_dims, auto_unbox = TRUE), error = function(e) NA_character_)
  list(n_samples = n_samples, total_features = total_features, B_dims_json = B_dims_json)
}
