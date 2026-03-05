# =============================================================================
# CROSS-VALIDATION UNMIXING FOR INCOMPLETE DATASETS
# =============================================================================

#' Cross-Validation Unmixing for Incomplete Datasets
#'
#' Performs unmixing using cross-validation to handle datasets where
#' not all samples have all data types. Fits B matrices on samples
#' with complete data, then uses solve_A_given_B to estimate A for
#' samples with partial coverage.
#'
#' @param data_list Named list of data matrices. Matrices can have
#'   different samples (row names). All samples across all blocks
#'   will be included in the output.
#' @param data_types Named character vector specifying data types.
#' @param K Number of sources (end-members).
#' @param min_blocks_for_B Minimum number of blocks a sample must have
#'   to be used for fitting B matrices (default: 2).
#' @param weighting Block weighting method (default: "inv_f").
#' @param method Solver method: "nnls" (default) or "ols".
#' @param verbose Print progress? (default: TRUE)
#' @param ... Additional arguments passed to unmix().
#'
#' @return A list of class "unmix_cv_result" containing:
#'   \item{A}{Full mixing matrix (all samples x K)}
#'   \item{A_type}{Character vector indicating estimation type per sample:
#'     "complete" (used in B fitting), "partial" (estimated via CV)}
#'   \item{B_list}{Source signature matrices (from complete samples)}
#'   \item{sample_coverage}{Data frame showing which blocks each sample has}
#'   \item{fit_complete}{Original unmix fit on complete samples}
#'   \item{n_complete}{Number of samples with complete coverage}
#'   \item{n_partial}{Number of samples with partial coverage}
#'
#' @details
#' The algorithm works in two stages:
#'
#' **Stage 1: Fit B on complete samples**
#' - Identify samples present in at least `min_blocks_for_B` blocks
#' - Run standard unmix() on these samples to get B matrices
#'
#' **Stage 2: Estimate A for partial samples**
#' - For each sample with partial coverage, use solve_A_given_B()
#'   on available blocks
#' - Combine estimates using inverse-variance weighting or averaging
#'
#' This approach maximizes sample utilization: instead of discarding
#' samples missing one proxy, their mixing proportions are estimated
#' from available data using the learned source signatures.
#'
#' @examples
#' \dontrun{
#' # Suppose we have:
#' # - 50 samples with DZ data
#' # - 40 samples with PT data (30 overlap with DZ)
#' # - 35 samples with HM data (25 overlap with DZ)
#'
#' result <- unmix_cv(
#'   data_list = list(DZ = DZ_mat, PT = PT_mat, HM = HM_mat),
#'   data_types = c(DZ = "continuous", PT = "simplex", HM = "simplex"),
#'   K = 3,
#'   min_blocks_for_B = 2
#' )
#'
#' # All 50+ unique samples now have A estimates
#' nrow(result$A)
#'
#' # See coverage
#' result$sample_coverage
#'
#' # Which samples used CV?
#' table(result$A_type)
#' }
#'
#' @seealso \code{\link{unmix}}, \code{\link{solve_A_given_B}}
#' @export
unmix_cv <- function(data_list,
                     data_types = NULL,
                     K,
                     min_blocks_for_B = 2,
                     weighting = "inv_f",
                     method = "nnls",
                     verbose = TRUE,
                     ...) {

  n_blocks <- length(data_list)
  block_names <- names(data_list)

  if (is.null(block_names)) {
    block_names <- paste0("Block", seq_len(n_blocks))
    names(data_list) <- block_names
  }

  # Set default data types
  if (is.null(data_types)) {
    data_types <- setNames(rep("continuous", n_blocks), block_names)
  }

  # ==========================================================================
  # 1. BUILD SAMPLE COVERAGE MAP
  # ==========================================================================

  all_samples <- unique(unlist(lapply(data_list, rownames)))
  n_total <- length(all_samples)

  if (verbose) {
    cat("=== Cross-Validation Unmixing ===\n")
    cat(sprintf("Total unique samples: %d\n", n_total))
    cat(sprintf("Data blocks: %d (%s)\n", n_blocks, paste(block_names, collapse = ", ")))
  }

  # Coverage matrix: which samples have which blocks
  coverage <- matrix(FALSE, nrow = n_total, ncol = n_blocks,
                     dimnames = list(all_samples, block_names))

  for (b in block_names) {
    samples_in_block <- rownames(data_list[[b]])
    coverage[samples_in_block, b] <- TRUE
  }

  n_blocks_per_sample <- rowSums(coverage)

  # Classify samples
  complete_samples <- all_samples[n_blocks_per_sample == n_blocks]
  partial_samples <- all_samples[n_blocks_per_sample >= min_blocks_for_B &
                                   n_blocks_per_sample < n_blocks]
  sparse_samples <- all_samples[n_blocks_per_sample < min_blocks_for_B]

  # Samples to use for B fitting
  samples_for_B <- all_samples[n_blocks_per_sample >= min_blocks_for_B]

  if (verbose) {
    cat(sprintf("\nSample coverage:\n"))
    cat(sprintf("  Complete (%d blocks): %d samples\n", n_blocks, length(complete_samples)))
    cat(sprintf("  Partial (>=%d blocks): %d samples\n", min_blocks_for_B, length(partial_samples)))
    cat(sprintf("  Sparse (<%d blocks): %d samples\n", min_blocks_for_B, length(sparse_samples)))
    cat(sprintf("\nUsing %d samples for B estimation\n", length(samples_for_B)))
  }

  if (length(samples_for_B) < K) {
    stop(sprintf(
      "Not enough samples with >= %d blocks for B estimation. Have %d, need >= %d (K)",
      min_blocks_for_B, length(samples_for_B), K
    ))
  }

  # ==========================================================================
  # 2. FIT B ON SAMPLES WITH SUFFICIENT COVERAGE
  # ==========================================================================

  if (verbose) cat("\n--- Stage 1: Fitting B matrices ---\n")

  # Find common samples across blocks that meet threshold
  # Use the intersection of samples that have at least min_blocks_for_B

  # For B fitting, we need to use samples that are in ALL selected blocks

  # But we want to maximize coverage, so we iterate

  # Strategy: use samples present in at least 2 blocks together
  # Build data_list for fitting using available samples per block

  # Simpler approach: find largest subset with complete coverage
  # among samples_for_B

  # Get samples present in all blocks (complete samples)
  if (length(complete_samples) >= K) {
    # Ideal case: use complete samples
    samples_fit <- complete_samples
  } else {
    # Fallback: find blocks with most overlap and use those
    # Use samples with maximum coverage
    max_coverage <- max(n_blocks_per_sample)
    samples_fit <- all_samples[n_blocks_per_sample == max_coverage]

    if (length(samples_fit) < K) {
      # Last resort: use all samples that meet minimum
      samples_fit <- samples_for_B
    }
  }

  # Determine which blocks to use for fitting
  # (those that have all samples_fit)
  blocks_for_fit <- block_names[sapply(block_names, function(b) {
    all(samples_fit %in% rownames(data_list[[b]]))
  })]

  if (length(blocks_for_fit) == 0) {
    stop("No blocks contain all fitting samples. Try reducing min_blocks_for_B.")
  }

  if (verbose) {
    cat(sprintf("Fitting with %d samples across %d blocks\n",
                length(samples_fit), length(blocks_for_fit)))
  }

  # Prepare data for fitting
  data_list_fit <- lapply(blocks_for_fit, function(b) {
    data_list[[b]][samples_fit, , drop = FALSE]
  })
  names(data_list_fit) <- blocks_for_fit

  data_types_fit <- data_types[blocks_for_fit]

  # Run unmix on complete data
  fit_complete <- unmix(
    data_list = data_list_fit,
    data_types = data_types_fit,
    K = K,
    weighting = weighting,
    verbose = verbose,
    ...
  )

  # ==========================================================================
  # 3. ESTIMATE A FOR ALL SAMPLES
  # ==========================================================================

  if (verbose) cat("\n--- Stage 2: Estimating A for all samples ---\n")

  # Initialize A matrix for all samples
  A_full <- matrix(NA_real_, nrow = n_total, ncol = K,
                   dimnames = list(all_samples, colnames(fit_complete$A)))

  A_type <- character(n_total)
  names(A_type) <- all_samples

  # For samples used in fitting: use original A
  A_full[samples_fit, ] <- fit_complete$A
  A_type[samples_fit] <- "complete"

  # For remaining samples: estimate via CV
  remaining_samples <- setdiff(all_samples, samples_fit)

  if (length(remaining_samples) > 0 && verbose) {
    cat(sprintf("Estimating A for %d additional samples via CV...\n",
                length(remaining_samples)))
  }

  # Get B matrices for all blocks (including those not in fit)
  # For blocks not in fit, we need to estimate B
  B_list_full <- list()

  for (b in block_names) {
    if (b %in% blocks_for_fit) {
      B_list_full[[b]] <- fit_complete$B_list[[b]]
    } else {
      # Estimate B for this block using fitted A
      samples_both <- intersect(samples_fit, rownames(data_list[[b]]))
      if (length(samples_both) >= K) {
        A_sub <- fit_complete$A[samples_both, , drop = FALSE]
        X_sub <- data_list[[b]][samples_both, , drop = FALSE]
        B_list_full[[b]] <- nnls_B_given_A(A_sub, X_sub)
        rownames(B_list_full[[b]]) <- colnames(fit_complete$A)
      } else {
        if (verbose) {
          cat(sprintf("  Warning: Cannot estimate B for block %s (too few overlapping samples)\n", b))
        }
        B_list_full[[b]] <- NULL
      }
    }
  }

  # Estimate A for remaining samples
  for (samp in remaining_samples) {
    # Get blocks available for this sample
    blocks_avail <- block_names[coverage[samp, ]]
    blocks_avail <- blocks_avail[!sapply(B_list_full[blocks_avail], is.null)]

    if (length(blocks_avail) == 0) {
      A_type[samp] <- "no_data"
      next
    }

    # Estimate A from each available block
    A_estimates <- list()

    for (b in blocks_avail) {
      X_samp <- data_list[[b]][samp, , drop = FALSE]
      B_block <- B_list_full[[b]]

      A_est <- solve_A_given_B(X_samp, B_block, method = method, verbose = FALSE)
      A_estimates[[b]] <- as.vector(A_est)
    }

    # Combine estimates (simple average for now)
    # TODO: could use inverse-variance weighting based on reconstruction error
    A_combined <- Reduce(`+`, A_estimates) / length(A_estimates)
    A_combined <- A_combined / sum(A_combined)  # Re-normalize

    A_full[samp, ] <- A_combined
    A_type[samp] <- "partial"
  }

  # ==========================================================================
  # 4. BUILD RESULT
  # ==========================================================================

  # Sample coverage data frame
  coverage_df <- as.data.frame(coverage)
  coverage_df$sample <- rownames(coverage_df)
  coverage_df$n_blocks <- n_blocks_per_sample
  coverage_df$A_type <- A_type
  coverage_df <- coverage_df[, c("sample", "n_blocks", "A_type", block_names)]

  result <- list(
    A = A_full,
    A_type = A_type,
    B_list = B_list_full,
    sample_coverage = coverage_df,
    fit_complete = fit_complete,
    n_complete = sum(A_type == "complete"),
    n_partial = sum(A_type == "partial"),
    n_no_data = sum(A_type == "no_data"),
    K = K,
    block_names = block_names,
    call = match.call()
  )

  class(result) <- "unmix_cv_result"

  if (verbose) {
    cat("\n=== Summary ===\n")
    cat(sprintf("Total samples with A estimates: %d\n", sum(!is.na(A_full[,1]))))
    cat(sprintf("  From complete fit: %d\n", result$n_complete))
    cat(sprintf("  From CV (partial): %d\n", result$n_partial))
    if (result$n_no_data > 0) {
      cat(sprintf("  No data available: %d\n", result$n_no_data))
    }
  }

  result
}


#' Print Method for unmix_cv_result
#' @export
print.unmix_cv_result <- function(x, ...) {
  cat("Cross-Validation Unmixing Results\n")
  cat("==================================\n\n")
  cat("Sources (K):", x$K, "\n")
  cat("Data blocks:", paste(x$block_names, collapse = ", "), "\n\n")
  cat("Sample estimates:\n")
  cat(sprintf("  Complete (joint fit): %d\n", x$n_complete))
  cat(sprintf("  Partial (CV):         %d\n", x$n_partial))
  if (x$n_no_data > 0) {
    cat(sprintf("  No data:              %d\n", x$n_no_data))
  }
  cat(sprintf("\nTotal: %d samples\n", nrow(x$A)))
  invisible(x)
}


#' Summary Method for unmix_cv_result
#' @export
summary.unmix_cv_result <- function(object, ...) {
  print(object)
  cat("\nSample coverage:\n")
  print(table(object$A_type))
  cat("\nMixing matrix A (first 6 rows):\n")
  print(head(object$A))
  invisible(object)
}






