#' Coupled Unmixing with Multiple Data Types
#'
#' Performs simultaneous unmixing of multiple data blocks (e.g., detrital zircon
#' ages, bulk petrology, heavy minerals) using non-negative least squares (NNLS)
#' with appropriate transformations for each data type.
#'
#' @param data_list Named list of data matrices. Each matrix should have samples
#'   as rows and features as columns. Names will be used to identify data blocks.
#' @param data_types Named character vector specifying the type of each data block.
#'   Options: "continuous" (default, e.g., KDE-discretized ages),
#'   "compositional" (e.g., mineral counts - will apply CLR transformation),
#'   "simplex" (already proportions, no transform),
#'   "clr" (already in CLR space).
#' @param K Integer, number of sources (end-members) to unmix.
#' @param weighting Global block weighting method to handle feature imbalance:
#'   "sqrt_inv_f" (default, weight = 1/sqrt(features)),
#'   "inv_f" (weight = 1/features, stronger correction),
#'   "none" (equal weights).
#' @param adaptive Sample-wise adaptive weighting:
#'   "none" (default, no adaptation),
#'   "entropy" (weight by estimation confidence per sample).
#' @param method Solver method: "nnls" (default) or "ols".
#' @param max_iter Maximum number of iterations (default: 2000).
#' @param tol Convergence tolerance for relative change in loss (default: 1e-8).
#' @param eta_A Step size for updating A (default: 0.5).
#' @param seed Random seed for reproducibility (default: 123).
#' @param verbose Logical, print progress? (default: FALSE).
#'
#' @return A list of class "unmix_result" containing:
#'   \item{A}{Mixing matrix (N x K), rows sum to 1}
#'   \item{B_list}{List of source signature matrices (K x F_i for each block)}
#'   \item{recon_list}{List of reconstructed data matrices}
#'   \item{data_types}{Character vector of data types used}
#'   \item{block_weights}{Global weights applied to each block}
#'   \item{sample_weights}{Sample-wise weights (if adaptive != "none")}
#'   \item{loss}{Final loss value}
#'   \item{loss_history}{Vector of loss at each iteration}
#'   \item{iters}{Number of iterations performed}
#'   \item{converged}{Logical, did algorithm converge?}
#'   \item{call}{The matched call}
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' N <- 50  # samples
#' K <- 3   # sources
#'
#' # Continuous data (e.g., detrital zircon ages)
#' DZ <- matrix(rnorm(N * 100), N, 100)
#'
#' # Compositional data (e.g., mineral counts)
#' BP_counts <- matrix(rpois(N * 10, lambda = 20), N, 10)
#'
#' # Unmix with default settings
#' result <- unmix(
#'   data_list = list(DZ = DZ, BP = BP_counts),
#'   data_types = c(DZ = "continuous", BP = "compositional"),
#'   K = 3
#' )
#'
#' # Unmix with strong feature correction
#' result <- unmix(
#'   data_list = list(DZ = DZ, BP = BP_counts),
#'   data_types = c(DZ = "continuous", BP = "simplex"),
#'   K = 3,
#'   weighting = "inv_f"
#' )
#'
#' # Unmix with adaptive weighting (for complementary blocks)
#' result <- unmix(
#'   data_list = list(DZ = DZ, BP = BP_counts),
#'   data_types = c(DZ = "continuous", BP = "simplex"),
#'   K = 3,
#'   weighting = "inv_f",
#'   adaptive = "entropy"
#' )
#'
#' # Examine results
#' print(result)
#' plot(result)
#' }
#'
#' @export
#' @importFrom compositions clr acomp
unmix <- function(data_list,
                  data_types = NULL,
                  K,
                  weighting = c("sqrt_inv_f", "inv_f", "none"),
                  adaptive = c("none", "entropy"),
                  method = c("nnls", "ols"),
                  max_iter = 2000,
                  tol = 1e-8,
                  eta_A = 0.5,
                  seed = 123,
                  verbose = FALSE) {

  # Match arguments
  weighting <- match.arg(weighting)
  adaptive <- match.arg(adaptive)
  method <- match.arg(method)

  # Input validation
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("data_list must be a non-empty list of matrices")
  }

  if (is.null(names(data_list))) {
    names(data_list) <- paste0("Block", seq_along(data_list))
  }

  # Check all matrices have same number of rows
  N_vals <- sapply(data_list, nrow)
  if (length(unique(N_vals)) > 1) {
    stop("All data matrices must have the same number of rows (samples)")
  }
  N <- N_vals[1]
  n_blocks <- length(data_list)

  if (K < 2 || K > N) {
    stop("K must be between 2 and number of samples")
  }

  # Set default data types
  if (is.null(data_types)) {
    data_types <- setNames(rep("continuous", length(data_list)), names(data_list))
  } else {
    if (length(data_types) != length(data_list)) {
      stop("data_types must have same length as data_list")
    }
    if (is.null(names(data_types))) {
      names(data_types) <- names(data_list)
    }
  }

  # Validate data types
  valid_types <- c("continuous", "compositional", "simplex", "clr")
  if (!all(data_types %in% valid_types)) {
    stop("data_types must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Transform data as needed
  X_list <- list()
  for (i in seq_along(data_list)) {
    name <- names(data_list)[i]
    type <- data_types[name]
    X <- data_list[[i]]

    if (type == "compositional") {
      # Apply CLR transformation
      if (any(X <= 0)) {
        X <- X + 1e-6
      }
      X_comp <- X / rowSums(X)
      X_clr <- as.matrix(compositions::clr(compositions::acomp(X_comp)))
      X_list[[name]] <- X_clr

    } else if (type == "clr") {
      X_list[[name]] <- X

    } else if (type == "simplex") {
      # Ensure rows sum to 1
      X_list[[name]] <- X / rowSums(X)

    } else {
      # Continuous data, use as-is
      X_list[[name]] <- X
    }
  }

  # Select solver functions
  solve_B <- if (method == "nnls") nnls_B_given_A else ols_B_given_A
  solve_A_rows <- if (method == "nnls") nnls_A_rows else ols_A_rows

  # Compute global block weights
  F_i <- sapply(X_list, ncol)
  global_weights <- switch(weighting,
                           "none" = rep(1, n_blocks),
                           "inv_f" = 1 / F_i,
                           "sqrt_inv_f" = 1 / sqrt(F_i)
  )
  global_weights[!is.finite(global_weights)] <- 1
  global_weights <- global_weights / sum(global_weights) * n_blocks
  names(global_weights) <- names(X_list)

  if (verbose) {
    cat("=== Unmix Settings ===\n")
    cat(sprintf("Samples: %d, Sources: %d, Blocks: %d\n", N, K, n_blocks))
    cat(sprintf("Weighting: %s, Adaptive: %s, Method: %s\n", weighting, adaptive, method))
    cat("Block weights:\n")
    for (nm in names(global_weights)) {
      cat(sprintf("  %s: %.4f (%d features)\n", nm, global_weights[nm], F_i[nm]))
    }
  }

  # Initialize
  set.seed(seed)
  A <- row_simplex(matrix(runif(N * K), N, K))

  # Initialize B for each block
  B_list <- lapply(X_list, function(X) solve_B(A, X))
  names(B_list) <- names(X_list)

  # Initialize sample weights
  W <- matrix(1, N, n_blocks)
  colnames(W) <- names(X_list)
  if (!is.null(rownames(data_list[[1]]))) {
    rownames(W) <- rownames(data_list[[1]])
  }

  # Track loss
  loss_history <- numeric(max_iter)
  prev_loss <- compute_weighted_loss(A, B_list, X_list, W, global_weights)
  loss_history[1] <- prev_loss

  converged <- FALSE

  # Main optimization loop
  for (t in 1:max_iter) {

    # Update B for each block
    for (i in seq_along(B_list)) {
      B_list[[i]] <- solve_B(A, X_list[[i]])
    }

    # Update sample weights (every 10 iterations if adaptive)
    if (adaptive == "entropy" && t %% 10 == 0) {
      W <- compute_entropy_weights(A, B_list, X_list, K)
    }

    # Update A jointly using all blocks with weights
    A_new <- solve_A_joint_weighted(B_list, X_list, W, global_weights, method)

    # Gradient step with momentum
    A <- (1 - eta_A) * A + eta_A * A_new
    A <- row_simplex(A)

    # Compute loss
    cur_loss <- compute_weighted_loss(A, B_list, X_list, W, global_weights)
    loss_history[t] <- cur_loss

    # Check convergence
    rel_change <- abs(prev_loss - cur_loss) / max(prev_loss, 1e-12)

    if (verbose && (t %% 100 == 0 || t <= 5)) {
      cat(sprintf("Iter %d: loss=%.6g, rel_change=%.3e\n", t, cur_loss, rel_change))
    }

    if (rel_change < tol) {
      converged <- TRUE
      if (verbose) {
        cat(sprintf("Converged at iteration %d\n", t))
      }
      loss_history <- loss_history[1:t]
      break
    }

    prev_loss <- cur_loss
  }

  # Final sample weights
  if (adaptive == "entropy") {
    W <- compute_entropy_weights(A, B_list, X_list, K)
  }

  # Reconstruct data
  recon_list <- lapply(B_list, function(B) A %*% B)

  # Add row/col names
  rownames(A) <- rownames(data_list[[1]])
  colnames(A) <- paste0("Source", 1:K)

  for (i in seq_along(B_list)) {
    rownames(B_list[[i]]) <- paste0("Source", 1:K)
    colnames(B_list[[i]]) <- colnames(X_list[[i]])
  }

  # Return results
  result <- list(
    A = A,
    B_list = B_list,
    recon_list = recon_list,
    X_list = X_list,
    data_types = data_types,
    data_names = names(data_list),
    block_weights = global_weights,
    sample_weights = W,
    loss = prev_loss,
    loss_history = loss_history,
    iters = ifelse(converged, t, max_iter),
    converged = converged,
    K = K,
    N = N,
    settings = list(
      weighting = weighting,
      adaptive = adaptive,
      method = method
    ),
    call = match.call()
  )

  class(result) <- "unmix_result"
  return(result)
}

#' Print Method for unmix_result Objects
#'
#' @param x An unmix_result object
#' @param ... Additional arguments (unused)
#' @export
print.unmix_result <- function(x, ...) {
  cat("Unmixing Results\n")
  cat("================\n\n")
  cat("Data blocks:", length(x$B_list), "\n")
  for (i in seq_along(x$data_names)) {
    name <- x$data_names[i]
    type <- x$data_types[name]
    ncol_i <- ncol(x$X_list[[name]])
    weight_i <- x$block_weights[name]
    cat(sprintf("  - %s: %s (%d features, weight=%.3f)\n", name, type, ncol_i, weight_i))
  }
  cat("\n")
  cat("Settings:\n")
  cat(sprintf("  Weighting: %s\n", x$settings$weighting))
  cat(sprintf("  Adaptive:  %s\n", x$settings$adaptive))
  cat(sprintf("  Method:    %s\n", x$settings$method))
  cat("\n")
  cat("Number of sources (K):", x$K, "\n")
  cat("Number of samples (N):", x$N, "\n")
  cat("Iterations:", x$iters, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Final loss:", format(x$loss, digits = 6), "\n")
  invisible(x)
}

#' Summary Method for unmix_result Objects
#'
#' @param object An unmix_result object
#' @param ... Additional arguments (unused)
#' @export
summary.unmix_result <- function(object, ...) {
  cat("Unmixing Summary\n")
  cat("================\n\n")
  print(object)
  cat("\nMixing Matrix A (first 6 samples):\n")
  print(head(object$A))

  cat("\nReconstruction Error by Block:\n")
  for (nm in names(object$B_list)) {
    recon <- object$A %*% object$B_list[[nm]]
    rmse <- sqrt(mean((object$X_list[[nm]] - recon)^2))
    cat(sprintf("  %s: RMSE = %.6f\n", nm, rmse))
  }

  invisible(object)
}
