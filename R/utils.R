#' @title Utility Functions for Source Unmixing
#' @name utils
#' @keywords internal


#' Row-wise Simplex Normalization
#'
#' Normalizes rows of a matrix to sum to 1 (simplex constraint).
#'
#' @param M Numeric matrix
#' @param eps Small value to avoid division by zero (default: 1e-12)
#' @return Matrix with rows normalized to sum to 1
#' @keywords internal
row_simplex <- function(M, eps = 1e-12) {
  M[M < 0] <- 0
  rs <- rowSums(M)
  rs[!is.finite(rs) | rs < eps] <- 1
  sweep(M, 1, rs, "/")
}

#' Inverse CLR Transformation (Row-wise)
#'
#' Applies inverse centered log-ratio transformation to each row of a matrix.
#'
#' @param M Matrix in CLR space (rows are log-ratios)
#' @param colnames_keep Optional column names to preserve
#' @return Matrix in simplex space (rows sum to 1)
#' @keywords internal
#' @importFrom compositions clrInv
clr_inv_rows <- function(M, colnames_keep = NULL) {
  out <- t(apply(M, 1, function(r) as.numeric(compositions::clrInv(r))))
  if (!is.null(colnames_keep)) colnames(out) <- colnames_keep
  out
}

# =============================================================================
# NNLS SOLVERS
# =============================================================================

#' NNLS for B Given A
#'
#' Solves for B in X ~ A %*% B using non-negative least squares,
#' column by column.
#'
#' @param A Matrix (N x K) - mixing proportions
#' @param X Matrix (N x F) - data
#' @return Matrix B (K x F) - source signatures
#' @keywords internal
#' @importFrom nnls nnls
nnls_B_given_A <- function(A, X) {
  K <- ncol(A)
  n_features <- ncol(X)
  B <- matrix(0, K, n_features)
  for (j in seq_len(n_features)) {
    B[, j] <- pmax(0, as.numeric(nnls::nnls(A, X[, j])$x))
  }
  B
}

#' NNLS for A Given B (Row-wise)
#'
#' Solves for A in X ~ A %*% B using non-negative least squares,
#' row by row.
#'
#' @param B Matrix (K x F) - source signatures
#' @param X Matrix (N x F) - data
#' @return Matrix A (N x K) - mixing proportions
#' @keywords internal
#' @importFrom nnls nnls
nnls_A_rows <- function(B, X) {
  K <- nrow(B)
  N <- nrow(X)
  Bt <- t(B)
  A <- matrix(0, N, K)
  for (i in seq_len(N)) {
    A[i, ] <- pmax(0, as.numeric(nnls::nnls(Bt, X[i, ])$x))
  }
  A
}

# =============================================================================
# OLS SOLVERS
# =============================================================================

#' OLS for B Given A
#'
#' Solves for B in X ~ A %*% B using ordinary least squares.
#'
#' @param A Matrix (N x K) - mixing proportions
#' @param X Matrix (N x F) - data
#' @return Matrix B (K x F) - source signatures
#' @keywords internal
ols_B_given_A <- function(A, X) {
  AtA <- t(A) %*% A
  AtA_inv <- tryCatch(
    solve(AtA),
    error = function(e) MASS::ginv(AtA)
  )
  AtA_inv %*% t(A) %*% X
}

#' OLS for A Given B (Row-wise)
#'
#' Solves for A in X ~ A %*% B using ordinary least squares,
#' row by row.
#'
#' @param B Matrix (K x F) - source signatures
#' @param X Matrix (N x F) - data
#' @return Matrix A (N x K) - mixing proportions
#' @keywords internal
ols_A_rows <- function(B, X) {
  K <- nrow(B)
  N <- nrow(X)
  Bt <- t(B)
  BtB <- t(Bt) %*% Bt
  BtB_inv <- tryCatch(
    solve(BtB),
    error = function(e) MASS::ginv(BtB)
  )
  A <- matrix(0, N, K)
  for (i in seq_len(N)) {
    A[i, ] <- as.vector(BtB_inv %*% t(Bt) %*% X[i, ])
  }
  A
}

# =============================================================================
# BLOCK NORMALIZATION
# =============================================================================

#' Normalize Data Block by Column
#'
#' Normalizes both B (source signatures) and X (data) by column-wise scaling.
#' This ensures fair weighting between different data blocks.
#'
#' @param B Matrix (K x F) - source signatures
#' @param X Matrix (N x F) - data
#' @return List with normalized Bn and Xn
#' @keywords internal
normalize_block_by_col <- function(B, X) {
  s <- sqrt(colSums(B^2))
  s[!is.finite(s) | s == 0] <- 1
  list(
    Bn = sweep(B, 2, s, "/"),
    Xn = sweep(X, 2, s, "/")
  )
}

# =============================================================================
# WEIGHTED SOLVERS
# =============================================================================

#' Solve for A Using Multiple Data Blocks with Weights
#'
#' Jointly solves for mixing matrix A using multiple data blocks
#' with global and sample-wise weights.
#'
#' @param B_list List of source signature matrices (one per data block)
#' @param X_list List of data matrices (one per data block)
#' @param W Matrix of sample weights (N x n_blocks)
#' @param global_weights Named vector of global block weights
#' @param method Solver: "nnls" or "ols"
#' @return Matrix A (N x K) with rows normalized to simplex
#' @keywords internal
solve_A_joint_weighted <- function(B_list, X_list, W, global_weights, method = "nnls") {
  N <- nrow(X_list[[1]])
  K <- nrow(B_list[[1]])

  # Normalize each block
  normalized_list <- mapply(normalize_block_by_col, B_list, X_list, SIMPLIFY = FALSE)

  A_new <- matrix(0, N, K)

  for (i in seq_len(N)) {
    B_parts <- X_parts <- list()
    for (b in seq_along(B_list)) {
      w <- sqrt(W[i, b] * global_weights[b])
      B_parts[[b]] <- normalized_list[[b]]$Bn * w
      X_parts[[b]] <- normalized_list[[b]]$Xn[i, ] * w
    }
    B_concat <- do.call(cbind, B_parts)
    x_concat <- unlist(X_parts)

    if (method == "nnls") {
      fit <- nnls::nnls(t(B_concat), x_concat)
      A_new[i, ] <- pmax(0, as.numeric(fit$x))
    } else {
      Bt <- t(B_concat)
      BtB_inv <- tryCatch(
        solve(t(Bt) %*% Bt),
        error = function(e) MASS::ginv(t(Bt) %*% Bt)
      )
      A_new[i, ] <- as.vector(BtB_inv %*% t(Bt) %*% x_concat)
    }
  }

  row_simplex(A_new)
}

#' Solve for A Using Multiple Data Blocks (Legacy - no weights)
#'
#' @param B_list List of source signature matrices
#' @param X_list List of data matrices
#' @return Matrix A (N x K) with rows normalized to simplex
#' @keywords internal
solve_A_joint_fair <- function(B_list, X_list) {
  W <- matrix(1, nrow(X_list[[1]]), length(B_list))
  global_weights <- rep(1, length(B_list))
  names(global_weights) <- names(B_list)
  solve_A_joint_weighted(B_list, X_list, W, global_weights, "nnls")
}

# =============================================================================
# ADAPTIVE WEIGHTING
# =============================================================================

#' Compute Entropy-Based Sample Weights
#'
#' Computes per-sample weights based on the entropy of single-block estimates.
#' Low entropy (confident estimate) gets higher weight.
#'
#' @param A Current mixing matrix (N x K)
#' @param B_list List of source signatures
#' @param X_list List of data matrices
#' @param K Number of sources
#' @param floor_weight Minimum weight (default: 0.01)
#' @return Matrix of weights (N x n_blocks)
#' @keywords internal
compute_entropy_weights <- function(A, B_list, X_list, K, floor_weight = 0.01) {
  N <- nrow(A)
  n_blocks <- length(B_list)

  W <- matrix(1, N, n_blocks)
  colnames(W) <- names(X_list)

  for (b in seq_along(B_list)) {
    # Estimate A using only this block
    A_single <- nnls_A_rows(B_list[[b]], X_list[[b]])
    A_single <- row_simplex(A_single)

    # Compute entropy per sample (low entropy = confident)
    entropy <- -rowSums(A_single * log(A_single + 1e-10))
    max_entropy <- log(K)

    # Weight inversely proportional to entropy
    W[, b] <- 1 / (entropy / max_entropy + 0.1)
  }

  W[W < floor_weight] <- floor_weight
  W <- W / rowSums(W) * n_blocks
  W
}

# =============================================================================
# LOSS FUNCTIONS
# =============================================================================

#' Loss Function for Coupled Unmixing
#'
#' Computes total reconstruction error across all data blocks.
#'
#' @param A Matrix (N x K) - mixing proportions
#' @param B_list List of source signature matrices
#' @param X_list List of data matrices
#' @return Scalar loss value (sum of squared errors)
#' @keywords internal
compute_loss <- function(A, B_list, X_list) {
  total_loss <- 0
  for (i in seq_along(B_list)) {
    resid <- X_list[[i]] - A %*% B_list[[i]]
    total_loss <- total_loss + sum(resid^2)
  }
  total_loss
}

#' Weighted Loss Function
#'
#' Computes weighted reconstruction error across all data blocks.
#'
#' @param A Matrix (N x K) - mixing proportions
#' @param B_list List of source signature matrices
#' @param X_list List of data matrices
#' @param W Matrix of sample weights (N x n_blocks)
#' @param global_weights Named vector of global block weights
#' @return Scalar loss value
#' @keywords internal
compute_weighted_loss <- function(A, B_list, X_list, W, global_weights) {
  total_loss <- 0
  for (b in seq_along(B_list)) {
    recon <- A %*% B_list[[b]]
    resid_sq <- (X_list[[b]] - recon)^2
    weighted_resid <- rowSums(resid_sq) * W[, b] * global_weights[b]
    total_loss <- total_loss + sum(weighted_resid)
  }
  total_loss
}

# =============================================================================
# DISTANCE FUNCTIONS
# =============================================================================

#' Hellinger Distance
#'
#' Computes Hellinger distance between two distributions.
#'
#' @param p First distribution (will be normalized)
#' @param q Second distribution (will be normalized)
#' @return Hellinger distance (0 to 1)
#' @export
hellinger_distance <- function(p, q) {
  p <- p / max(sum(p), 1e-10)
  q <- q / max(sum(q), 1e-10)
  sqrt(0.5 * sum((sqrt(p) - sqrt(q))^2))
}

# =============================================================================
# LOCAL MIXING PROPORTIONS
# =============================================================================

#' Solve for A Given Fixed B (Local Mixing Proportions)
#'
#' Computes sample-wise mixing proportions (A) for a single data block
#' using fixed source signatures (B) from a previous unmixing fit.
#' This allows analyzing how mixing proportions vary across different
#' data types while maintaining consistent source definitions.
#'
#' @param X Matrix (N x F) of observed data for one block.
#'   Rows are samples, columns are features.
#' @param B Matrix (K x F) of source signatures for this block.
#'   Typically obtained from `fit$B_list$BlockName`.
#' @param method Solver method: "nnls" (default) or "ols".
#' @param verbose Logical, print diagnostic info? (default: FALSE)
#'
#' @return Matrix A_local (N x K) with mixing proportions.
#'   Rows sum to 1 (simplex constraint).
#'
#' @details
#' Given the linear mixing model X ≈ A %*% B, this function solves
#' for A row-by-row using non-negative least squares (NNLS):
#'
#'   x_i ≈ t(B) %*% a_i
#'
#' Each row a_i is then normalized to sum to 1.
#'
#' **Use case**: After running `unmix()` with multiple blocks to get
#' a global A and block-specific B matrices, use this function to
#' compute "local" mixing proportions for each block individually.
#' Comparing A_global vs A_local reveals which samples have consistent
#' mixing across proxies vs which show proxy-dependent variation.
#'
#' @examples
#' \dontrun{
#' # 1. Run joint unmixing
#' fit <- unmix(
#'   data_list = list(DZ = DZ_mat, HM = HM_mat, PT = PT_mat),
#'   data_types = c(DZ = "continuous", HM = "simplex", PT = "simplex"),
#'   K = 3
#' )
#'
#' # 2. Compute local A for each block
#' A_local_DZ <- solve_A_given_B(fit$X_list$DZ, fit$B_list$DZ)
#' A_local_HM <- solve_A_given_B(fit$X_list$HM, fit$B_list$HM)
#'
#' # 3. Compare with global
#' mae_DZ <- mean(abs(fit$A - A_local_DZ))
#' }
#'
#' @seealso \code{\link{compute_A_local}}, \code{\link{unmix_cv}}
#' @export
#' @importFrom nnls nnls
#' @importFrom MASS ginv
solve_A_given_B <- function(X, B, method = c("nnls", "ols"), verbose = FALSE) {

  method <- match.arg(method)

  X <- as.matrix(X)
  B <- as.matrix(B)

  N   <- nrow(X)
  F_x <- ncol(X)
  K   <- nrow(B)
  F_b <- ncol(B)

  if (verbose) {
    cat(sprintf("  X: %d samples x %d features\n", N, F_x))
    cat(sprintf("  B: %d sources x %d features\n", K, F_b))
  }

  if (F_x != F_b) {
    stop(sprintf(
      "Feature dimension mismatch: X has %d columns, B has %d columns",
      F_x, F_b
    ))
  }

  # x_i ≈ A[i,] %*% B  =>  x_i ≈ t(B) %*% a_i
  Bt <- t(B)  # F x K

  A_local <- matrix(NA_real_, nrow = N, ncol = K)

  for (i in seq_len(N)) {
    x_i <- X[i, ]

    if (method == "nnls") {
      fit_i <- nnls::nnls(Bt, x_i)
      a_i   <- coef(fit_i)
    } else {
      # OLS solution
      BtB <- crossprod(Bt)
      BtB_inv <- tryCatch(
        solve(BtB),
        error = function(e) MASS::ginv(BtB)
      )
      a_i <- as.vector(BtB_inv %*% crossprod(Bt, x_i))
      a_i[a_i < 0] <- 0
    }

    # Normalize to simplex
    s <- sum(a_i)
    if (s > 1e-10) {
      a_i <- a_i / s
    } else {
      a_i <- rep(1 / K, K)
    }

    A_local[i, ] <- a_i
  }

  rownames(A_local) <- rownames(X)
  colnames(A_local) <- rownames(B)

  A_local
}


#' Compute A_local for All Blocks
#'
#' Computes local mixing proportions for all data blocks in an unmix result.
#'
#' @param fit An unmix_result object from \code{\link{unmix}}.
#' @param blocks Character vector of block names to process.
#'   Default: all blocks in fit.
#' @param method Solver method: "nnls" (default) or "ols".
#' @param verbose Logical, print progress? (default: TRUE)
#'
#' @return A list of class "A_local_result" containing:
#'   \item{A_global}{The global mixing matrix from the original fit}
#'   \item{A_local}{Named list of A_local matrices, one per block}
#'   \item{comparison}{Data frame with MAE by block}
#'   \item{samples}{Sample names}
#'
#' @examples
#' \dontrun{
#' fit <- unmix(data_list, data_types, K = 3)
#' result <- compute_A_local(fit)
#'
#' # Access local A for specific block
#' result$A_local$DZ
#'
#' # See comparison
#' result$comparison
#' }
#'
#' @export
compute_A_local <- function(fit, blocks = NULL, method = "nnls", verbose = TRUE) {

  if (!inherits(fit, "unmix_result")) {
    stop("fit must be an unmix_result object")
  }

  if (is.null(blocks)) {
    blocks <- names(fit$B_list)
  }

  A_global <- fit$A
  A_local_list <- list()

  for (block in blocks) {
    if (verbose) {
      cat(sprintf("Computing A_local for %s:\n", block))
    }

    X_block <- fit$X_list[[block]]
    B_block <- fit$B_list[[block]]

    if (is.null(X_block) || is.null(B_block)) {
      warning(sprintf("Block '%s' not found in fit, skipping.", block))
      next
    }

    A_local_list[[block]] <- solve_A_given_B(X_block, B_block,
                                             method = method,
                                             verbose = verbose)
  }

  # Compute MAE comparison
  comparison <- data.frame(
    block = names(A_local_list),
    MAE = sapply(A_local_list, function(A_loc) {
      mean(abs(A_global - A_loc))
    }),
    RMSE = sapply(A_local_list, function(A_loc) {
      sqrt(mean((A_global - A_loc)^2))
    }),
    max_diff = sapply(A_local_list, function(A_loc) {
      max(abs(A_global - A_loc))
    }),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n=== A_local vs A_global comparison ===\n")
    print(comparison)
  }

  result <- list(
    A_global = A_global,
    A_local = A_local_list,
    comparison = comparison,
    samples = rownames(A_global),
    K = fit$K
  )

  class(result) <- "A_local_result"
  result
}


#' Compare A_global vs A_local
#'
#' Computes per-block and per-sample MAE between global and local
#' mixing proportions.
#'
#' @param fit An unmix_result object.
#' @param A_local_list Output from \code{\link{compute_A_local}}.
#'   If NULL, will be computed internally.
#'
#' @return A list with:
#'   \item{by_block}{Data frame with MAE per block}
#'   \item{by_sample}{Data frame with MAE per sample per block}
#'
#' @export
compare_A_global_local <- function(fit, A_local_list = NULL) {

  if (is.null(A_local_list)) {
    A_local_list <- compute_A_local(fit, verbose = FALSE)$A_local
  }

  A_global <- fit$A

  # MAE by block
  mae_by_block <- data.frame(
    block = names(A_local_list),
    MAE = sapply(A_local_list, function(A_loc) {
      mean(abs(A_global - A_loc))
    }),
    stringsAsFactors = FALSE
  )

  # MAE by sample
  mae_by_sample <- do.call(rbind, lapply(names(A_local_list), function(block) {
    A_loc <- A_local_list[[block]]
    data.frame(
      sample = rownames(A_global),
      block = block,
      MAE = rowMeans(abs(A_global - A_loc)),
      stringsAsFactors = FALSE
    )
  }))

  list(
    by_block = mae_by_block,
    by_sample = mae_by_sample
  )
}
