#' @title Utility Functions for Source Unmixing
#' @name utils
#' @keywords internal
NULL

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
  F <- ncol(X)
  B <- matrix(0, K, F)
  for (j in seq_len(F)) {
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
  B <- AtA_inv %*% t(A) %*% X
  B
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
  n_blocks <- length(B_list)

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
