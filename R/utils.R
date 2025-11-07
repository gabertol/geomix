#' @title Utility Functions for Provenance Unmixing
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

#' Solve for A Using Multiple Data Blocks
#' 
#' Jointly solves for mixing matrix A using multiple data blocks
#' with fair normalization.
#' 
#' @param B_list List of source signature matrices (one per data block)
#' @param X_list List of data matrices (one per data block)
#' @return Matrix A (N x K) with rows normalized to simplex
#' @keywords internal
solve_A_joint_fair <- function(B_list, X_list) {
  # Normalize each block
  normalized_list <- mapply(normalize_block_by_col, B_list, X_list, SIMPLIFY = FALSE)
  
  # Concatenate normalized blocks
  B_joint <- do.call(cbind, lapply(normalized_list, `[[`, "Bn"))
  X_joint <- do.call(cbind, lapply(normalized_list, `[[`, "Xn"))
  
  # Solve NNLS
  A_hat <- nnls_A_rows(B_joint, X_joint)
  
  # Project to simplex
  row_simplex(A_hat)
}

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
