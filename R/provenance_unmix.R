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
#'   "clr" (already in CLR space).
#' @param K Integer, number of sources (end-members) to unmix.
#' @param max_iter Maximum number of iterations (default: 2000).
#' @param tol Convergence tolerance for relative change in loss (default: 1e-8).
#' @param eta_A Step size for updating A (default: 1.0, full update).
#' @param seed Random seed for reproducibility (default: 123).
#' @param verbose Logical, print progress? (default: FALSE).
#'
#' @return A list of class "unmix_result" containing:
#'   \item{A}{Mixing matrix (N x K), rows sum to 1}
#'   \item{B_list}{List of source signature matrices (K x F_i for each block)}
#'   \item{recon_list}{List of reconstructed data matrices}
#'   \item{data_types}{Character vector of data types used}
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
#' # Unmix
#' result <- unmix(
#'   data_list = list(DZ = DZ, BP = BP_counts),
#'   data_types = c(DZ = "continuous", BP = "compositional"),
#'   K = 3
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
                              max_iter = 2000, 
                              tol = 1e-8,
                              eta_A = 1.0, 
                              seed = 123, 
                              verbose = FALSE) {
  
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
  valid_types <- c("continuous", "compositional", "clr")
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
      # Add small pseudocount if needed
      if (any(X <= 0)) {
        X <- X + 1e-6
      }
      # Close to simplex
      X_comp <- X / rowSums(X)
      # CLR transform
      X_clr <- as.matrix(compositions::clr(compositions::acomp(X_comp)))
      X_list[[name]] <- X_clr
      
    } else if (type == "clr") {
      # Already in CLR space
      X_list[[name]] <- X
      
    } else {
      # Continuous data, use as-is
      X_list[[name]] <- X
    }
  }
  
  # Initialize
  set.seed(seed)
  A <- row_simplex(matrix(runif(N * K), N, K))
  
  # Initialize B for each block
  B_list <- lapply(X_list, function(X) nnls_B_given_A(A, X))
  
  # Track loss
  loss_history <- numeric(max_iter)
  prev_loss <- compute_loss(A, B_list, X_list)
  loss_history[1] <- prev_loss
  
  converged <- FALSE
  
  # Main optimization loop
  for (t in 1:max_iter) {
    
    # Update B for each block
    for (i in seq_along(B_list)) {
      B_list[[i]] <- nnls_B_given_A(A, X_list[[i]])
    }
    
    # Update A jointly using all blocks
    A_new <- solve_A_joint_fair(B_list, X_list)
    
    # Gradient step
    A <- (1 - eta_A) * A + eta_A * A_new
    A <- row_simplex(A)
    
    # Compute loss
    cur_loss <- compute_loss(A, B_list, X_list)
    loss_history[t] <- cur_loss
    
    # Check convergence
    rel_change <- abs(prev_loss - cur_loss) / max(prev_loss, 1e-12)
    
    if (verbose && (t %% 10 == 0)) {
      message(sprintf("Iter %d: loss=%.6g, rel_change=%.3e", t, cur_loss, rel_change))
    }
    
    if (rel_change < tol) {
      converged <- TRUE
      if (verbose) {
        message(sprintf("Converged at iteration %d", t))
      }
      loss_history <- loss_history[1:t]
      break
    }
    
    prev_loss <- cur_loss
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
    loss = prev_loss,
    loss_history = loss_history,
    iters = ifelse(converged, t, max_iter),
    converged = converged,
    K = K,
    N = N,
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
    cat(sprintf("  - %s: %s (%d features)\n", name, type, ncol_i))
  }
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
  invisible(object)
}
