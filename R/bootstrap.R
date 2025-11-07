#' Bootstrap Uncertainty Estimation for Provenance Unmixing
#' 
#' Performs non-parametric bootstrap to estimate uncertainty in mixing 
#' proportions (A matrix) and source signatures (B matrices).
#' 
#' @param data_list Named list of data matrices (same as provenance_unmix)
#' @param data_types Character vector specifying data types
#' @param K Number of sources
#' @param n_boot Number of bootstrap iterations (default: 100)
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param parallel Logical, use parallel processing? (default: FALSE)
#' @param n_cores Number of cores if parallel=TRUE (default: detectCores()-1)
#' @param verbose Logical, print progress? (default: TRUE)
#' @param ... Additional arguments passed to provenance_unmix
#' 
#' @return A list of class "provunmix_boot" containing:
#'   \item{fit_original}{Original fit on full dataset}
#'   \item{A_boot}{Array (n_boot x N x K) of bootstrapped A matrices}
#'   \item{B_boot_list}{List of arrays (n_boot x K x F) for each data block}
#'   \item{A_mean}{Mean of bootstrapped A estimates}
#'   \item{A_lower}{Lower confidence bound for A}
#'   \item{A_upper}{Upper confidence bound for A}
#'   \item{A_se}{Standard error for A}
#'   \item{conf_level}{Confidence level used}
#'   \item{n_boot}{Number of bootstrap iterations}
#' 
#' @examples
#' \dontrun{
#' # Bootstrap uncertainty estimation
#' boot_result <- bootstrap_provunmix(
#'   data_list = list(DZ = DZ, BP = BP_counts),
#'   data_types = c(DZ = "continuous", BP = "compositional"),
#'   K = 3,
#'   n_boot = 100
#' )
#' 
#' # Plot confidence intervals
#' plot(boot_result, type = "uncertainty")
#' }
#' 
#' @export
bootstrap_provunmix <- function(data_list,
                                 data_types = NULL,
                                 K,
                                 n_boot = 100,
                                 conf_level = 0.95,
                                 parallel = FALSE,
                                 n_cores = NULL,
                                 verbose = TRUE,
                                 ...) {
  
  if (verbose) {
    message(sprintf("Starting bootstrap with %d iterations...", n_boot))
  }
  
  # Original fit
  fit_orig <- provenance_unmix(
    data_list = data_list,
    data_types = data_types,
    K = K,
    verbose = FALSE,
    ...
  )
  
  N <- fit_orig$N
  n_blocks <- length(data_list)
  
  # Storage for bootstrap results
  A_boot <- array(NA, dim = c(n_boot, N, K))
  B_boot_list <- vector("list", n_blocks)
  names(B_boot_list) <- names(data_list)
  
  for (i in seq_along(B_boot_list)) {
    F_i <- ncol(fit_orig$B_list[[i]])
    B_boot_list[[i]] <- array(NA, dim = c(n_boot, K, F_i))
  }
  
  # Bootstrap loop
  bootstrap_iter <- function(b) {
    # Resample with replacement
    idx <- sample(1:N, N, replace = TRUE)
    
    # Create resampled data
    data_boot <- lapply(data_list, function(X) X[idx, , drop = FALSE])
    
    # Fit on bootstrap sample
    fit_boot <- tryCatch({
      provenance_unmix(
        data_list = data_boot,
        data_types = data_types,
        K = K,
        verbose = FALSE,
        ...
      )
    }, error = function(e) {
      if (verbose) message(sprintf("Bootstrap iteration %d failed: %s", b, e$message))
      return(NULL)
    })
    
    if (is.null(fit_boot)) {
      return(NULL)
    }
    
    # Extract results
    result <- list(
      A = fit_boot$A,
      B_list = fit_boot$B_list
    )
    
    return(result)
  }
  
  # Run bootstrap
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    if (verbose) {
      message(sprintf("Using %d cores for parallel processing", n_cores))
    }
    
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    # Export necessary objects
    parallel::clusterExport(cl, c("data_list", "data_types", "K", "provenance_unmix"),
                           envir = environment())
    
    boot_results <- parallel::parLapply(cl, 1:n_boot, bootstrap_iter)
    
  } else {
    boot_results <- vector("list", n_boot)
    for (b in 1:n_boot) {
      if (verbose && b %% 10 == 0) {
        message(sprintf("Bootstrap iteration %d/%d", b, n_boot))
      }
      boot_results[[b]] <- bootstrap_iter(b)
    }
  }
  
  # Store results
  for (b in seq_along(boot_results)) {
    if (!is.null(boot_results[[b]])) {
      A_boot[b, , ] <- boot_results[[b]]$A
      for (i in seq_along(B_boot_list)) {
        B_boot_list[[i]][b, , ] <- boot_results[[b]]$B_list[[i]]
      }
    }
  }
  
  # Remove failed iterations
  valid_iters <- !apply(A_boot, 1, function(x) all(is.na(x)))
  A_boot <- A_boot[valid_iters, , , drop = FALSE]
  for (i in seq_along(B_boot_list)) {
    B_boot_list[[i]] <- B_boot_list[[i]][valid_iters, , , drop = FALSE]
  }
  
  n_valid <- sum(valid_iters)
  if (verbose) {
    message(sprintf("Completed %d/%d successful bootstrap iterations", n_valid, n_boot))
  }
  
  # Compute statistics
  alpha <- 1 - conf_level
  
  A_mean <- apply(A_boot, c(2, 3), mean, na.rm = TRUE)
  A_se <- apply(A_boot, c(2, 3), sd, na.rm = TRUE)
  A_lower <- apply(A_boot, c(2, 3), quantile, probs = alpha/2, na.rm = TRUE)
  A_upper <- apply(A_boot, c(2, 3), quantile, probs = 1 - alpha/2, na.rm = TRUE)
  
  # Add names
  rownames(A_mean) <- rownames(fit_orig$A)
  colnames(A_mean) <- colnames(fit_orig$A)
  rownames(A_se) <- rownames(A_lower) <- rownames(A_upper) <- rownames(A_mean)
  colnames(A_se) <- colnames(A_lower) <- colnames(A_upper) <- colnames(A_mean)
  
  # Return results
  result <- list(
    fit_original = fit_orig,
    A_boot = A_boot,
    B_boot_list = B_boot_list,
    A_mean = A_mean,
    A_lower = A_lower,
    A_upper = A_upper,
    A_se = A_se,
    conf_level = conf_level,
    n_boot = n_valid,
    call = match.call()
  )
  
  class(result) <- "provunmix_boot"
  return(result)
}

#' Print Method for provunmix_boot Objects
#' 
#' @param x A provunmix_boot object
#' @param ... Additional arguments (unused)
#' @export
print.provunmix_boot <- function(x, ...) {
  cat("Bootstrap Provenance Unmixing Results\n")
  cat("======================================\n\n")
  cat("Number of bootstrap iterations:", x$n_boot, "\n")
  cat("Confidence level:", x$conf_level * 100, "%\n")
  cat("\nOriginal fit:\n")
  print(x$fit_original)
  cat("\nBootstrap statistics (A matrix):\n")
  cat("Mean standard error:", mean(x$A_se), "\n")
  cat("Range of standard errors:", range(x$A_se), "\n")
  invisible(x)
}

#' Summary Method for provunmix_boot Objects
#' 
#' @param object A provunmix_boot object
#' @param ... Additional arguments (unused)
#' @export
summary.provunmix_boot <- function(object, ...) {
  cat("Bootstrap Provenance Unmixing Summary\n")
  cat("======================================\n\n")
  print(object)
  
  cat("\n\nMixing Matrix A with Confidence Intervals (first 3 samples, first 3 sources):\n")
  for (i in 1:min(3, nrow(object$A_mean))) {
    cat(sprintf("\nSample %s:\n", rownames(object$A_mean)[i]))
    for (k in 1:min(3, ncol(object$A_mean))) {
      cat(sprintf("  Source%d: %.3f (95%% CI: %.3f - %.3f, SE: %.3f)\n",
                  k,
                  object$fit_original$A[i, k],
                  object$A_lower[i, k],
                  object$A_upper[i, k],
                  object$A_se[i, k]))
    }
  }
  
  invisible(object)
}
