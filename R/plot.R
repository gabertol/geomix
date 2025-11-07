#' Plot Methods for Provenance Unmixing Results
#' 
#' Visualize provenance unmixing results including mixing proportions,
#' source signatures, and uncertainty estimates.
#' 
#' @param x A provunmix or provunmix_boot object
#' @param type Type of plot: "mixing" (default), "sources", "loss", "uncertainty"
#' @param ... Additional plotting parameters
#' 
#' @export
plot.provunmix <- function(x, type = c("mixing", "sources", "loss"), ...) {
  type <- match.arg(type)
  
  if (type == "mixing") {
    plot_mixing_matrix(x, ...)
  } else if (type == "sources") {
    plot_sources(x, ...)
  } else if (type == "loss") {
    plot_loss_history(x, ...)
  }
}

#' @export
plot.provunmix_boot <- function(x, type = c("mixing", "uncertainty", "sources"), ...) {
  type <- match.arg(type)
  
  if (type == "mixing") {
    plot_mixing_matrix(x$fit_original, ...)
  } else if (type == "uncertainty") {
    plot_uncertainty(x, ...)
  } else if (type == "sources") {
    plot_sources(x$fit_original, ...)
  }
}

#' Plot Mixing Matrix
#' 
#' @param x A provunmix object
#' @param samples Optional vector of sample indices to plot (default: all)
#' @param ... Additional parameters
#' @keywords internal
plot_mixing_matrix <- function(x, samples = NULL, ...) {
  A <- x$A
  
  if (is.null(samples)) {
    samples <- 1:nrow(A)
  }
  
  A_plot <- A[samples, , drop = FALSE]
  
  # Stacked barplot
  par(mar = c(5, 4, 4, 8), xpd = TRUE)
  
  colors <- grDevices::rainbow(x$K, alpha = 0.7)
  
  bp <- barplot(t(A_plot), 
                col = colors,
                border = NA,
                ylab = "Mixing Proportion",
                xlab = "Sample",
                main = "Source Mixing Proportions",
                names.arg = rownames(A_plot),
                las = 2,
                ...)
  
  legend("topright", 
         inset = c(-0.2, 0),
         legend = colnames(A),
         fill = colors,
         bty = "n")
  
  invisible(bp)
}

#' Plot Source Signatures
#' 
#' @param x A provunmix object
#' @param block Which data block to plot (default: 1)
#' @param ... Additional parameters
#' @keywords internal
plot_sources <- function(x, block = 1, ...) {
  
  if (block > length(x$B_list)) {
    stop("block index out of range")
  }
  
  B <- x$B_list[[block]]
  block_name <- x$data_names[block]
  
  # Heatmap of source signatures
  par(mar = c(5, 4, 4, 2))
  
  image(1:ncol(B), 1:nrow(B), t(B),
        col = grDevices::hcl.colors(50, "YlOrRd", rev = TRUE),
        xlab = "Feature",
        ylab = "Source",
        main = paste("Source Signatures:", block_name),
        axes = FALSE,
        ...)
  
  axis(1, at = 1:ncol(B), labels = colnames(B), las = 2, cex.axis = 0.8)
  axis(2, at = 1:nrow(B), labels = rownames(B), las = 1)
  
  invisible(NULL)
}

#' Plot Loss History
#' 
#' @param x A provunmix object
#' @param log Logical, use log scale for y-axis? (default: TRUE)
#' @param ... Additional parameters
#' @keywords internal
plot_loss_history <- function(x, log = TRUE, ...) {
  par(mar = c(5, 4, 4, 2))
  
  if (log) {
    plot(x$loss_history, 
         type = "l",
         log = "y",
         xlab = "Iteration",
         ylab = "Loss (log scale)",
         main = "Convergence History",
         col = "blue",
         lwd = 2,
         ...)
  } else {
    plot(x$loss_history, 
         type = "l",
         xlab = "Iteration",
         ylab = "Loss",
         main = "Convergence History",
         col = "blue",
         lwd = 2,
         ...)
  }
  
  grid()
  invisible(NULL)
}

#' Plot Uncertainty Estimates
#' 
#' @param x A provunmix_boot object
#' @param samples Optional vector of sample indices to plot (default: first 10)
#' @param ... Additional parameters
#' @keywords internal
plot_uncertainty <- function(x, samples = NULL, ...) {
  
  if (is.null(samples)) {
    samples <- 1:min(10, nrow(x$A_mean))
  }
  
  A_orig <- x$fit_original$A[samples, , drop = FALSE]
  A_lower <- x$A_lower[samples, , drop = FALSE]
  A_upper <- x$A_upper[samples, , drop = FALSE]
  
  K <- ncol(A_orig)
  n_samples <- length(samples)
  
  par(mfrow = c(1, K), mar = c(5, 4, 4, 2))
  
  for (k in 1:K) {
    # Error bars
    x_pos <- 1:n_samples
    plot(x_pos, A_orig[, k],
         ylim = c(0, 1),
         pch = 19,
         xlab = "Sample",
         ylab = "Proportion",
         main = paste("Source", k),
         xaxt = "n",
         ...)
    
    axis(1, at = x_pos, labels = rownames(A_orig), las = 2, cex.axis = 0.8)
    
    # Add error bars
    arrows(x_pos, A_lower[, k], x_pos, A_upper[, k],
           angle = 90, code = 3, length = 0.05, col = "gray50")
    
    grid()
  }
  
  par(mfrow = c(1, 1))
  invisible(NULL)
}
