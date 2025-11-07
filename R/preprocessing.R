#' Data Preprocessing Functions for Provenance Analysis
#' 
#' Functions to prepare raw data for provenance unmixing, including
#' KDE discretization for detrital zircon ages and compositional data
#' preparation for bulk petrology.
#' 
#' @name preprocessing
NULL

#' Build KDE-Discretized Block from Detrital Zircon Ages
#' 
#' Converts raw detrital zircon age measurements into discretized probability
#' density functions (PDFs) via kernel density estimation (KDE). This is the
#' recommended format for the "continuous" data type in provenance_unmix().
#' 
#' @param ages_df Data frame with columns: sample (character/factor) and 
#'   one or more numeric columns containing age measurements (e.g., age_concordia,
#'   ti_temp, eu_eu, etc.). Each row is one grain measurement.
#' @param age_vars Character vector of column names to discretize. If NULL,
#'   uses all numeric columns except 'sample'.
#' @param n_points Integer, number of bins for discretization (default: 129).
#'   Higher = more resolution, but more memory.
#' @param probs Numeric vector of length 2, quantiles to use for grid bounds
#'   (default: c(0.03, 0.97)). This removes extreme outliers.
#' @param bw_method Bandwidth selection method: "SJ" (Sheather-Jones, default)
#'   or "nrd0" (Scott's rule). SJ is generally better but slower.
#' @param bw_multiplier Bandwidth multiplier (default: 1.0). Values > 1 smooth
#'   more, < 1 smooth less.
#' @param samples_order Optional character vector specifying sample order in
#'   output matrix. If NULL, uses alphabetical order.
#' 
#' @return A list with components:
#'   \item{DZ_mat}{Matrix (N x F*n_points) of discretized PDFs, where N is
#'     number of samples and F is number of features (age_vars)}
#'   \item{samples}{Character vector of sample names (row names of DZ_mat)}
#'   \item{age_vars}{Character vector of feature names}
#'   \item{grids}{Named list of evaluation grids (one per feature)}
#'   \item{n_points}{Number of bins per feature}
#' 
#' @details
#' For each sample and each age variable:
#' 1. Compute bandwidth using specified method (SJ or nrd0)
#' 2. Create evaluation grid from quantiles (probs) of all data
#' 3. Evaluate KDE at grid points
#' 4. Normalize to unit area (trapezoidal rule)
#' 
#' The resulting matrix has F*n_points columns, with columns ordered as:
#' [var1_bin1, var1_bin2, ..., var1_binN, var2_bin1, ..., varF_binN]
#' 
#' @examples
#' \dontrun{
#' # Example data
#' dz_raw <- data.frame(
#'   sample = rep(c("S1", "S2", "S3"), each = 100),
#'   age_concordia = c(rnorm(100, 500, 50), rnorm(100, 1000, 100), rnorm(100, 1500, 150)),
#'   ti_temp = rnorm(300, 750, 50),
#'   th_u = rlnorm(300, 0, 0.5)
#' )
#' 
#' # Build KDE block
#' dz_block <- build_dz_kde_block(
#'   ages_df = dz_raw,
#'   age_vars = c("age_concordia", "ti_temp", "th_u"),
#'   n_points = 129
#' )
#' 
#' # Use in provenance_unmix
#' result <- provenance_unmix(
#'   data_list = list(DZ = dz_block$DZ_mat, BP = bp_counts),
#'   data_types = c(DZ = "continuous", BP = "compositional"),
#'   K = 3
#' )
#' }
#' 
#' @export
build_dz_kde_block <- function(ages_df,
                                age_vars = NULL,
                                n_points = 129,
                                probs = c(0.03, 0.97),
                                bw_method = c("SJ", "nrd0"),
                                bw_multiplier = 1.0,
                                samples_order = NULL) {
  
  bw_method <- match.arg(bw_method)
  
  # Validate input
  if (!"sample" %in% names(ages_df)) {
    stop("ages_df must have a 'sample' column")
  }
  
  # Determine age variables
  if (is.null(age_vars)) {
    age_vars <- setdiff(names(ages_df), "sample")
    age_vars <- age_vars[sapply(ages_df[age_vars], is.numeric)]
    if (length(age_vars) == 0) {
      stop("No numeric columns found (besides 'sample')")
    }
    message("Using age variables: ", paste(age_vars, collapse = ", "))
  }
  
  # Check all age_vars exist
  missing_vars <- setdiff(age_vars, names(ages_df))
  if (length(missing_vars) > 0) {
    stop("Variables not found in ages_df: ", paste(missing_vars, collapse = ", "))
  }
  
  # Get unique samples
  samples <- unique(ages_df$sample)
  if (!is.null(samples_order)) {
    if (!all(samples_order %in% samples)) {
      stop("Some samples in samples_order not found in data")
    }
    samples <- samples_order
  } else {
    samples <- sort(samples)
  }
  
  N <- length(samples)
  F <- length(age_vars)
  
  # Helper: get evaluation grid
  get_evalpoints <- function(x, n = n_points, probs = c(0.03, 0.97)) {
    x <- x[is.finite(x)]
    if (length(x) < 2) {
      return(seq(0, 1, length.out = n))
    }
    q <- stats::quantile(x, probs = probs, na.rm = TRUE)
    if (!is.finite(q[1]) || !is.finite(q[2]) || q[1] >= q[2]) {
      q <- range(x, na.rm = TRUE)
      if (!is.finite(q[1]) || !is.finite(q[2]) || q[1] == q[2]) {
        q <- c(0, 1)
      }
    }
    seq(q[1], q[2], length.out = n)
  }
  
  # Helper: bandwidth selection
  default_bandwidth <- function(x, method = "SJ", multiplier = 1.0) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(NA_real_)
    
    bw <- tryCatch({
      if (method == "SJ") {
        stats::bw.SJ(x, method = "ste")
      } else {
        stats::bw.nrd0(x)
      }
    }, error = function(e) NA_real_)
    
    if (!is.finite(bw) || bw <= 0) {
      iqr <- stats::IQR(x, na.rm = TRUE)
      bw <- if (is.finite(iqr) && iqr > 0) iqr / 1.34 else sd(x, na.rm = TRUE)
      if (!is.finite(bw) || bw <= 0) bw <- 1
    }
    
    multiplier * bw
  }
  
  # Helper: discretize KDE
  discretize_kde <- function(x, evalpoints, bw = NULL, bw_method = "SJ", bw_multiplier = 1.0) {
    x <- x[is.finite(x)]
    if (length(x) < 2) {
      return(rep(0, length(evalpoints)))
    }
    
    if (is.null(bw) || !is.finite(bw) || bw <= 0) {
      bw <- default_bandwidth(x, method = bw_method, multiplier = bw_multiplier)
    }
    
    dens <- stats::density(x, bw = bw, 
                          from = min(evalpoints), 
                          to = max(evalpoints),
                          n = length(evalpoints), 
                          kernel = "gaussian")
    dens$y
  }
  
  # Build evaluation grids (one per feature)
  grids <- list()
  for (var in age_vars) {
    all_vals <- ages_df[[var]]
    grids[[var]] <- get_evalpoints(all_vals, n = n_points, probs = probs)
  }
  
  # Build DZ matrix
  DZ_mat <- matrix(0, nrow = N, ncol = F * n_points)
  colnames_vec <- character(F * n_points)
  
  col_idx <- 1
  for (f in seq_along(age_vars)) {
    var <- age_vars[f]
    grid <- grids[[var]]
    
    # Compute bandwidth on all data (for consistency)
    all_vals <- ages_df[[var]]
    global_bw <- default_bandwidth(all_vals, method = bw_method, multiplier = bw_multiplier)
    
    # Discretize each sample
    for (i in seq_along(samples)) {
      samp <- samples[i]
      vals <- ages_df[[var]][ages_df$sample == samp]
      
      y <- discretize_kde(vals, grid, bw = global_bw, 
                         bw_method = bw_method, 
                         bw_multiplier = bw_multiplier)
      
      y[!is.finite(y)] <- 0
      
      # Normalize to unit area (trapezoidal rule)
      dx <- diff(grid)
      area <- sum((y[-1] + y[-length(y)]) * 0.5 * dx)
      if (is.finite(area) && area > 0) {
        y <- y / area
      }
      
      DZ_mat[i, col_idx:(col_idx + n_points - 1)] <- y
    }
    
    # Column names
    colnames_vec[col_idx:(col_idx + n_points - 1)] <- paste0(var, "_", 1:n_points)
    col_idx <- col_idx + n_points
  }
  
  rownames(DZ_mat) <- samples
  colnames(DZ_mat) <- colnames_vec
  
  # Return
  list(
    DZ_mat = DZ_mat,
    samples = samples,
    age_vars = age_vars,
    grids = grids,
    n_points = n_points,
    bw_method = bw_method,
    bw_multiplier = bw_multiplier
  )
}

#' Build Compositional Block from Bulk Petrology Counts
#' 
#' Converts raw mineral counts into compositional data ready for unmixing.
#' Optionally applies centered log-ratio (CLR) transformation or leaves in
#' simplex space (proportions).
#' 
#' @param counts_df Data frame with columns: sample (character/factor) and
#'   one or more numeric columns containing mineral counts or categorical
#'   mineral names. If categorical (e.g., "Qtz", "Fsp"), function will count
#'   occurrences per sample.
#' @param mineral_cols Character vector of column names to use as minerals.
#'   If NULL, uses all columns except 'sample'.
#' @param apply_clr Logical, apply CLR transformation? (default: FALSE).
#'   If FALSE, returns simplex (proportions). If TRUE, returns CLR coordinates.
#' @param impute Small value to add before CLR to avoid log(0) (default: 1e-6).
#'   Only used if apply_clr = TRUE.
#' @param samples_order Optional character vector specifying sample order.
#' @param min_count_per_sample Minimum total count per sample (default: 10).
#'   Samples with fewer counts are removed with a warning.
#' 
#' @return A list with components:
#'   \item{BP_mat}{Matrix (N x D) of compositional data (proportions if apply_clr=FALSE,
#'     CLR coordinates if apply_clr=TRUE)}
#'   \item{samples}{Character vector of sample names}
#'   \item{minerals}{Character vector of mineral names}
#'   \item{total_counts}{Integer vector of total counts per sample}
#'   \item{is_clr}{Logical, was CLR applied?}
#' 
#' @details
#' If counts_df has one row per grain (long format), function will aggregate by sample.
#' If counts_df has one row per sample with counts in columns (wide format), uses directly.
#' 
#' CLR transformation: clr(x) = log(x) - mean(log(x))
#' Converts from D-1 dimensional simplex to D-dimensional Euclidean space.
#' 
#' @examples
#' \dontrun{
#' # Example 1: Long format (one row per grain)
#' bp_raw <- data.frame(
#'   sample = rep(c("S1", "S2", "S3"), c(200, 180, 220)),
#'   mineral = sample(c("Qtz", "Fsp", "Lith", "Musc"), 600, replace = TRUE,
#'                   prob = c(0.4, 0.3, 0.2, 0.1))
#' )
#' 
#' bp_block <- build_bp_block(
#'   counts_df = bp_raw,
#'   mineral_cols = "mineral",  # column with mineral names
#'   apply_clr = FALSE  # keep as proportions
#' )
#' 
#' # Example 2: Wide format (one row per sample, counts in columns)
#' bp_wide <- data.frame(
#'   sample = c("S1", "S2", "S3"),
#'   Qtz = c(80, 70, 90),
#'   Fsp = c(60, 55, 65),
#'   Lith = c(40, 35, 45),
#'   Musc = c(20, 20, 20)
#' )
#' 
#' bp_block <- build_bp_block(
#'   counts_df = bp_wide,
#'   mineral_cols = c("Qtz", "Fsp", "Lith", "Musc"),
#'   apply_clr = TRUE  # transform to CLR
#' )
#' 
#' # Use in provenance_unmix
#' result <- provenance_unmix(
#'   data_list = list(DZ = dz_mat, BP = bp_block$BP_mat),
#'   data_types = c(DZ = "continuous", 
#'                  BP = ifelse(bp_block$is_clr, "clr", "compositional")),
#'   K = 3
#' )
#' }
#' 
#' @export
#' @importFrom compositions clr acomp
build_bp_block <- function(counts_df,
                            mineral_cols = NULL,
                            apply_clr = FALSE,
                            impute = 1e-6,
                            samples_order = NULL,
                            min_count_per_sample = 10) {
  
  # Validate input
  if (!"sample" %in% names(counts_df)) {
    stop("counts_df must have a 'sample' column")
  }
  
  # Determine mineral columns
  if (is.null(mineral_cols)) {
    mineral_cols <- setdiff(names(counts_df), "sample")
    message("Using mineral columns: ", paste(mineral_cols, collapse = ", "))
  }
  
  # Check if long format (one column with mineral names) or wide format (counts in columns)
  if (length(mineral_cols) == 1 && is.character(counts_df[[mineral_cols[1]]]) || 
      is.factor(counts_df[[mineral_cols[1]]])) {
    
    # Long format: one row per grain
    message("Detected long format (one row per grain). Aggregating by sample...")
    
    mineral_col_name <- mineral_cols[1]
    
    # Count by sample and mineral
    counts_wide <- counts_df %>%
      dplyr::group_by(sample, !!rlang::sym(mineral_col_name)) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(mineral_col_name), 
                        values_from = n, 
                        values_fill = 0)
    
    counts_df <- as.data.frame(counts_wide)
    mineral_cols <- setdiff(names(counts_df), "sample")
    
  } else {
    # Wide format: already have counts
    message("Detected wide format (counts in columns).")
  }
  
  # Extract samples and counts
  samples <- counts_df$sample
  counts_mat <- as.matrix(counts_df[, mineral_cols, drop = FALSE])
  rownames(counts_mat) <- samples
  
  # Check minimum counts
  total_counts <- rowSums(counts_mat)
  low_count_samples <- samples[total_counts < min_count_per_sample]
  
  if (length(low_count_samples) > 0) {
    warning(sprintf("Removing %d samples with < %d counts: %s",
                   length(low_count_samples), 
                   min_count_per_sample,
                   paste(low_count_samples, collapse = ", ")))
    
    keep_idx <- total_counts >= min_count_per_sample
    counts_mat <- counts_mat[keep_idx, , drop = FALSE]
    samples <- samples[keep_idx]
    total_counts <- total_counts[keep_idx]
  }
  
  # Apply sample order if specified
  if (!is.null(samples_order)) {
    if (!all(samples_order %in% samples)) {
      stop("Some samples in samples_order not found in data")
    }
    idx <- match(samples_order, samples)
    counts_mat <- counts_mat[idx, , drop = FALSE]
    samples <- samples_order
    total_counts <- total_counts[idx]
  }
  
  # Convert to proportions (close to simplex)
  props_mat <- counts_mat / rowSums(counts_mat)
  
  # Apply CLR if requested
  if (apply_clr) {
    # Add pseudocount
    props_mat_pseudo <- props_mat + impute
    props_mat_pseudo <- props_mat_pseudo / rowSums(props_mat_pseudo)
    
    # CLR transform
    BP_mat <- as.matrix(compositions::clr(compositions::acomp(props_mat_pseudo)))
    rownames(BP_mat) <- samples
    colnames(BP_mat) <- mineral_cols
  } else {
    BP_mat <- props_mat
  }
  
  # Return
  list(
    BP_mat = BP_mat,
    samples = samples,
    minerals = mineral_cols,
    total_counts = total_counts,
    is_clr = apply_clr
  )
}

#' Prepare Multiple Data Blocks for Unmixing
#' 
#' Convenience function to prepare multiple data blocks at once.
#' Combines build_dz_kde_block() and build_bp_block() with alignment checking.
#' 
#' @param dz_raw Data frame with detrital zircon ages (optional)
#' @param bp_raw Data frame with bulk petrology counts (optional)
#' @param hm_raw Data frame with heavy mineral counts (optional)
#' @param other_raw Named list of additional data frames (optional)
#' @param dz_vars Character vector of age variables for DZ
#' @param bp_vars Character vector of mineral columns for BP
#' @param hm_vars Character vector of mineral columns for HM
#' @param n_points_dz Number of bins for DZ discretization (default: 129)
#' @param apply_clr_bp Apply CLR to BP? (default: FALSE, will be done in unmixing)
#' @param apply_clr_hm Apply CLR to HM? (default: FALSE)
#' @param verbose Print progress messages? (default: TRUE)
#' 
#' @return A list ready for provenance_unmix():
#'   \item{data_list}{Named list of data matrices}
#'   \item{data_types}{Named character vector of data types}
#'   \item{samples}{Common sample names (aligned)}
#'   \item{metadata}{List with preprocessing details for each block}
#' 
#' @examples
#' \dontrun{
#' # Prepare all data at once
#' prepared <- prepare_data_blocks(
#'   dz_raw = my_dz_ages,
#'   bp_raw = my_bp_counts,
#'   dz_vars = c("age_concordia", "ti_temp"),
#'   bp_vars = c("Qtz", "Fsp", "Lith", "Musc"),
#'   n_points_dz = 129
#' )
#' 
#' # Unmix directly
#' result <- provenance_unmix(
#'   data_list = prepared$data_list,
#'   data_types = prepared$data_types,
#'   K = 3
#' )
#' }
#' 
#' @export
prepare_data_blocks <- function(dz_raw = NULL,
                                bp_raw = NULL,
                                hm_raw = NULL,
                                other_raw = NULL,
                                dz_vars = NULL,
                                bp_vars = NULL,
                                hm_vars = NULL,
                                n_points_dz = 129,
                                apply_clr_bp = FALSE,
                                apply_clr_hm = FALSE,
                                verbose = TRUE) {
  
  data_list <- list()
  data_types <- character()
  metadata_list <- list()
  
  # Process DZ
  if (!is.null(dz_raw)) {
    if (verbose) message("Processing detrital zircon data...")
    
    dz_block <- build_dz_kde_block(
      ages_df = dz_raw,
      age_vars = dz_vars,
      n_points = n_points_dz
    )
    
    data_list[["DZ"]] <- dz_block$DZ_mat
    data_types["DZ"] <- "continuous"
    metadata_list[["DZ"]] <- dz_block
  }
  
  # Process BP
  if (!is.null(bp_raw)) {
    if (verbose) message("Processing bulk petrology data...")
    
    bp_block <- build_bp_block(
      counts_df = bp_raw,
      mineral_cols = bp_vars,
      apply_clr = apply_clr_bp
    )
    
    data_list[["BP"]] <- bp_block$BP_mat
    data_types["BP"] <- ifelse(apply_clr_bp, "clr", "compositional")
    metadata_list[["BP"]] <- bp_block
  }
  
  # Process HM
  if (!is.null(hm_raw)) {
    if (verbose) message("Processing heavy mineral data...")
    
    hm_block <- build_bp_block(
      counts_df = hm_raw,
      mineral_cols = hm_vars,
      apply_clr = apply_clr_hm
    )
    
    data_list[["HM"]] <- hm_block$BP_mat
    data_types["HM"] <- ifelse(apply_clr_hm, "clr", "compositional")
    metadata_list[["HM"]] <- hm_block
  }
  
  # Process other blocks
  if (!is.null(other_raw) && length(other_raw) > 0) {
    for (name in names(other_raw)) {
      if (verbose) message(sprintf("Processing %s data...", name))
      # Assume continuous for now (user can override)
      data_list[[name]] <- as.matrix(other_raw[[name]])
      data_types[name] <- "continuous"
    }
  }
  
  # Check alignment
  if (length(data_list) > 0) {
    all_samples <- lapply(data_list, rownames)
    common_samples <- Reduce(intersect, all_samples)
    
    if (length(common_samples) == 0) {
      stop("No common samples found across all data blocks!")
    }
    
    n_samples_per_block <- sapply(all_samples, length)
    if (any(n_samples_per_block != length(common_samples))) {
      warning(sprintf(
        "Not all samples present in all blocks. Using %d common samples.",
        length(common_samples)
      ))
      
      # Subset to common samples
      for (i in seq_along(data_list)) {
        data_list[[i]] <- data_list[[i]][common_samples, , drop = FALSE]
      }
    }
    
    if (verbose) {
      message(sprintf("Data preparation complete: %d samples, %d blocks",
                     length(common_samples), length(data_list)))
    }
    
    return(list(
      data_list = data_list,
      data_types = data_types,
      samples = common_samples,
      metadata = metadata_list
    ))
    
  } else {
    stop("No data blocks provided!")
  }
}
