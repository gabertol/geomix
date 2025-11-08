#' Data Preprocessing Functions for Source Analysis
#'
#' Functions to prepare raw data for source unmixing, including
#' KDE discretization for detrital zircon ages and compositional data
#' preparation for bulk petrology.
#' 
#' @name preprocessing
NULL

#' Build KDE-Discretized Block from Continuous Distributional Data
#'
#' Converts continuous distributional measurements into discretized probability
#' density functions (PDFs) via kernel density estimation (KDE). This is the
#' recommended format for the "continuous" data type in unmix().
#'
#' @param data_df Data frame with columns: sample (character/factor) and
#'   one or more numeric columns containing continuous measurements (e.g.,
#'   detrital zircon ages, apatite fission track ages, trace element ratios, etc.).
#'   Each row is one grain/measurement.
#' @param vars Character vector of column names to discretize. If NULL,
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
#'   \item{data_mat}{Matrix (N x F*n_points) of discretized PDFs, where N is
#'     number of samples and F is number of features (vars)}
#'   \item{samples}{Character vector of sample names (row names of data_mat)}
#'   \item{vars}{Character vector of feature names}
#'   \item{grids}{Named list of evaluation grids (one per feature)}
#'   \item{n_points}{Number of bins per feature}
#'   \item{data_type}{Character, type for use in unmix ("continuous")}
#'
#' @details
#' For each sample and each variable:
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
#' # Example 1: Detrital zircon ages
#' dz_raw <- data.frame(
#'   sample = rep(c("S1", "S2", "S3"), each = 100),
#'   age_concordia = c(rnorm(100, 500, 50), rnorm(100, 1000, 100), rnorm(100, 1500, 150)),
#'   ti_temp = rnorm(300, 750, 50),
#'   th_u = rlnorm(300, 0, 0.5)
#' )
#'
#' kde_block <- build_kde_block(
#'   data_df = dz_raw,
#'   vars = c("age_concordia", "ti_temp", "th_u"),
#'   n_points = 129
#' )
#'
#' # Example 2: Apatite fission track ages
#' apt_raw <- data.frame(
#'   sample = rep(c("S1", "S2", "S3"), each = 80),
#'   age_ft = rnorm(240, 100, 20),
#'   age_u_pb = rnorm(240, 500, 50)
#' )
#'
#' apt_block <- build_kde_block(
#'   data_df = apt_raw,
#'   vars = c("age_ft", "age_u_pb")
#' )
#'
#' # Use in unmix
#' result <- unmix(
#'   data_list = list(DZ = kde_block$data_mat, APT = apt_block$data_mat),
#'   data_types = c(DZ = "continuous", APT = "continuous"),
#'   K = 3
#' )
#' }
#'
#' @export
build_kde_block <- function(data_df,
                            vars = NULL,
                            n_points = 129,
                            probs = c(0.03, 0.97),
                            bw_method = c("SJ", "nrd0"),
                            bw_multiplier = 1.0,
                            samples_order = NULL) {

  bw_method <- match.arg(bw_method)

  # Validate input
  if (!"sample" %in% names(data_df)) {
    stop("data_df must have a 'sample' column")
  }

  # Determine variables
  if (is.null(vars)) {
    vars <- setdiff(names(data_df), "sample")
    vars <- vars[sapply(data_df[vars], is.numeric)]
    if (length(vars) == 0) {
      stop("No numeric columns found (besides 'sample')")
    }
    message("Using variables: ", paste(vars, collapse = ", "))
  }

  # Check all vars exist
  missing_vars <- setdiff(vars, names(data_df))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data_df: ", paste(missing_vars, collapse = ", "))
  }

  # Get unique samples
  samples <- unique(data_df$sample)
  if (!is.null(samples_order)) {
    if (!all(samples_order %in% samples)) {
      stop("Some samples in samples_order not found in data")
    }
    samples <- samples_order
  } else {
    samples <- sort(samples)
  }
  
  N <- length(samples)
  F <- length(vars)

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
  for (var in vars) {
    all_vals <- data_df[[var]]
    grids[[var]] <- get_evalpoints(all_vals, n = n_points, probs = probs)
  }

  # Build data matrix
  data_mat <- matrix(0, nrow = N, ncol = F * n_points)
  colnames_vec <- character(F * n_points)

  col_idx <- 1
  for (f in seq_along(vars)) {
    var <- vars[f]
    grid <- grids[[var]]

    # Compute bandwidth on all data (for consistency)
    all_vals <- data_df[[var]]
    global_bw <- default_bandwidth(all_vals, method = bw_method, multiplier = bw_multiplier)

    # Discretize each sample
    for (i in seq_along(samples)) {
      samp <- samples[i]
      vals <- data_df[[var]][data_df$sample == samp]

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

      data_mat[i, col_idx:(col_idx + n_points - 1)] <- y
    }

    # Column names
    colnames_vec[col_idx:(col_idx + n_points - 1)] <- paste0(var, "_", 1:n_points)
    col_idx <- col_idx + n_points
  }

  rownames(data_mat) <- samples
  colnames(data_mat) <- colnames_vec

  # Return
  list(
    data_mat = data_mat,
    samples = samples,
    vars = vars,
    grids = grids,
    n_points = n_points,
    bw_method = bw_method,
    bw_multiplier = bw_multiplier,
    data_type = "continuous"
  )
}

#' Build Point Counting Block from Discrete Counts
#'
#' Converts raw point counting data (discrete mineral/component counts) into
#' compositional data ready for unmixing. This function is specifically for
#' COUNTING data (integers representing number of grains counted), not
#' data that are already in proportion/composition form.
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
#'   \item{data_mat}{Matrix (N x D) of compositional data (proportions if apply_clr=FALSE,
#'     CLR coordinates if apply_clr=TRUE)}
#'   \item{samples}{Character vector of sample names}
#'   \item{components}{Character vector of component names}
#'   \item{total_counts}{Integer vector of total counts per sample}
#'   \item{is_clr}{Logical, was CLR applied?}
#'   \item{data_type}{Character, type for use in unmix}
#'
#' @details
#' **Input Format:**
#'
#' Point counting data are DISCRETE COUNTS (integers) from counting exercises.
#' Examples: counted 100 grains, 40 were Quartz, 30 Feldspar, 30 Lithics.
#'
#' If counts_df has one row per grain (long format), function will aggregate by sample.
#' If counts_df has one row per sample with counts in columns (wide format), uses directly.
#'
#' This function will:
#' 1. Aggregate counts if needed
#' 2. Filter samples with low counts
#' 3. Convert to proportions (closure to simplex: divide by row sum)
#' 4. Optionally apply CLR transformation
#'
#' CLR transformation: clr(x) = log(x) - mean(log(x))
#' Converts from D-1 dimensional simplex to D-dimensional Euclidean space.
#'
#' @examples
#' \dontrun{
#' # Example 1: Long format (one row per grain counted)
#' bp_raw <- data.frame(
#'   sample = rep(c("S1", "S2", "S3"), c(200, 180, 220)),
#'   mineral = sample(c("Qtz", "Fsp", "Lith", "Musc"), 600, replace = TRUE,
#'                   prob = c(0.4, 0.3, 0.2, 0.1))
#' )
#'
#' bp_block <- build_point_counting_block(
#'   counts_df = bp_raw,
#'   mineral_cols = "mineral",  # column with mineral names
#'   apply_clr = FALSE  # keep as proportions
#' )
#'
#' # Example 2: Wide format (one row per sample, counts in columns)
#' bp_wide <- data.frame(
#'   sample = c("S1", "S2", "S3"),
#'   Qtz = c(80, 70, 90),      # these are COUNTS (integers)
#'   Fsp = c(60, 55, 65),
#'   Lith = c(40, 35, 45),
#'   Musc = c(20, 20, 20)
#' )
#'
#' bp_block <- build_point_counting_block(
#'   counts_df = bp_wide,
#'   mineral_cols = c("Qtz", "Fsp", "Lith", "Musc"),
#'   apply_clr = TRUE  # transform to CLR
#' )
#'
#' # Use in unmix
#' result <- unmix(
#'   data_list = list(DZ = dz_mat, BP = bp_block$data_mat),
#'   data_types = c(DZ = "continuous", BP = bp_block$data_type),
#'   K = 3
#' )
#' }
#'
#' @export
#' @importFrom compositions clr acomp
build_point_counting_block <- function(counts_df,
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
    data_mat = BP_mat,
    samples = samples,
    components = mineral_cols,
    total_counts = total_counts,
    is_clr = apply_clr,
    data_type = ifelse(apply_clr, "clr", "compositional")
  )
}

#' Build Compositional Block from Data Already in Compositional Form
#'
#' Processes data that are ALREADY compositional (i.e., proportions/compositions
#' that sum to 1) for use in source unmixing. This is for data that are
#' already closed to the simplex, NOT for raw counts.
#'
#' @param comp_df Data frame with columns: sample (character/factor) and
#'   one or more numeric columns containing compositional data (proportions).
#'   Values should sum to ~1.0 for each sample (small deviations are normalized).
#' @param comp_cols Character vector of column names to use as components.
#'   If NULL, uses all numeric columns except 'sample'.
#' @param apply_clr Logical, apply CLR transformation? (default: FALSE).
#'   If FALSE, returns simplex (proportions). If TRUE, returns CLR coordinates.
#' @param impute Small value to add before CLR to avoid log(0) (default: 1e-6).
#'   Only used if apply_clr = TRUE.
#' @param samples_order Optional character vector specifying sample order.
#' @param check_closure Logical, check if rows sum to ~1? (default: TRUE).
#'   If TRUE and row sums deviate significantly from 1, issues a warning.
#' @param closure_tolerance Tolerance for closure check (default: 0.01).
#'   Rows must sum to 1 ± tolerance.
#'
#' @return A list with components:
#'   \item{data_mat}{Matrix (N x D) of compositional data (proportions if apply_clr=FALSE,
#'     CLR coordinates if apply_clr=TRUE)}
#'   \item{samples}{Character vector of sample names}
#'   \item{components}{Character vector of component names}
#'   \item{row_sums}{Numeric vector of row sums (should be ~1.0)}
#'   \item{is_clr}{Logical, was CLR applied?}
#'   \item{data_type}{Character, type for use in unmix}
#'
#' @details
#' **Input Format:**
#'
#' Compositional data are already in PROPORTION form (values between 0 and 1, sum to 1).
#' Examples:
#' - Modal mineralogy from XRD (already normalized to 100%)
#' - Geochemical oxides (e.g., SiO2=0.65, Al2O3=0.15, ..., sum=1.00)
#' - Any data already closed to the simplex
#'
#' This function will:
#' 1. Validate that data are in [0, 1] range
#' 2. Optionally check that rows sum to ~1.0
#' 3. Re-normalize to simplex (closure) to ensure exact sum=1
#' 4. Optionally apply CLR transformation
#'
#' **Difference from build_point_counting_block():**
#'
#' - `build_point_counting_block()`: Input = COUNTS (integers: 40, 30, 30, ...)
#' - `build_compositional_block()`: Input = PROPORTIONS (decimals: 0.40, 0.30, 0.30, ...)
#'
#' @examples
#' \dontrun{
#' # Example: Modal mineralogy data (already proportions)
#' modal_data <- data.frame(
#'   sample = c("S1", "S2", "S3"),
#'   Qtz = c(0.40, 0.35, 0.45),     # these are PROPORTIONS (sum to 1)
#'   Fsp = c(0.30, 0.35, 0.30),
#'   Lith = c(0.20, 0.20, 0.15),
#'   Musc = c(0.10, 0.10, 0.10)
#' )
#'
#' comp_block <- build_compositional_block(
#'   comp_df = modal_data,
#'   comp_cols = c("Qtz", "Fsp", "Lith", "Musc"),
#'   apply_clr = TRUE,
#'   check_closure = TRUE
#' )
#'
#' # Use in unmix
#' result <- unmix(
#'   data_list = list(DZ = dz_mat, Modal = comp_block$data_mat),
#'   data_types = c(DZ = "continuous", Modal = comp_block$data_type),
#'   K = 3
#' )
#' }
#'
#' @export
#' @importFrom compositions clr acomp
build_compositional_block <- function(comp_df,
                                       comp_cols = NULL,
                                       apply_clr = FALSE,
                                       impute = 1e-6,
                                       samples_order = NULL,
                                       check_closure = TRUE,
                                       closure_tolerance = 0.01) {

  # Validate input
  if (!"sample" %in% names(comp_df)) {
    stop("comp_df must have a 'sample' column")
  }

  # Determine compositional columns
  if (is.null(comp_cols)) {
    comp_cols <- setdiff(names(comp_df), "sample")
    comp_cols <- comp_cols[sapply(comp_df[comp_cols], is.numeric)]
    if (length(comp_cols) == 0) {
      stop("No numeric columns found (besides 'sample')")
    }
    message("Using compositional columns: ", paste(comp_cols, collapse = ", "))
  }

  # Extract samples and data
  samples <- comp_df$sample
  comp_mat <- as.matrix(comp_df[, comp_cols, drop = FALSE])
  rownames(comp_mat) <- samples

  # Validate compositional data
  if (any(comp_mat < 0, na.rm = TRUE)) {
    stop("Compositional data contains negative values. Check your input.")
  }

  if (any(comp_mat > 1.1, na.rm = TRUE)) {
    warning("Some values > 1.1 detected. Are you sure this is compositional data (not counts)?")
  }

  # Check closure (row sums should be ~1)
  row_sums <- rowSums(comp_mat)
  if (check_closure) {
    deviations <- abs(row_sums - 1.0)
    if (any(deviations > closure_tolerance, na.rm = TRUE)) {
      max_dev <- max(deviations, na.rm = TRUE)
      warning(sprintf(
        "Some samples deviate from closure (sum != 1). Max deviation: %.4f. Re-normalizing to simplex.",
        max_dev
      ))
    }
  }

  # Apply sample order if specified
  if (!is.null(samples_order)) {
    if (!all(samples_order %in% samples)) {
      stop("Some samples in samples_order not found in data")
    }
    idx <- match(samples_order, samples)
    comp_mat <- comp_mat[idx, , drop = FALSE]
    samples <- samples_order
    row_sums <- row_sums[idx]
  }

  # Close to simplex (ensure exact sum = 1)
  comp_mat <- comp_mat / rowSums(comp_mat)

  # Apply CLR if requested
  if (apply_clr) {
    # Add pseudocount
    comp_mat_pseudo <- comp_mat + impute
    comp_mat_pseudo <- comp_mat_pseudo / rowSums(comp_mat_pseudo)

    # CLR transform
    result_mat <- as.matrix(compositions::clr(compositions::acomp(comp_mat_pseudo)))
    rownames(result_mat) <- samples
    colnames(result_mat) <- comp_cols
  } else {
    result_mat <- comp_mat
  }

  # Return
  list(
    data_mat = result_mat,
    samples = samples,
    components = comp_cols,
    row_sums = row_sums,
    is_clr = apply_clr,
    data_type = ifelse(apply_clr, "clr", "compositional")
  )
}

#' Prepare Multiple Data Blocks for Unmixing
#'
#' Convenience function to prepare multiple data blocks at once.
#' Combines KDE, point counting, and compositional data preparation with
#' alignment checking. Supports multiple blocks of the same type via named lists.
#'
#' @param kde_raw Data frame or named list of data frames with KDE data
#'   (e.g., detrital zircon ages). Optional.
#' @param point_counting_raw Data frame or named list of data frames with
#'   point counting data (discrete counts). Optional.
#' @param compositional_raw Data frame or named list of data frames with
#'   compositional data (already proportions). Optional.
#' @param other_raw Named list of additional data frames (optional)
#' @param kde_vars Character vector or named list of character vectors for KDE variables.
#'   If a single vector, applies to all KDE blocks. If named list, must match names in kde_raw.
#' @param point_counting_vars Character vector or named list for point counting variables.
#' @param compositional_vars Character vector or named list for compositional variables.
#' @param n_points_kde Number of bins for KDE discretization (default: 129).
#'   Can be a single value or named list matching kde_raw.
#' @param apply_clr_point_counting Apply CLR to point counting? (default: FALSE).
#'   Can be single logical or named list.
#' @param apply_clr_compositional Apply CLR to compositional? (default: FALSE).
#'   Can be single logical or named list.
#' @param verbose Print progress messages? (default: TRUE)
#' 
#' @return A list ready for unmix():
#'   \item{data_list}{Named list of data matrices}
#'   \item{data_types}{Named character vector of data types}
#'   \item{samples}{Common sample names (aligned)}
#'   \item{metadata}{List with preprocessing details for each block}
#' 
#' @examples
#' \dontrun{
#' # Example 1: Single block per type
#' prepared <- prepare_data_blocks(
#'   kde_raw = my_dz_ages,
#'   point_counting_raw = my_mineral_counts,
#'   kde_vars = c("age_concordia", "ti_temp"),
#'   point_counting_vars = c("Qtz", "Fsp", "Lith", "Musc"),
#'   n_points_kde = 129
#' )
#'
#' # Example 2: Multiple KDE blocks
#' prepared <- prepare_data_blocks(
#'   kde_raw = list(
#'     Zircon = zircon_ages,
#'     Apatite = apatite_ages
#'   ),
#'   kde_vars = list(
#'     Zircon = c("age_concordia", "ti_temp"),
#'     Apatite = c("age_u_pb", "age_fission_track")
#'   )
#' )
#'
#' # Example 3: Multiple compositional blocks
#' prepared <- prepare_data_blocks(
#'   compositional_raw = list(
#'     MajorOxides = major_oxides,
#'     TraceElements = trace_elements,
#'     ModalMin = modal_mineralogy
#'   ),
#'   compositional_vars = list(
#'     MajorOxides = c("SiO2", "Al2O3", "FeO", "MgO"),
#'     TraceElements = c("Zr", "Y", "Nb", "La", "Ce"),
#'     ModalMin = c("Qtz", "Fsp", "Mica", "Amph")
#'   ),
#'   apply_clr_compositional = TRUE
#' )
#'
#' # Example 4: Mixed - 2 KDE + 3 Compositional
#' prepared <- prepare_data_blocks(
#'   kde_raw = list(
#'     Zircon = zircon_ages,
#'     Apatite = apatite_ages
#'   ),
#'   compositional_raw = list(
#'     MajorOxides = major_oxides,
#'     TraceElements = trace_elements,
#'     ModalMin = modal_mineralogy
#'   ),
#'   kde_vars = list(
#'     Zircon = c("age_concordia"),
#'     Apatite = c("age_u_pb")
#'   ),
#'   compositional_vars = list(
#'     MajorOxides = c("SiO2", "Al2O3", "FeO"),
#'     TraceElements = c("Zr", "Y", "Nb"),
#'     ModalMin = c("Qtz", "Fsp", "Mica")
#'   )
#' )
#'
#' # Unmix directly
#' result <- unmix(
#'   data_list = prepared$data_list,
#'   data_types = prepared$data_types,
#'   K = 3
#' )
#' }
#' 
#' @export
prepare_data_blocks <- function(kde_raw = NULL,
                                point_counting_raw = NULL,
                                compositional_raw = NULL,
                                other_raw = NULL,
                                kde_vars = NULL,
                                point_counting_vars = NULL,
                                compositional_vars = NULL,
                                n_points_kde = 129,
                                apply_clr_point_counting = FALSE,
                                apply_clr_compositional = FALSE,
                                verbose = TRUE) {

  data_list <- list()
  data_types <- character()
  metadata_list <- list()

  # Helper function to get parameter for specific block
  get_param <- function(param, block_name, default = NULL) {
    if (is.null(param)) return(default)
    if (is.list(param) && !is.data.frame(param)) {
      if (block_name %in% names(param)) {
        return(param[[block_name]])
      } else {
        stop(sprintf("Block '%s' not found in parameter list", block_name))
      }
    } else {
      return(param)
    }
  }

  # Process KDE data
  if (!is.null(kde_raw)) {
    # Check if single data.frame or list of data.frames
    if (is.data.frame(kde_raw)) {
      # Single block
      if (verbose) message("Processing KDE data...")

      kde_block <- build_kde_block(
        data_df = kde_raw,
        vars = kde_vars,
        n_points = n_points_kde
      )

      data_list[["KDE"]] <- kde_block$data_mat
      data_types["KDE"] <- "continuous"
      metadata_list[["KDE"]] <- kde_block

    } else if (is.list(kde_raw)) {
      # Multiple blocks
      if (is.null(names(kde_raw))) {
        stop("kde_raw list must have names")
      }

      for (block_name in names(kde_raw)) {
        if (verbose) message(sprintf("Processing KDE data: %s...", block_name))

        kde_block <- build_kde_block(
          data_df = kde_raw[[block_name]],
          vars = get_param(kde_vars, block_name, NULL),
          n_points = get_param(n_points_kde, block_name, 129)
        )

        data_list[[block_name]] <- kde_block$data_mat
        data_types[block_name] <- "continuous"
        metadata_list[[block_name]] <- kde_block
      }
    } else {
      stop("kde_raw must be a data.frame or a named list of data.frames")
    }
  }

  # Process Point Counting data
  if (!is.null(point_counting_raw)) {
    # Check if single data.frame or list of data.frames
    if (is.data.frame(point_counting_raw)) {
      # Single block
      if (verbose) message("Processing point counting data (discrete counts)...")

      pc_block <- build_point_counting_block(
        counts_df = point_counting_raw,
        mineral_cols = point_counting_vars,
        apply_clr = apply_clr_point_counting
      )

      data_list[["PointCounting"]] <- pc_block$data_mat
      data_types["PointCounting"] <- pc_block$data_type
      metadata_list[["PointCounting"]] <- pc_block

    } else if (is.list(point_counting_raw)) {
      # Multiple blocks
      if (is.null(names(point_counting_raw))) {
        stop("point_counting_raw list must have names")
      }

      for (block_name in names(point_counting_raw)) {
        if (verbose) message(sprintf("Processing point counting data: %s...", block_name))

        pc_block <- build_point_counting_block(
          counts_df = point_counting_raw[[block_name]],
          mineral_cols = get_param(point_counting_vars, block_name, NULL),
          apply_clr = get_param(apply_clr_point_counting, block_name, FALSE)
        )

        data_list[[block_name]] <- pc_block$data_mat
        data_types[block_name] <- pc_block$data_type
        metadata_list[[block_name]] <- pc_block
      }
    } else {
      stop("point_counting_raw must be a data.frame or a named list of data.frames")
    }
  }

  # Process Compositional data
  if (!is.null(compositional_raw)) {
    # Check if single data.frame or list of data.frames
    if (is.data.frame(compositional_raw)) {
      # Single block
      if (verbose) message("Processing compositional data (already proportions)...")

      comp_block <- build_compositional_block(
        comp_df = compositional_raw,
        comp_cols = compositional_vars,
        apply_clr = apply_clr_compositional
      )

      data_list[["Compositional"]] <- comp_block$data_mat
      data_types["Compositional"] <- comp_block$data_type
      metadata_list[["Compositional"]] <- comp_block

    } else if (is.list(compositional_raw)) {
      # Multiple blocks
      if (is.null(names(compositional_raw))) {
        stop("compositional_raw list must have names")
      }

      for (block_name in names(compositional_raw)) {
        if (verbose) message(sprintf("Processing compositional data: %s...", block_name))

        comp_block <- build_compositional_block(
          comp_df = compositional_raw[[block_name]],
          comp_cols = get_param(compositional_vars, block_name, NULL),
          apply_clr = get_param(apply_clr_compositional, block_name, FALSE)
        )

        data_list[[block_name]] <- comp_block$data_mat
        data_types[block_name] <- comp_block$data_type
        metadata_list[[block_name]] <- comp_block
      }
    } else {
      stop("compositional_raw must be a data.frame or a named list of data.frames")
    }
  }

  # Process HM (deprecated, map to point counting)
  if (!is.null(hm_raw)) {
    if (verbose) message("Processing heavy mineral data (using point counting)...")

    hm_block <- build_point_counting_block(
      counts_df = hm_raw,
      mineral_cols = hm_vars,
      apply_clr = ifelse(is.null(apply_clr_hm), FALSE, apply_clr_hm)
    )

    data_list[["HM"]] <- hm_block$data_mat
    data_types["HM"] <- hm_block$data_type
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
