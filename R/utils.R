
# -------------------------------------------------------------------------
# Utility functions for quantifying and visualising spatial dispersion
# -------------------------------------------------------------------------


get_centroids <- function() {
  
  # Generates the reference centroid coordinates for the survey design.
  #
  # Arguments:
  #   None
  #
  # Returns:
  #   A matrix with 10 rows and 2 columns (x, y) representing segment centroids.
  #   The units are metres.
  
  # x coordinates: Two parallel transects at x=250 and x=750
  x_coords <- rep(c(250, 750), each = 5)
  
  # y coordinates: 5 segments of 200m each (midpoints: 100, 300, ..., 900)
  # Transect 1 goes up (100-900), Transect 2 goes down (900-100)
  y_seq <- seq(100, 900, by = 200)
  y_coords <- c(y_seq, rev(y_seq))
  
  cbind(x = x_coords, y = y_coords)
}


inv_obs <- function(counts, centroids = get_centroids()) {
  
  # Reverses the observation binning process by reassigning birds to their 
  # original spatial coverage areas.
  #
  # Arguments:
  #      counts: Integer vector of length 10 representing bird counts per segment.
  #   centroids: Matrix (10x2) of segment centers. Defaults to MPPL design.
  #
  # Returns:
  #   A numeric matrix (Nx2) of reconstructed point coordinates.
  #
  # Details:
  #   This function performs the inverse of the MPPL observation process. 
  #   In the field, observations are aggregated into 200m-long segments and 
  #   include birds detected within 100m bands on either side of the transect. 
  #   Therefore, a single count represents birds located anywhere within a 
  #   200m x 200m area centered on the segment centroid.
  #
  #   This function reconstructs the fine-scale distribution by uniformly 
  #   distributing the observed count 'n' across this specific 200m x 200m 
  #   spatial domain (Centroid +/- 100m in both X and Y dimensions).
  
  # Validate input
  if (length(counts) != nrow(centroids)) {
    stop("Length of 'counts' must match the number of rows in 'centroids'.")
  }
  
  # Vectorized expansion: Repeat segment indices based on counts
  indices <- rep(seq_along(counts), times = counts)
  
  # Map indices to reference coordinates
  xy <- centroids[indices, , drop = FALSE]
  
  # Reassign points to their specific survey area.
  # The uniform distribution reflects the assumption that a bird was equally 
  # likely to be anywhere within the specific 200x200m survey unit 
  # (Segment length: 200m; Total band width: 100m left + 100m right = 200m).
  n <- nrow(xy)
  if (n > 0) {
    xy[, 1] <- xy[, 1] + runif(n, min = -100, max = 100)
    xy[, 2] <- xy[, 2] + runif(n, min = -100, max = 100)
  }
  
  xy
}


nnd <- function(xy) {
  
  # Calculates the mean Nearest Neighbour Distance (NND) for a point pattern.
  #
  # Arguments:
  #   xy: A numeric matrix or data frame of coordinates.
  #
  # Returns:
  #   A single numeric value representing the mean distance to the nearest neighbour.
  #   Returns NA if fewer than 2 points are provided.
  
  xy <- as.matrix(xy)
  n <- nrow(xy)
  
  if (n < 2) return(NA_real_)
  
  # Calculate full distance matrix
  d_mat <- as.matrix(dist(xy))
  
  # Set diagonal to Inf so a point is not its own nearest neighbour
  diag(d_mat) <- Inf
  
  # Find min distance for each row (point) and take the mean
  mean(apply(d_mat, 1, min))
}


centroid <- function(xy) {
  
  # Calculates the geometric center (centroid) of a cloud of points.
  #
  # Arguments:
  #   xy: A numeric matrix or data frame of coordinates.
  #
  # Returns:
  #   A named numeric vector corresponding to the column means (x, y).
  
  colMeans(xy)
}


sdd <- function(xy, na.rm = FALSE) {
  
  # Calculates the Standard Distance Deviation (SDD) for a set of 2D points.
  #
  # Arguments:
  #      xy: Numeric matrix or data frame with exactly 2 columns.
  #   na.rm: Logical. If TRUE, removes rows with NA values.
  #
  # Returns:
  #   A numeric scalar representing the standard distance deviation.
  #
  # Details:
  #   SDD is the spatial equivalent of standard deviation. It measures the 
  #   degree to which points are concentrated around the mean center.
  
  # Input validation
  if (is.data.frame(xy)) xy <- as.matrix(xy)
  stopifnot(is.numeric(xy), ncol(xy) == 2, is.logical(na.rm))
  
  if (na.rm) xy <- xy[complete.cases(xy), , drop = FALSE]
  
  n <- nrow(xy)
  if (n < 2) return(0)
  
  # Efficiently center the data using scale()
  # scale(..., scale=FALSE) subtracts the column means
  centered <- scale(xy, center = TRUE, scale = FALSE)
  
  # Calculate SDD: sqrt( sum(dist^2) / n )
  # Summing the squares of the centered matrix is equivalent to sum of squared Euclidean distances
  sqrt(sum(centered^2) / n)
}


plot_obs <- function(xy, centroids = get_centroids()) {
  
  # Visualizes the spatial distribution of observations relative to survey segments.
  #
  # Arguments:
  #          xy: Numeric matrix of point coordinates.
  #   centroids: Matrix of segment centers (default: MPPL design).
  #
  # Returns:
  #   NULL (Invisibly). Generates a base R plot.
  
  # Save graphical parameters and restore on exit (scoping safety)
  op <- par(mar = c(2, 2, 2, 2))
  on.exit(par(op))
  
  n <- nrow(xy)
  
  # Plot points
  plot(xy, 
       xlim = c(0, 1000), 
       ylim = c(0, 1000), 
       asp = 1, 
       col = adjustcolor("blue", alpha.f = max(0.1, 20 / n)), # Prevent alpha > 1 or too small
       axes = FALSE, 
       xlab = "", 
       ylab = "")
  
  # Add context: Transect lines
  abline(v = c(250, 750), col = "red")
  
  # Add context: Segment markers
  seg_markers <- centroids
  seg_markers[, 2] <- seg_markers[, 2] + 100
  points(seg_markers, pch = "-", cex = 1.5, col = "red")
  
  # Plot boundary
  rect(0, 0, 1000, 1000)
  
  invisible(NULL)
}


csr <- function(n, n_segments = 10, n_sim = 1) {
  
  # Generates counts under Complete Spatial Randomness (CSR).
  #
  # Arguments:
  #            n: Total number of individuals to distribute.
  #   n_segments: Number of spatial bins (default 10).
  #        n_sim: Number of simulations to run.
  #
  # Returns:
  #   A matrix where each column is a simulation of counts per segment.
  
  rmultinom(n = n_sim, size = n, prob = rep(1 / n_segments, n_segments))
}


nsdi <- function(counts, n_sim = 100, metric_fn = sdd) {
  
  # Calculates the Normalised Spatial Dispersion Index (NSDI).
  #
  # Arguments:
  #      counts: Integer vector of observed counts per segment.
  #       n_sim: Integer. Number of iterations for reconstruction and CSR generation.
  #   metric_fn: Function to calculate spatial statistic (default: sdd).
  #
  # Returns:
  #   A numeric scalar: Ratio of observed metric to mean random metric.
  #
  # Details:
  #   NSDI compares the observed spatial dispersion to that expected under 
  #   Complete Spatial Randomness (CSR).
  #
  #   1. Observed Value: To calculate the metric for the observed data, the 
  #      function performs 'n_sim' iterations of the inverse observation process 
  #      (inv_obs). This repeatedly reassigns individuals to specific coordinates 
  #      within their detection zones, effectively integrating over the 
  #      uncertainty of the exact locations.
  #
  #   2. Expected Value: It generates 'n_sim' random count vectors (CSR) and 
  #      applies the same inverse mapping to them.
  #
  #   Interpretation:
  #     NSDI < 1: Clustering (aggregated distribution).
  #     NSDI > 1: Dispersion (regular distribution).
  #     NSDI ~ 1: Random distribution.
  
  # 1. Observed Statistic (Mean of reconstructed realizations)
  # We replicate the inverse mapping of the OBSERVED counts 'n_sim' times.
  # This averages out the variation inherent in reassigning birds within 
  # their 200x200m survey units.
  stat_obs_vec <- vapply(
    seq_len(n_sim), 
    function(i) metric_fn(inv_obs(counts)), 
    numeric(1)
  )
  mean_obs <- mean(stat_obs_vec)
  
  # 2. Random (CSR) Statistic
  # Generate random counts under the null hypothesis (equal probability per segment)
  total_n <- sum(counts)
  counts_rand_mat <- csr(n = total_n, n_segments = length(counts), n_sim = n_sim)
  
  # For each random count vector, apply the inverse observation mapping 
  # to convert counts to coordinates, then calculate the metric.
  stat_rand_vec <- apply(counts_rand_mat, 2, function(x) metric_fn(inv_obs(x)))
  mean_rand <- mean(stat_rand_vec)
  
  # Avoid division by zero
  if (mean_rand == 0) return(NA_real_)
  
  mean_obs / mean_rand
}


disp <- function(df, ...) {
  # Calculate spatial dispersion for every combination of 'plot_id' and 'year'.
  
  # Arguments:
  #    df: data.frame containing variables: id_year, id, year, observer_id, n
  #   ...: arguments passed to `nsdi()`
  #   
  # Returns: 
  #   data.frame with above mentioned variables, distribution of observations 
  #   across segments (n_1:n_10), and the variable 'dispersion'.
  
  yr <- sort(unique(df$year))
  id <- sort(unique(df$id))
  
  nd <- expand.grid(id = id, year = yr, segment = 1:10)
  nd <- transform(nd, id_year = interaction(id, year, sep = "_"))
  nd <- transform(nd, id_year_seg = interaction(id_year, segment, sep = "_"))
  nd <- merge(nd, df[c("id_year_seg", "n")], by = "id_year_seg", all.x = TRUE)
  nd$n[is.na(nd$n)] <- 0
  
  n_by_segment <- aggregate(n ~ id_year, nd, I)
  n_by_segment <- data.frame(id_year = n_by_segment$id_year, n_by_segment$n)
  names(n_by_segment)[-1] <- paste("n", 1:10, sep = "_")
  
  disp_by_idy <- aggregate(n ~ id_year, nd, nsdi, ...)
  # disp_by_idy <- aggregate(n ~ id_year, nd, nsdi)
  names(disp_by_idy)[2] <- "dispersion"
  
  n_by_idy <- aggregate(n ~ id_year + id + year + observer_id, df, sum)
  merge(merge(n_by_idy, n_by_segment), disp_by_idy)
}


term_plot <- function(model, formula, 
                      group_levels = NULL, 
                      transform_x = NULL, 
                      col_main = "black",
                      col_pal = c("blue", "black", "red"),
                      ...) {
  
  # Visualizes main effects or interactions from regression models.
  #
  # Arguments:
  #          model: A fitted model object (e.g., lm, glm, gam).
  #        formula: A one-sided formula (e.g., ~ x or ~ x + group) or a 
  #                 character vector specifying variable names.
  #   group_levels: Optional vector of specific levels/values for the grouping
  #                 variable. Defaults to min and 95th percentile if NULL.
  #    transform_x: Optional function to transform the x-axis for plotting 
  #                 (e.g., exp, expm1). Defaults to NULL (no transformation).
  #       col_main: Color for the main effect line/interval (used if input 
  #                 has 1 variable).
  #        col_pal: Vector of colors for interaction levels (used if input 
  #                 has 2 variables).
  #            ...: Additional graphical parameters passed to base plot().
  #
  # Returns:
  #    Invisibly returns the prediction data frame containing the grid, 
  #    predictions (fit), confidence intervals (lwr, upr), and plotting x-values.
  #
  # Details:
  #    The function constructs a prediction grid by holding nuisance variables at
  #    their mean (numeric) or reference level (factors). For 'mgcv::gam' objects,
  #    it identifies random effects (bs="re") and excludes them from prediction
  #    by passing their exact term labels to the 'exclude' argument.
  
  # 1. Parse Input (Formula or Character)
  # -------------------------------------
  if (inherits(formula, "formula")) {
    plot_vars <- all.vars(formula)
  } else if (is.character(formula)) {
    plot_vars <- formula
  } else {
    stop("Argument 'formula' must be a formula or a character vector.")
  }
  
  n_vars <- length(plot_vars)
  
  if (n_vars == 1) {
    # Main effect mode
    x_name <- plot_vars[1]
    g_name <- NULL
  } else if (n_vars == 2) {
    # Interaction mode
    x_name <- plot_vars[1]
    g_name <- plot_vars[2]
  } else {
    stop("Must specify exactly 1 or 2 variables.")
  }
  
  # 2. Extract Data and Setup Grid
  # ------------------------------
  dat <- model.frame(model)
  
  req_vars <- if (is.null(g_name)) x_name else c(x_name, g_name)
  if (!all(req_vars %in% names(dat))) {
    stop("Variables specified not found in model frame.")
  }
  
  # Create smooth sequence for x-axis
  x_seq <- seq(min(dat[[x_name]], na.rm = TRUE), 
               max(dat[[x_name]], na.rm = TRUE), 
               length.out = 100)
  
  # Setup grid based on mode
  if (is.null(g_name)) {
    # Main Effect
    pred_grid <- data.frame(x_seq)
    names(pred_grid) <- x_name
    
    # Internal dummy grouping for unified loop
    pred_grid$.group_dummy <- factor("Main Effect")
    plot_g_name <- ".group_dummy"
    
  } else {
    # Interaction
    if (is.null(group_levels)) {
      group_levels <- stats::quantile(dat[[g_name]], probs = c(0, 0.95), na.rm = TRUE)
    }
    pred_grid <- expand.grid(x_seq, group_levels)
    names(pred_grid) <- c(x_name, g_name)
    
    plot_g_name <- g_name
  }
  
  # 3. Handle Nuisance Variables
  # ----------------------------
  model_terms <- terms(model)
  all_model_vars <- all.vars(model_terms)
  
  resp_idx <- attr(model_terms, "response")
  resp_name <- if (resp_idx > 0) all_model_vars[resp_idx] else NULL
  
  vars_to_exclude <- c(resp_name, x_name, "k")
  if (!is.null(g_name)) vars_to_exclude <- c(vars_to_exclude, g_name)
  
  nuisance_vars <- setdiff(all_model_vars, vars_to_exclude)
  
  for (v in nuisance_vars) {
    if (!v %in% names(dat)) next 
    
    if (is.numeric(dat[[v]])) {
      pred_grid[[v]] <- mean(dat[[v]], na.rm = TRUE)
    } else if (is.factor(dat[[v]])) {
      pred_grid[[v]] <- levels(dat[[v]])[1]
    } else {
      pred_grid[[v]] <- stats::na.omit(dat[[v]])[1]
    }
  }
  
  # 4. Predict
  # ----------
  pred_args <- list(object = model, newdata = pred_grid, type = "link", se.fit = TRUE)
  
  # Strict Random Effect Handling for GAMs
  if (inherits(model, "gam") && !is.null(model$smooth)) {
    is_re <- vapply(model$smooth, inherits, logical(1), what = "random.effect")
    
    if (any(is_re)) {
      # Extract exact term labels for exclusion
      re_terms <- vapply(model$smooth[is_re], function(s) s$label, character(1))
      pred_args$exclude <- re_terms
    }
  }
  
  pred_link <- do.call(stats::predict, pred_args)
  
  # 5. Calculate CI and Back-transform
  # ----------------------------------
  fam <- stats::family(model)
  linkinv <- fam$linkinv
  crit_val <- stats::qnorm(0.975)
  
  fit_link <- pred_link$fit
  upr_link <- fit_link + (crit_val * pred_link$se.fit)
  lwr_link <- fit_link - (crit_val * pred_link$se.fit)
  
  pred_grid$fit <- linkinv(fit_link)
  pred_grid$upr <- linkinv(upr_link)
  pred_grid$lwr <- linkinv(lwr_link)
  
  # 6. Transform for Plotting
  # -------------------------
  if (is.function(transform_x)) {
    pred_grid$x_plot <- transform_x(pred_grid[[x_name]])
  } else {
    pred_grid$x_plot <- pred_grid[[x_name]]
  }
  
  # 7. Plotting Logic
  # -----------------
  u_groups <- sort(unique(pred_grid[[plot_g_name]]))
  n_groups <- length(u_groups)
  
  if (is.null(g_name)) {
    plot_cols <- rep(col_main, n_groups)
  } else {
    plot_cols <- rep_len(col_pal, n_groups)
  }
  
  poly_cols <- grDevices::adjustcolor(plot_cols, alpha.f = 0.2)
  ltys <- seq_len(n_groups)
  
  dots <- list(...)
  if (is.null(dots$ylab)) dots$ylab <- if (!is.null(resp_name)) resp_name else "Response"
  if (is.null(dots$xlab)) dots$xlab <- x_name
  if (is.null(dots$las))  dots$las  <- 1
  
  if (is.null(dots$ylim)) {
    dots$ylim <- range(c(pred_grid$lwr, pred_grid$upr), na.rm = TRUE)
  }
  
  base_args <- list(x = pred_grid$x_plot, y = pred_grid$fit, type = "n")
  
  do.call(graphics::plot, c(base_args, dots))
  
  for (i in seq_along(u_groups)) {
    sub_dat <- pred_grid[pred_grid[[plot_g_name]] == u_groups[i], ]
    
    graphics::polygon(x = c(sub_dat$x_plot, rev(sub_dat$x_plot)),
                      y = c(sub_dat$lwr, rev(sub_dat$upr)),
                      col = poly_cols[i], border = NA)
    
    graphics::lines(sub_dat$x_plot, sub_dat$fit, 
                    col = plot_cols[i], lty = ltys[i], lwd = 2)
  }
  
  return(invisible(pred_grid))
}


format_pval <- function(p, threshold = 0.001, digits = 3, sig_level = 0.05, show_stars = FALSE) {
  
  # Formats p-values for plots with prefixes, bolding, and optional stars.
  #
  # Arguments:
  #   p:          Numeric vector of p-values.
  #   threshold:  Numeric cutoff for displaying '<' (default: 0.001).
  #   digits:     Integer for decimal places (default: 3).
  #   sig_level:  Numeric cutoff for bold formatting (default: 0.05).
  #   show_stars: Logical. If TRUE, adds significance stars (* <.05, ** <.01, *** <.001).
  #
  # Returns:
  #   An expression vector suitable for use in 'text()' or 'legend()'.
  #
  # Details:
  #   Uses Base R's `symnum` function to generate standard significance stars
  #   if `show_stars` is TRUE. Returns an expression vector to support
  #   formatting (bold) via plotmath.
  
  # 1. Input validation
  if (!is.numeric(p)) stop("'p' must be a numeric vector.")
  
  # 2. Pre-format numeric strings (e.g., "0.042")
  val_str <- formatC(p, digits = digits, format = "f")
  thresh_str <- formatC(threshold, digits = digits, format = "f")
  
  # 3. Build expression vector
  expr_list <- lapply(seq_along(p), function(i) {
    
    val <- p[i]
    
    # Handle NA explicitly
    if (is.na(val)) return(bquote("NA"))
    
    # --- Step A: Build the text label ---
    # Check if value is below the display threshold (e.g. < 0.001)
    if (val < threshold) {
      lbl <- paste0("p < ", thresh_str)
    } else {
      lbl <- paste0("p = ", val_str[i])
    }
    
    # --- Step B: Add Stars (Optional) ---
    if (show_stars) {
      # symnum is a robust Base R tool for mapping numbers to symbols
      star_txt <- symnum(val, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                         symbols = c("***", "**", "*", ""))
      
      # Convert strictly to character (symnum returns a special class)
      star_txt <- as.character(star_txt)
      
      # Append to label with a space if a star exists
      if (nchar(star_txt) > 0) {
        lbl <- paste(lbl, star_txt)
      }
    }
    
    # --- Step C: Apply Bolding ---
    # If significant, wrap the whole string in bold()
    if (val < sig_level) {
      return(bquote(bold(.(lbl))))
    } else {
      return(bquote(.(lbl)))
    }
  })
  
  # Return as a strictly subsettable expression vector
  as.expression(expr_list)
}
