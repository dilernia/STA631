#################
# Title: Visualize Optimal Regression Estimates in 3D
# Author: Andrew DiLernia
# Date: 09/08/2026
# Purpose: Visualize optimal simple linear regression estimates in 3D
#          (SSE surface or log-likelihood surface)
#################

#' Plot an interactive 3D SSE or log-likelihood surface for simple linear
#' regression
#'
#' Takes a provided dataset, fits OLS (equivalent to the MLE under normal errors), 
#' and builds an interactive plotly surface of either SSE or the (profile) 
#' log-likelihood as a function of the intercept and slope, with the exact OLS/MLE 
#' optimum marked for the provided data.
#'
#' @param data A data frame or tibble containing 'x' and 'y' columns.
#' @param grid_n Number of grid points along each axis (intercept, slope).
#' @param intercept_pad,slope_pad Half-widths of the grid around the OLS
#'   estimate, in the intercept and slope directions respectively.
#' @param marker_offset_pct Fraction of the surface's z-range used to nudge
#'   the optimum marker away from the surface (upward), so its bottom vertex
#'   approximates sitting on the surface rather than being centered on it.
#' @param surface Which surface to plot: "sse" or "likelihood".
#'
#' @return A plotly htmlwidget.
plot_regression_surface <- function(data = NULL,
                                    grid_n = 80,
                                    intercept_pad = 10,
                                    slope_pad = 5,
                                    marker_offset_pct = 0.03,
                                    surface = c("sse", "likelihood")) {
  
  if(is.null(data)) {
    # 1. Simulate data from SLR model -------------------------------------------
    set.seed(1994)
    
    # True parameters
    true_intercept_slope <- c(4, 0.6)
    n_values <- 15
    
    data <- tibble::tibble(
      x = seq(1, 10, length.out = n_values),
      y = true_intercept_slope[1] + true_intercept_slope[2] * x + stats::rnorm(n_values, mean = 0, sd = 2)
    )
  }
  
  surface <- match.arg(surface)
  
  # Extract variables for modeling directly from the provided dataset
  x <- data$x
  y <- data$y
  n_obs <- length(x)
  
  # 2. Fit OLS (which is also the MLE under normal errors) for reference --
  # This ensures the optimum matches the estimator values for the provided data
  fit <- stats::lm(y ~ x)
  ols_intercept <- stats::coef(fit)[1]
  ols_slope     <- stats::coef(fit)[2]
  
  # 3. Define the SSE function and, from it, the profile log-likelihood ---
  sse <- function(b0, b1) {
    sum((y - (b0 + b1 * x))^2)
  }
  
  loglik <- function(b0, b1) {
    s <- sse(b0, b1)
    -n_obs / 2 * log(2 * pi) - n_obs / 2 * log(s / n_obs) - n_obs / 2
  }
  
  z_fun <- if (surface == "sse") sse else loglik
  
  # 4. Build a grid of (intercept, slope) values around the OLS solution --
  intercept_seq <- seq(ols_intercept - intercept_pad, ols_intercept + intercept_pad, length.out = grid_n)
  slope_seq     <- seq(ols_slope - slope_pad, ols_slope + slope_pad, length.out = grid_n)
  
  z_grid_matrix <- outer(intercept_seq, slope_seq, Vectorize(z_fun))
  z_matrix <- t(z_grid_matrix)
  
  # 5. Build the interactive surface ---------------------------------------
  z_axis_title <- if (surface == "sse") "SSE" else "Log-Likelihood"
  plot_title <- if (surface == "sse") {
    "SSE Surface for Simple Linear Regression"
  } else {
    "Log-Likelihood Surface for Simple Linear Regression"
  }
  
  fig <- plotly::plot_ly(
    x = ~intercept_seq,
    y = ~slope_seq,
    z = ~z_matrix,
    type = "surface",
    colorscale = "Viridis",
    colorbar = list(title = z_axis_title),
    contours = list(
      z = list(show = TRUE, usecolormap = TRUE, project = list(z = TRUE))
    ),
    hovertemplate = paste0(
      "Intercept: %{x:.3f}<br>",
      "Slope: %{y:.3f}<br>",
      z_axis_title, ": %{z:.3f}<extra></extra>"
    )
  )
  
  fig <- plotly::layout(
    fig,
    title = plot_title,
    scene = list(
      xaxis = list(title = "Intercept (\u03B2\u2080)"),
      yaxis = list(title = "Slope (\u03B2\u2081)"),
      zaxis = list(title = z_axis_title),
      camera = list(eye = list(x = 1.5, y = -1.5, z = 0.8))
    )
  )
  
  # 6. Mark the optimum (the OLS / MLE solution) on the surface -----------
  z_range <- diff(range(z_matrix))
  marker_z_offset <- marker_offset_pct * z_range
  opt_z <- z_fun(ols_intercept, ols_slope)
  marker_label <- if (surface == "sse") "SSE minimum" else "Likelihood maximum"
  
  fig <- plotly::add_trace(
    fig,
    x = ols_intercept,
    y = ols_slope,
    z = opt_z + marker_z_offset,
    type = "scatter3d",
    mode = "markers",
    marker = list(color = "orange", size = 8, symbol = "diamond"),
    name = marker_label,
    showlegend = TRUE,
    hovertemplate = paste0(
      "Intercept: %{x:.3f}<br>",
      "Slope: %{y:.3f}<br>",
      z_axis_title, ": ", sprintf("%.3f", opt_z), "<extra>", marker_label, "</extra>"
    )
  )
  
  fig
}

# Examples ------------------------------------------------------------------
# plot_regression_surface(surface = "sse")
# plot_regression_surface(surface = "likelihood")