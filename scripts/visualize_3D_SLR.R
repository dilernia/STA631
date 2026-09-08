#################
# Title: Visualize Optimal Regression Estimates in 3D
# Author: Andrew DiLernia
# Date: 09/08/2026
# Purpose: Visualize optimal simple linear regression estimates in 3D
#################

# This code was generated "collaboratively" with Claude.

library(plotly)

#' Plot an interactive 3D SSE surface for simple linear regression
#'
#' Simulates data from a known linear model, fits OLS (equivalent to the MLE
#' under normal errors), and builds an interactive plotly surface of SSE as a
#' function of the intercept and slope, with the OLS/MLE minimum marked.
#'
#' @param n Sample size for the simulated data.
#' @param true_intercept True intercept used to simulate y.
#' @param true_slope True slope used to simulate y.
#' @param x_min,x_max Range from which x is drawn (uniformly).
#' @param error_sd Standard deviation of the simulated normal error term.
#' @param grid_n Number of grid points along each axis (intercept, slope).
#' @param intercept_pad,slope_pad Half-widths of the grid around the OLS
#'   estimate, in the intercept and slope directions respectively.
#' @param marker_offset_pct Fraction of the SSE range used to nudge the
#'   minimum marker upward, so its bottom vertex approximates sitting on the
#'   surface rather than being centered on it.
#' @param seed Random seed for reproducibility.
#'
#' @return A plotly htmlwidget.
plot_sse_surface <- function(n = 50,
                             true_intercept = 2,
                             true_slope = 1.5,
                             x_min = 3,
                             x_max = 7,
                             error_sd = 2,
                             grid_n = 80,
                             intercept_pad = 10,
                             slope_pad = 5,
                             marker_offset_pct = 0.03,
                             seed = 123) {
  
  # 1. Simulate data for a simple linear regression -----------------------
  set.seed(seed)
  x <- runif(n, x_min, x_max)  # narrower spread gives the slope less
  # leverage over SSE, so the intercept
  # matters relatively more
  y <- true_intercept + true_slope * x + rnorm(n, sd = error_sd)
  
  # 2. Fit OLS (which is also the MLE under normal errors) for reference --
  fit <- lm(y ~ x)
  ols_intercept <- coef(fit)[1]
  ols_slope     <- coef(fit)[2]
  
  # 3. Define the SSE (sum of squared errors) function ---------------------
  #    SSE(b0, b1) = sum((y_i - (b0 + b1 * x_i))^2)
  sse <- function(b0, b1) {
    sum((y - (b0 + b1 * x))^2)
  }
  
  # 4. Build a grid of (intercept, slope) values around the OLS solution --
  intercept_seq <- seq(ols_intercept - intercept_pad, ols_intercept + intercept_pad, length.out = grid_n)
  slope_seq     <- seq(ols_slope - slope_pad, ols_slope + slope_pad, length.out = grid_n)
  
  # outer() gives a matrix with rows indexed by intercept_seq, cols by slope_seq
  sse_matrix <- outer(intercept_seq, slope_seq, Vectorize(sse))
  
  # plotly wants z[i, j] to correspond to x[j], y[i], so transpose:
  # rows -> slope_seq (y-axis), cols -> intercept_seq (x-axis)
  z_matrix <- t(sse_matrix)
  
  # 5. Build the interactive surface ---------------------------------------
  fig <- plot_ly(
    x = ~intercept_seq,
    y = ~slope_seq,
    z = ~z_matrix,
    type = "surface",
    colorscale = "Viridis",
    colorbar = list(title = "SSE"),
    contours = list(
      z = list(show = TRUE, usecolormap = TRUE, project = list(z = TRUE))
    ),
    hovertemplate = paste(
      "Intercept: %{x:.3f}<br>",
      "Slope: %{y:.3f}<br>",
      "SSE: %{z:.3f}<extra></extra>"
    )
  ) %>%
    layout(
      title = "SSE Surface for Simple Linear Regression",
      scene = list(
        xaxis = list(title = "Intercept (\u03B2\u2080)"),
        yaxis = list(title = "Slope (\u03B2\u2081)"),
        zaxis = list(title = "SSE"),
        camera = list(eye = list(x = 1.5, y = -1.5, z = 0.8))
      )
    )
  
  # 6. Mark the minimum (the OLS / MLE solution) on the surface -----------
  # plotly 3D markers are always centered on their (x, y, z) coordinate —
  # there's no built-in way to anchor a different part of the marker to a
  # point. To make the *bottom* of the diamond appear to sit on the surface,
  # we nudge its z value up by a small offset (a fraction of the overall SSE
  # range). This is an approximation: the exact visual alignment will drift
  # slightly as you zoom, since marker size is fixed in screen pixels, not
  # data units.
  z_range <- diff(range(z_matrix))
  marker_z_offset <- marker_offset_pct * z_range
  min_sse <- sse(ols_intercept, ols_slope)
  
  fig <- fig %>%
    add_trace(
      x = ols_intercept,
      y = ols_slope,
      z = min_sse + marker_z_offset,
      type = "scatter3d",
      mode = "markers",
      marker = list(color = "orange", size = 8, symbol = "diamond"),
      name = "OLS / MLE minimum",
      showlegend = TRUE,
      hovertemplate = paste(
        "Intercept: %{x:.3f}<br>",
        "Slope: %{y:.3f}<br>",
        "SSE: ", sprintf("%.3f", min_sse), "<extra>OLS / MLE minimum</extra>"
      )
    )
  
  fig
}