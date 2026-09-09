#################
# Title: Visualize Optimal Regression Estimates in 3D
# Author: Andrew DiLernia
# Date: 09/08/2026
# Purpose: Visualize optimal simple linear regression estimates in 3D
#          (SSE surface or log-likelihood surface)
#################

# This code was generated "collaboratively" with Claude.

library(plotly)

#' Plot an interactive 3D SSE or log-likelihood surface for simple linear
#' regression
#'
#' Simulates data from a known linear model, fits OLS (equivalent to the MLE
#' under normal errors), and builds an interactive plotly surface of either
#' SSE or the (profile) log-likelihood as a function of the intercept and
#' slope, with the OLS/MLE optimum marked.
#'
#' @param n Sample size for the simulated data.
#' @param true_intercept True intercept used to simulate y.
#' @param true_slope True slope used to simulate y.
#' @param x_min,x_max Range from which x is drawn (uniformly).
#' @param error_sd Standard deviation of the simulated normal error term.
#' @param grid_n Number of grid points along each axis (intercept, slope).
#' @param intercept_pad,slope_pad Half-widths of the grid around the OLS
#'   estimate, in the intercept and slope directions respectively.
#' @param marker_offset_pct Fraction of the surface's z-range used to nudge
#'   the optimum marker away from the surface (upward), so its bottom vertex
#'   approximates sitting on the surface rather than being centered on it.
#' @param surface Which surface to plot: "sse" (sum of squared errors, a
#'   bowl with a minimum at the OLS estimate) or "likelihood" (the profile
#'   log-likelihood, concentrating out sigma at each (b0, b1) grid point, a
#'   dome with a maximum at the MLE). Under normality, SSE and the profile
#'   log-likelihood are monotonic transforms of one another, so both surfaces
#'   are optimized at the same (intercept, slope) point -- only the shape and
#'   the direction of the optimum (min vs. max) differ.
#' @param seed Random seed for reproducibility.
#'
#' @return A plotly htmlwidget.
plot_regression_surface <- function(n = 50,
                                    true_intercept = 2,
                                    true_slope = 1.5,
                                    x_min = 3,
                                    x_max = 7,
                                    error_sd = 2,
                                    grid_n = 80,
                                    intercept_pad = 10,
                                    slope_pad = 5,
                                    marker_offset_pct = 0.03,
                                    surface = c("sse", "likelihood"),
                                    seed = 123) {
  
  surface <- match.arg(surface)
  
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
  
  # 3. Define the SSE function and, from it, the profile log-likelihood ---
  #    SSE(b0, b1) = sum((y_i - (b0 + b1 * x_i))^2)
  sse <- function(b0, b1) {
    sum((y - (b0 + b1 * x))^2)
  }
  
  #    Concentrated (profile) normal log-likelihood: at each (b0, b1), sigma
  #    is set to its own MLE, sigma_hat = sqrt(SSE(b0, b1) / n), and plugged
  #    back in. This collapses the 3-parameter (b0, b1, sigma) likelihood
  #    down to a function of (b0, b1) alone, matching the SSE surface's
  #    axes.
  #      ll(b0, b1) = -n/2 * log(2*pi) - n/2 * log(SSE(b0, b1) / n) - n/2
  loglik <- function(b0, b1) {
    s <- sse(b0, b1)
    -n / 2 * log(2 * pi) - n / 2 * log(s / n) - n / 2
  }
  
  z_fun <- if (surface == "sse") sse else loglik
  
  # 4. Build a grid of (intercept, slope) values around the OLS solution --
  intercept_seq <- seq(ols_intercept - intercept_pad, ols_intercept + intercept_pad, length.out = grid_n)
  slope_seq     <- seq(ols_slope - slope_pad, ols_slope + slope_pad, length.out = grid_n)
  
  # outer() gives a matrix with rows indexed by intercept_seq, cols by slope_seq
  z_grid_matrix <- outer(intercept_seq, slope_seq, Vectorize(z_fun))
  
  # plotly wants z[i, j] to correspond to x[j], y[i], so transpose:
  # rows -> slope_seq (y-axis), cols -> intercept_seq (x-axis)
  z_matrix <- t(z_grid_matrix)
  
  # 5. Build the interactive surface ---------------------------------------
  z_axis_title <- if (surface == "sse") "SSE" else "Log-Likelihood"
  plot_title <- if (surface == "sse") {
    "SSE Surface for Simple Linear Regression"
  } else {
    "Log-Likelihood Surface for Simple Linear Regression"
  }
  
  fig <- plot_ly(
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
  ) %>%
    layout(
      title = plot_title,
      scene = list(
        xaxis = list(title = "Intercept (\u03B2\u2080)"),
        yaxis = list(title = "Slope (\u03B2\u2081)"),
        zaxis = list(title = z_axis_title),
        camera = list(eye = list(x = 1.5, y = -1.5, z = 0.8))
      )
    )
  
  # 6. Mark the optimum (the OLS / MLE solution) on the surface -----------
  # plotly 3D markers are always centered on their (x, y, z) coordinate —
  # there's no built-in way to anchor a different part of the marker to a
  # point. To make the *bottom* of the diamond appear to sit on the surface,
  # we nudge its z value up by a small offset (a fraction of the overall
  # z-range). This is an approximation: the exact visual alignment will
  # drift slightly as you zoom, since marker size is fixed in screen pixels,
  # not data units. This nudge is applied the same way for both surfaces:
  # SSE's optimum is a minimum (a bowl) and the log-likelihood's optimum is
  # a maximum (a dome), but in both cases the marker sits at the highest
  # local point on the surface, so an upward nudge keeps its bottom vertex
  # anchored there.
  z_range <- diff(range(z_matrix))
  marker_z_offset <- marker_offset_pct * z_range
  opt_z <- z_fun(ols_intercept, ols_slope)
  marker_label <- if (surface == "sse") "SSE minimum" else "Likelihood maximum"
  
  fig <- fig %>%
    add_trace(
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