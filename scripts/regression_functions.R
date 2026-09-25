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
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a standard linear model
#' plot_regression_surface(surface = "sse")
#' plot_regression_surface(surface = "likelihood")
#' }
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

#' Create an Interactive 3D Regression Surface Plot
#'
#' @description
#' Generates an interactive 3D scatter plot and a fitted regression plane using `plotly`.
#' The function dynamically extracts the data, variable names, and estimated coefficients 
#' directly from a fitted two-predictor linear model. The resulting plot title automatically displays 
#' the estimated regression equation, while axis labels and hover text are dynamically 
#' populated using the model's variable names.
#'
#' @param model An object of class \code{lm} representing a multiple linear regression
#'   model. The model must contain exactly one continuous response variable and two 
#'   numeric predictor variables.
#'
#' @return A \code{plotly} object containing the interactive 3D visualization.
#'
#' @export
#'
#' @import plotly
#'
#' @examples
#' # Fit a model predicting MPG based on Weight and Horsepower using the mtcars dataset
#' my_model <- lm(mpg ~ wt + hp, data = mtcars)
#'
#' # Generate the interactive Plotly 3D visualization
#' plot_regression_plane(my_model)
plot_regression_plane <- function(model) {
  
  # 1. Extract variable names dynamically from the model formula
  vars <- all.vars(formula(model))
  if (length(vars) != 3) {
    stop("This function requires a model with exactly one response and two predictors.")
  }
  
  y_name  <- vars[1]
  x1_name <- vars[2]
  x2_name <- vars[3]
  
  # 2. Extract the exact data dataframe used to fit the model
  df <- model$model
  
  # 3. Extract coefficients and construct the dynamic equation title
  b <- coef(model)
  # Uses %.2f for the intercept and %+.2f to force +/- signs for the slopes
  eq_string <- sprintf("Predicted %s = %.2f %+.2f*%s %+.2f*%s", 
                       y_name, b[1], b[2], x1_name, b[3], x2_name)
  
  # 4. Create a structured grid of values for the regression plane
  axis_x1 <- seq(min(df[[x1_name]]), max(df[[x1_name]]), length.out = 30)
  axis_x2 <- seq(min(df[[x2_name]]), max(df[[x2_name]]), length.out = 30)
  
  plane_grid <- expand.grid(x1 = axis_x1, x2 = axis_x2)
  # Rename columns to perfectly match what predict() expects
  colnames(plane_grid) <- c(x1_name, x2_name)
  
  # Predict values and transpose matrix to align with Plotly's X/Y mapping
  plane_grid$y_hat <- predict(model, newdata = plane_grid)
  plane_z <- t(matrix(plane_grid$y_hat, nrow = length(axis_x1), ncol = length(axis_x2)))
  
  # 5. Build dynamic hover templates using the extracted variable names
  hover_scatter <- sprintf("%s: %%{x:.2f}<br>%s: %%{y:.2f}<br>%s: %%{z:.2f}<extra></extra>", 
                           x1_name, x2_name, y_name)
  hover_plane <- sprintf("%s: %%{x:.2f}<br>%s: %%{y:.2f}<br>Predicted %s: %%{z:.2f}<extra></extra>", 
                         x1_name, x2_name, y_name)
  
  # 6. Generate the Plotly 3D Visualization
  p <- plot_ly() %>%
    
    # Add the 3D scatter points
    add_markers(
      x = df[[x1_name]], 
      y = df[[x2_name]], 
      z = df[[y_name]], 
      marker = list(size = 4, color = "blue", opacity = 0.6),
      name = "Actual Data",
      hovertemplate = hover_scatter
    ) %>%
    
    # Add the regression surface
    add_surface(
      x = axis_x1, 
      y = axis_x2, 
      z = plane_z, 
      opacity = 0.5,
      colorscale = "Viridis",
      showscale = FALSE,
      name = "Regression Plane",
      hovertemplate = hover_plane
    ) %>%
    
    # Customize layout with dynamic labels and forced cube aspect ratio
    layout(
      title = list(text = eq_string, font = list(size = 15)),
      scene = list(
        xaxis = list(title = x1_name),
        yaxis = list(title = x2_name),
        zaxis = list(title = y_name),
        camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2)),
        aspectmode = "cube"
      )
    )
  
  return(p)
}

#' Plot Variance Inflation Factor (VIF) Values
#'
#' @description
#' Calculates and plots the Variance Inflation Factor (VIF) or Generalized VIF (GVIF) 
#' for a fitted regression model. It generates a lollipop chart of the VIF values ranked 
#' in descending order, featuring reference lines at 5 and 10 to help visually identify 
#' potential multicollinearity.
#'
#' @param modFit A fitted model object (e.g., \code{lm}, \code{glm}) supported by \code{car::vif()}.
#'
#' @return A \code{ggplot} object displaying the calculated VIF or GVIF values.
#'
#' @importFrom car vif
#' @importFrom dplyr mutate
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 ggplot aes geom_segment geom_vline scale_x_continuous expansion geom_point labs
#' @importFrom magrittr %>%
#' @importFrom stringr str_wrap
#' @importFrom tibble tibble
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a standard linear model
#' model <- lm(mpg ~ cyl + disp + hp + wt, data = mtcars)
#' 
#' # Generate and print the VIF plot
#' vif_plot(model)
#' }
vif_plot <- function(modFit) {
  
  vifs <- car::vif(modFit)
  
  if("GVIF^(1/(2*Df))" %in% colnames(vifs)) {
    vifs <- vifs[, "GVIF^(1/(2*Df))"]
    ggTitle <- "Generalized variance inflation factors"
    xLab <- "GVIF"
  } else {
    ggTitle <- "Variance inflation factors"
    xLab <- "VIF"
  }
  
  preds <- names(vifs)
  
  vifGG <- tibble(Predictor = preds,
                  VIF = vifs) %>% 
    dplyr::mutate(Predictor = fct_reorder(Predictor, -VIF)) %>% 
    ggplot(aes(x = VIF, y = Predictor)) + 
    geom_segment(aes(xend = VIF, x = 0, 
                     yend = Predictor, y = Predictor)) +
    geom_vline(xintercept = 5, linetype = "dotted") +
    geom_vline(xintercept = 10, linetype = "dotted") +
    scale_x_continuous(limits = c(0, max(vifs)*1.1),
                       expand = expansion(mult = c(0, 0.10))) +
    geom_point(size = 3, color = "steelblue") +
    labs(title = ggTitle, x = xLab,
         caption = stringr::str_wrap(paste0(xLab, " values calculated via the car package: https://search.r-project.org/CRAN/refmans/car/html/vif.html")))
  
  return(suppressWarnings(print(vifGG)))
}
