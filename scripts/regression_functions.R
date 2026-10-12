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
                                    surface = base::c("sse", "likelihood")) {
  
  if(base::is.null(data)) {
    # 1. Simulate data from SLR model -------------------------------------------
    base::set.seed(1994)
    
    # True parameters
    true_intercept_slope <- base::c(4, 0.6)
    n_values <- 15
    
    data <- tibble::tibble(
      x = base::seq(1, 10, length.out = n_values),
      y = true_intercept_slope[1] + true_intercept_slope[2] * x + stats::rnorm(n_values, mean = 0, sd = 2)
    )
  }
  
  surface <- base::match.arg(surface)
  
  # Extract variables for modeling directly from the provided dataset
  x <- data$x
  y <- data$y
  n_obs <- base::length(x)
  
  # 2. Fit OLS (which is also the MLE under normal errors) for reference --
  # This ensures the optimum matches the estimator values for the provided data
  fit <- stats::lm(y ~ x)
  ols_intercept <- stats::coef(fit)[1]
  ols_slope     <- stats::coef(fit)[2]
  
  # 3. Define the SSE function and, from it, the profile log-likelihood ---
  sse <- function(b0, b1) {
    base::sum((y - (b0 + b1 * x))^2)
  }
  
  loglik <- function(b0, b1) {
    s <- sse(b0, b1)
    -n_obs / 2 * base::log(2 * base::pi) - n_obs / 2 * base::log(s / n_obs) - n_obs / 2
  }
  
  z_fun <- if (surface == "sse") sse else loglik
  
  # 4. Build a grid of (intercept, slope) values around the OLS solution --
  intercept_seq <- base::seq(ols_intercept - intercept_pad, ols_intercept + intercept_pad, length.out = grid_n)
  slope_seq     <- base::seq(ols_slope - slope_pad, ols_slope + slope_pad, length.out = grid_n)
  
  z_grid_matrix <- base::outer(intercept_seq, slope_seq, base::Vectorize(z_fun))
  z_matrix <- base::t(z_grid_matrix)
  
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
    colorbar = base::list(title = z_axis_title),
    contours = base::list(
      z = base::list(show = TRUE, usecolormap = TRUE, project = base::list(z = TRUE))
    ),
    hovertemplate = base::paste0(
      "Intercept: %{x:.3f}<br>",
      "Slope: %{y:.3f}<br>",
      z_axis_title, ": %{z:.3f}<extra></extra>"
    )
  )
  
  fig <- plotly::layout(
    fig,
    title = plot_title,
    scene = base::list(
      xaxis = base::list(title = "Intercept (\u03B2\u2080)"),
      yaxis = base::list(title = "Slope (\u03B2\u2081)"),
      zaxis = base::list(title = z_axis_title),
      camera = base::list(eye = base::list(x = 1.5, y = -1.5, z = 0.8))
    )
  )
  
  # 6. Mark the optimum (the OLS / MLE solution) on the surface -----------
  z_range <- base::diff(base::range(z_matrix))
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
    marker = base::list(color = "orange", size = 8, symbol = "diamond"),
    name = marker_label,
    showlegend = TRUE,
    hovertemplate = base::paste0(
      "Intercept: %{x:.3f}<br>",
      "Slope: %{y:.3f}<br>",
      z_axis_title, ": ", base::sprintf("%.3f", opt_z), "<extra>", marker_label, "</extra>"
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
#' my_model <- stats::lm(mpg ~ wt + hp, data = mtcars)
#'
#' # Generate the interactive Plotly 3D visualization
#' plot_regression_plane(my_model)
plot_regression_plane <- function(model) {
  
  # 1. Extract variable names dynamically from the model formula
  vars <- base::all.vars(stats::formula(model))
  if (base::length(vars) != 3) {
    base::stop("This function requires a model with exactly one response and two predictors.")
  }
  
  y_name  <- vars[1]
  x1_name <- vars[2]
  x2_name <- vars[3]
  
  # 2. Extract the exact data dataframe used to fit the model
  df <- model$model
  
  # 3. Extract coefficients and construct the dynamic equation title
  b <- stats::coef(model)
  
  # Uses %.2f for the intercept and %+.2f to force +/- signs for the slopes
  eq_string <- base::sprintf("Predicted %s = %.2f %+.2f*%s %+.2f*%s", 
                             y_name, b[1], b[2], x1_name, b[3], x2_name)
  
  # 4. Create a structured grid of values for the regression plane
  axis_x1 <- base::seq(base::min(df[[x1_name]]), base::max(df[[x1_name]]), length.out = 30)
  axis_x2 <- base::seq(base::min(df[[x2_name]]), base::max(df[[x2_name]]), length.out = 30)
  
  plane_grid <- base::expand.grid(x1 = axis_x1, x2 = axis_x2)
  # Rename columns to perfectly match what predict() expects
  base::colnames(plane_grid) <- base::c(x1_name, x2_name)
  
  # Predict values and transpose matrix to align with Plotly's X/Y mapping
  plane_grid$y_hat <- stats::predict(model, newdata = plane_grid)
  plane_z <- base::t(base::matrix(plane_grid$y_hat, nrow = base::length(axis_x1), ncol = base::length(axis_x2)))
  
  # 5. Build dynamic hover templates using the extracted variable names
  hover_scatter <- base::sprintf("%s: %%{x:.2f}<br>%s: %%{y:.2f}<br>%s: %%{z:.2f}<extra></extra>", 
                                 x1_name, x2_name, y_name)
  hover_plane <- base::sprintf("%s: %%{x:.2f}<br>%s: %%{y:.2f}<br>Predicted %s: %%{z:.2f}<extra></extra>", 
                               x1_name, x2_name, y_name)
  
  # 6. Generate the Plotly 3D Visualization
  p <- plotly::plot_ly() %>%
    
    # Add the 3D scatter points
    plotly::add_markers(
      x = df[[x1_name]], 
      y = df[[x2_name]], 
      z = df[[y_name]], 
      marker = base::list(size = 4, color = "blue", opacity = 0.6),
      name = "Actual Data",
      hovertemplate = hover_scatter
    ) %>%
    
    # Add the regression surface
    plotly::add_surface(
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
    plotly::layout(
      title = base::list(text = eq_string, font = base::list(size = 15)),
      scene = base::list(
        xaxis = base::list(title = x1_name),
        yaxis = base::list(title = x2_name),
        zaxis = base::list(title = y_name),
        camera = base::list(eye = base::list(x = 1.5, y = 1.5, z = 1.2)),
        aspectmode = "cube"
      )
    )
  
  base::return(p)
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
#' model <- stats::lm(mpg ~ cyl + disp + hp + wt, data = mtcars)
#' 
#' # Generate and print the VIF plot
#' vif_plot(model)
#' }
vif_plot <- function(modFit) {
  
  vifs <- car::vif(modFit)
  
  if("GVIF^(1/(2*Df))" %in% base::colnames(vifs)) {
    vifs <- vifs[, "GVIF^(1/(2*Df))"]
    ggTitle <- "Generalized variance inflation factors"
    xLab <- "GVIF"
  } else {
    ggTitle <- "Variance inflation factors"
    xLab <- "VIF"
  }
  
  preds <- base::names(vifs)
  
  vifGG <- tibble::tibble(Predictor = preds,
                          VIF = vifs) %>% 
    dplyr::mutate(Predictor = forcats::fct_reorder(Predictor, -VIF)) %>% 
    ggplot2::ggplot(ggplot2::aes(x = VIF, y = Predictor)) + 
    ggplot2::geom_segment(ggplot2::aes(xend = VIF, x = 0, 
                                       yend = Predictor, y = Predictor)) +
    ggplot2::geom_vline(xintercept = 5, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 10, linetype = "dotted") +
    ggplot2::scale_x_continuous(limits = base::c(0, base::max(vifs)*1.1),
                                expand = ggplot2::expansion(mult = base::c(0, 0.10))) +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::labs(title = ggTitle, x = xLab,
                  caption = stringr::str_wrap(base::paste0(xLab, " values calculated via the car package: https://search.r-project.org/CRAN/refmans/car/html/vif.html")))
  
  base::return(base::suppressWarnings(base::print(vifGG)))
}

#' Create a 3D Contour or Surface Plot
#'
#' This function takes a three-column data frame, interpolates the data over a grid
#' using the `akima` package, and generates either a static filled contour plot using
#' `ggplot2` or an interactive 3D surface plot using `plotly`. 
#'
#' @details
#' Interpolating to plot contours is based on techniques discussed at:
#' \url{https://stackoverflow.com/questions/65873211/empty-contour-plot-in-ggplot}
#'
#' @param data_3d A data frame containing at least three columns representing x, y, and z coordinates.
#'   Missing values are automatically removed.
#' @param dupes A character string indicating how to handle duplicate coordinates.
#'   Passed to the `duplicate` argument of \code{\link[akima]{interp}}. Defaults to \code{"mean"}.
#' @param plot_type A character string specifying the output format. Must be either 
#'   \code{"static"} (returns a \code{ggplot2} object) or \code{"interactive"} (returns a 
#'   \code{plotly} object). Defaults to \code{"static"}.
#'
#' @return A \code{ggplot} object or a \code{plotly} widget depending on the \code{plot_type} argument.
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr rename pull
#' @importFrom akima interp
#' @importFrom ggplot2 ggplot aes geom_contour_filled geom_point labs theme_bw theme element_text element_blank
#' @importFrom plotly plot_ly layout colorbar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(x = rnorm(100), y = rnorm(100), z = runif(100))
#' 
#' # Generate a static ggplot2 contour plot
#' plot_3d(df, plot_type = "static")
#' 
#' # Generate an interactive plotly 3D surface
#' plot_3d(df, plot_type = "interactive")
#' }
plot_3d <- function(data_3d, dupes = "mean", plot_type = c("static", "interactive")) {
  
  plot_type <- match.arg(plot_type)
  
  data_3d_clean <- data_3d |> tidyr::drop_na() |> 
    dplyr::rename("x" = 1, "y" = 2, "z" = 3)
  
  suppressWarnings(grid <- akima::interp(dplyr::pull(data_3d_clean, 1), 
                                         dplyr::pull(data_3d_clean, 2),
                                         dplyr::pull(data_3d_clean, 3),
                                         duplicate = dupes))
  
  my_vars <- colnames(data_3d)
  
  if (plot_type == "static") {
    
    grid_df <- data.frame(x = rep(grid$x, ncol(grid$z)), 
                          y = rep(grid$y, each = nrow(grid$z)), 
                          z = as.numeric(grid$z))
    
    p <- grid_df |>
      ggplot2::ggplot(ggplot2::aes(x = x, y = y, z = z)) +
      ggplot2::geom_contour_filled(ggplot2::aes(x = x, 
                                                y = y, 
                                                z = z)) + 
      ggplot2::geom_point(data = data_3d_clean, colour="white",pch=21, 
                          fill = "black", size = 1.5) +
      ggplot2::labs(x = my_vars[1], y = my_vars[2],
                    fill = my_vars[3],
                    title = "Three-dimensional distribution") +
      ggplot2::theme_bw() + 
      ggplot2::theme(text = ggplot2::element_text(face = "bold"), 
                     panel.grid = ggplot2::element_blank(),
                     legend.position = "bottom")
    
    return(p)
    
  } else if (plot_type == "interactive") {
    
    # Construct a dynamic hover template using the original variable names
    ht <- paste0(
      my_vars[1], ": %{x}<br>",
      my_vars[2], ": %{y}<br>",
      my_vars[3], ": %{z}<extra></extra>" 
    )
    
    p <- plotly::plot_ly(
      x = ~grid$x, 
      y = ~grid$y, 
      z = ~grid$z, 
      type = "surface",
      hovertemplate = ht
    ) |>
      plotly::colorbar(title = my_vars[3]) |>
      plotly::layout(
        title = "Three-dimensional distribution",
        scene = list(
          xaxis = list(title = my_vars[1]),
          yaxis = list(title = my_vars[2]),
          zaxis = list(title = my_vars[3])
        )
      )
    
    return(p)
  }
}

#' Create Empirical Logit Plots
#'
#' This function generates empirical logit plots for quantitative and 
#' categorical predictors from a fitted logistic regression model. 
#' Quantitative predictors are binned, and empirical log odds are plotted 
#' as points. Categorical predictors are plotted as bar charts, with the 
#' bars sorted by their height (the log odds).
#'
#' @param model_fit A fitted logistic regression model object (e.g., from \code{glm()}).
#' @param n_bins An integer specifying the target number of bins for quantitative predictors. Default is 4.
#'
#' @return A \code{ggplot} object if there is a single predictor, or a named list of \code{ggplot} objects if there are multiple predictors.
#' 
#' @import tidyverse
#' @export
empirical_logit_plot <- function(model_fit,
                                 n_bins = 4) {
  
  library(tidyverse)
  
  my_data <- model_fit$model
  response <- colnames(my_data)[1]
  
  # Reference category for the response variable
  ref_cat <- model_fit$data[[response]][as.logical(model_fit$y)][1]
  
  # Making response a dummy variable
  my_data[[response]] <- model_fit[["y"]]
  
  quant_predictors <- colnames(my_data)[-c(1)][which(sapply(my_data, FUN = is.numeric)[-c(1)])]
  cat_predictors <- colnames(my_data)[-c(1)][which(sapply(my_data, FUN = function(x){!is.numeric(x)})[-c(1)])]
  
  # Check for complete separation
  if(length(unique(my_data[[response]])) == 1) {
    stop("Complete separation in data. Cannot create empirical logit plot.")
  }
  
  # Instantiating list for empirical logit plots
  emp_plots <- list()
  
  # Creating empirical logit plots for quantitative predictors
  for(predictor in quant_predictors) {
    
    # Create a dataframe with binned predictor values and the proportion of positive responses in each bin
    agg_data <- my_data |>
      dplyr::group_by(bin = ggplot2::cut_interval(.data[[predictor]], n = n_bins)) |>
      dplyr::summarize(prop_response = mean(.data[[response]])) |>
      dplyr::ungroup() |>
      dplyr::mutate(log_odds = log(prop_response / (1 - prop_response)))
    
    nans <- any(is.infinite(agg_data$log_odds))
    n_bins_new <- n_bins
    
    # Display warning if NaN values
    if(nans) {
      while(nans & n_bins_new >= 1) {
        # Create a dataframe with binned predictor values and the proportion of positive responses in each bin
        agg_data <- my_data |>
          dplyr::group_by(bin = cut_interval(.data[[predictor]], n = n_bins_new - 1)) |>
          dplyr::summarize(prop_response = mean(.data[[response]])) |>
          dplyr::ungroup() |>
          dplyr::mutate(log_odds = log(prop_response / (1 - prop_response)))
        
        # Check for NaN values
        nans <- any(is.infinite(agg_data$log_odds))
        
        # Reduce number of bins by 1
        n_bins_new <- n_bins_new - 1
      }
      warning(paste0("Complete separation for ", predictor, 
                     " based on ",
                     n_bins, " bins."))
    }
    
    used_bins <- nrow(agg_data)
    
    if(n_bins > used_bins) {
      message(paste0("Number of bins reduced to ", 
                     used_bins, " instead of ", n_bins, " for ",
                     predictor, "."))
    }
    
    # Create the empirical logit plot
    emp_plots[[predictor]] <- ggplot2::ggplot(agg_data, ggplot2::aes(x = bin, y = log_odds)) +
      ggplot2::geom_point() +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
      ggplot2::labs(x = predictor,
                    y = paste0("Log odds of ", response, " = ", ref_cat),
                    title = "Empirical Logit Plot") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  
  # Creating empirical logit plots for categorical predictors
  for(predictor in cat_predictors) {
    agg_data <- my_data |>
      dplyr::group_by(!!as.name(predictor)) |>
      dplyr::summarize(prop_response  = mean(!!as.name(response)),
                       log_odds = log(prop_response / (1 - prop_response)))
    
    if(!any(is.infinite(agg_data$log_odds))) {
      emp_plots[[predictor]] <- agg_data |>
        ggplot2::ggplot(ggplot2::aes(x = reorder(!!as.name(predictor), log_odds), y = log_odds)) +
        ggplot2::geom_col() +
        ggplot2::labs(x = predictor,
                      y = paste0("Log odds of ", response, " = ", ref_cat),
                      title = "Empirical Logit Plot") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
    } else {
      warning(paste0("Complete separation for levels of ", predictor, 
                     "."))
    }
  }
  
  if(length(emp_plots) > 1) {
    return(emp_plots[colnames(my_data)[-c(1)]])
  } else {
    return(emp_plots[[1]])
  }
}
