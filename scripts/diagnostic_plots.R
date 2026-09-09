#################
# Title: Log Transformation
# Author: Andrew DiLernia
# Date: 09/01/2021
# Purpose: Display use of log transformation for right-skewed distributions
#################

library(tidyverse)
library(scales)
library(shadowtext)
library(ggfortify)
library(patchwork)
library(car)

# Generating right-skewed values
set.seed(1994)

setwd("/Users/dilerand/Library/CloudStorage/GoogleDrive-asdilernia@gmail.com/My Drive/dilerand@mail.gvsu.edu 2022-11-21 21:19/GVSU/DSA 220/Activities/05-Simple-Linear-Regression/graphics")

# Creating custom ggplot theme
# Original source code from https://github.com/alex23lemm/theme_fivethirtyeight/blob/master/theme_fivethirtyeight.R
sas_scatter <- function(myGG = NULL) {
  myGG + geom_point(shape = 1, color = "#05379b") + 
    theme_bw() %+replace%
    theme(panel.grid = element_blank(), 
          text = element_text(face = "bold"))
}

# Testing theme
N <- 100
set.seed(1994)
myGG <- data.frame(X = rnorm(N), Y = rnorm(N)) |> 
  ggplot(aes(x = X, y = Y)) + 
  labs(title = "Example SAS Scatter Plot")

sas_scatter(myGG)

# Creating raw and log transformed survival times
survivalTimes <- data.frame(Survival = rexp(n = 1000, rate = 1 / 400) + 7) |> 
  mutate(LogSurvival = log(Survival)) |> 
  pivot_longer(cols = everything(), values_to = "Time",
               names_to = "Type") |> 
  mutate(Type = factor(ifelse(Type == "Survival", "Survival in Days",
                       "Log(Survival in Days)"), 
                       levels = c("Survival in Days", "Log(Survival in Days)")))

# Creating scatter plot for log survival times
survivalTimes |> ggplot(aes(x = Time)) + 
  geom_histogram(fill = "#90B0D9", color = "black") + 
  labs(x = "Survival time", y = "Frequency",
       title = "Survival Times: Raw and Log Transformed") +
  facet_grid(. ~ Type, scales = "free_x") + ggthemes::theme_clean() +
  theme(text = element_text(face = "bold"))

# Residual by Fitted Values Plot ------------------------------------------

# Generating homoskedastic and independent errors
set.seed(1994)
N <- 300
slopeParam <- 6
sinFreq <- 2
residuals_data <- data.frame(Fitted.Value = 1:N,
      HomoInd = rnorm(n = N),
      HeteroFanInd = map_dbl(.x = 1:N,
      .f = function(x){rnorm(n = 1, mean = 0, sd = 1.006^x)}),
      HeteroFunnelInd = map_dbl(.x = 1:N,
                                .f = function(x){rnorm(n = 1, mean = 0, sd = 1.006^(-x))}),
      HeteroBowInd = map_dbl(.x = 1:N,
                                .f = function(x){rnorm(n = 1, mean = 0, sd = 1.00008^((x-N/2)^2))}),
      HomoQuadDep = 0.50*slopeParam*(1:N - N/2)^2 / (N - N/2)^2,
      HomoPosLinDep = slopeParam*1:N / N,
      HomoSinDep = sin(seq(0, sinFreq*2*pi, length.out = N))) |> 
  mutate(FanPosLinDep = HomoPosLinDep + HeteroFanInd,
         HomoInvQuadDep = -HomoQuadDep,
         FunnelPosLinDep = HomoPosLinDep + HeteroFunnelInd,
         HomoNegLinDep = -HomoPosLinDep + HomoInd,
         HomoPosLinDep = HomoPosLinDep + HomoInd,
         FanQuadDep = HomoQuadDep + HeteroFanInd,
         FunnelQuadDep = HomoQuadDep + HeteroFunnelInd,
         HomoQuadDep = HomoQuadDep + HomoInd,
         BowSinDep = HomoSinDep + HeteroBowInd,
         HomoSinDep = HomoSinDep + HomoInd) |> 
  pivot_longer(cols = -Fitted.Value, names_to = "Type",
               values_to = "Residual") |> 
  mutate(Homogeneity = ifelse(str_detect(Type, pattern = "Homo"), 
                           "Homoskedastic", "Heteroskedastic"),
         Dependence = ifelse(str_detect(Type, pattern = "Ind"), 
                           "Independent", "Dependent")) |>
  group_by(Type) |> mutate(Residual = scale(Residual, scale = FALSE))

# Function for plotting residuals
plot_residuals <- function(types, errors = residuals_data) {
  errors <- errors |> 
    filter(Type %in% types)
  
  label_data <- errors |> 
    group_by(Homogeneity, Dependence) |> 
    summarize(Residual = 0.98*max(Residual),
              Fitted.Value = 1) |> 
    ungroup() |> 
    mutate(Label = toupper(letters)[1:4],
      Residual = case_when(Homogeneity == "Homoskedastic" ~ 
                        sum((Homogeneity == "Homoskedastic")*Residual) / 2,
                        TRUE ~ sum((Homogeneity == "Heteroskedastic")*Residual) / 2))
    
  errors |> 
    ggplot(aes(x = Fitted.Value, y = Residual)) + 
    geom_hline(yintercept = 0, color = "gray", linetype = "dotted") + 
    geom_point(color = "#56B4E9", alpha = 0.6, size = 1) + 
    facet_grid(Homogeneity ~ Dependence, scales = "free_y") + 
    scale_y_continuous(breaks = 0) +
    geom_shadowtext(data = label_data, color = "white", 
              mapping = aes(x = Fitted.Value,
                            y = Residual,
                            label = Label)) +
    labs(x = "Fitted values",
         y = "Residuals") + theme_bw() + 
    theme(panel.grid = element_blank(), 
          text = element_text(face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

# Plotting Linear Dependencies
residual_by_predicted_linear <- plot_residuals(c("HeteroFanInd", "HomoInd", "FanPosLinDep", "HomoPosLinDep"))
ggsave(filename = "residual_by_predicted_plot_linear.png",
       width = 1605, height = 996, 
       units = "px", plot = residual_by_predicted_linear)

# Plotting Quadratic Dependencies
residual_by_predicted_quadratic <- plot_residuals(c("HomoInd", "HeteroFunnelInd", "HomoQuadDep", "FunnelQuadDep"))
ggsave(filename = "residual_by_predicted_plot_quadratic.png",
       width = 1605, height = 996, 
       units = "px", plot = residual_by_predicted_quadratic)

# Plotting Oscillating Dependencies
residual_by_predicted_sin <- plot_residuals(c("HomoInd", "HeteroBowInd", "HomoSinDep", "BowSinDep"))
ggsave(filename = "residual_by_predicted_plot_sin.png",
       width = 1605, height = 996, 
       units = "px", plot = residual_by_predicted_sin)

# Histograms and QQ Plots -------------------------------------------------

# Generating data
set.seed(1994)
N <- 120
normal_data <- data.frame(Residual = c(rnorm(n = N), 
                                      rexp(n = N, rate = 2),
                                      -rexp(n = N, rate = 2), 
                                      rnorm(n = N/2, mean = -2.5), 
                                      rnorm(n = N/2, mean = 2.5),
                                      rcauchy(n = N/2, scale = 0.01), rnorm(n = N/2),
                                      round(rnorm(n = N, sd = 2))),
                         Type = c(rep("Normal", N), rep("Right-skewed", N),
                                  rep("Left-skewed", N), rep("Bimodal", N),
                         rep("Peaked", N), rep("Discrete", N))) |> 
  group_by(Type) |> mutate(Residual = scale(Residual)) |> ungroup() |> 
  mutate(Type = fct_relevel(Type, "Normal", "Right-skewed", "Left-skewed",
                            "Bimodal", "Peaked", "Discrete"))

# Modifying cauchy to not be so extreme
normal_data <- normal_data |> 
  mutate(Residual = case_when(Type == "Peaked" ~ Residual / sqrt(abs(Residual)),
                              TRUE ~ Residual))

# Function for plotting histograms of residuals
plot_histograms <- function(plot_data = normal_data |> filter(Type %in% c("Normal", "Bimodal", "Peaked")),
                        letter_set = 1, letters_x = -2.3, letters_y = 0.80) {
  
  label_data <- plot_data |> 
    group_by(Type) |> 
    summarize(Density = letters_y) |> 
    ungroup() |> 
    mutate(Residual = 0.80*min(plot_data$Residual), 
           Label = toupper(letters)[1:3 + 3*(letter_set-1)])
  
  plot_data |> 
    ggplot(aes(x = Residual)) + 
    geom_histogram(aes(y = after_stat(density)), fill = "#56B4E9", color = "white") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "#009E73") +
    facet_grid(. ~ Type, scales = "free_y") + 
    geom_shadowtext(data = label_data, color = "white", 
                    mapping = aes(x = letters_x,
                                  y = Density,label = Label)) +
    scale_x_continuous(breaks = seq(-6, 6, by = 2),
                       limits = range(normal_data$Residual), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, .01))) +
    labs(y = "Density") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          text = element_text(face = "bold"))
}

# Plotting histograms
ggsave(plot = normal_data |> filter(Type %in% c("Normal", "Bimodal", "Peaked")) |> 
         plot_histograms(letter_set = 1,letters_y = 1.6),
       filename = "histograms_bimodal_peaked.png",
       width = 2000, height = 996, 
       units = "px")

ggsave(plot = normal_data |> filter(Type %in% c("Left-skewed", "Right-skewed", "Discrete")) |> 
         plot_histograms(letter_set = 1, letters_y = 0.82),
       filename = "histograms_skewed_discrete.png",
       width = 2000, height = 996, 
       units = "px")

# Function for plotting QQ plots of residuals
plot_qq <- function(plot_data = normal_data, letter_set = 1) {
  
  label_data <- plot_data |> 
    group_by(Type) |> 
    summarize(Residual = max(plot_data$Residual), 
              xpos = 0.95*min(plot_data$Residual)) |> 
    ungroup() |> 
    mutate(Label = toupper(letters)[1:3 + 3*(letter_set-1)])
  
  plot_data |> ggplot(aes(sample = Residual)) + 
    stat_smooth(aes(x = Residual, y = Residual), method=lm, 
                se=FALSE, formula=y~x-1, color = "grey", fullrange=TRUE,
                linetype = "dotted", linewidth = 0.7) +
    stat_qq(color = "#56B4E9", alpha = 0.6, size = 1) +
    facet_grid(. ~ Type) + 
    geom_shadowtext(data = label_data, color = "white", 
                    mapping = aes(x = xpos,
                                  y = Residual,label = Label)) +
    scale_x_continuous(limits = range(plot_data$Residual)) +
    scale_y_continuous(limits = range(plot_data$Residual)) +
    labs(x = "Theoretical quantiles",
         y = "Standardized residuals") + 
    theme_bw() + theme(panel.grid = element_blank(), 
                       text = element_text(face = "bold"))
}

# Plotting QQ plots
ggsave(plot = normal_data |> filter(Type %in% c("Normal", "Bimodal", "Peaked")) |> 
         plot_qq(letter_set = 1),
       filename = "qqs_bimodal_peaked.png",
       width = 2000, height = 996, 
       units = "px")

ggsave(plot = normal_data |> filter(Type %in% c("Left-skewed", "Right-skewed", "Discrete")) |> 
         plot_qq(letter_set = 1),
       filename = "qqs_skewed_discrete.png",
       width = 2000, height = 996, 
       units = "px")

# Making plots for exam
exam_folder <- "/Users/dilerand/Library/CloudStorage/GoogleDrive-asdilernia@gmail.com/My Drive/dilerand@mail.gvsu.edu 2022-11-21 21:19/GVSU/DSA 220/Exams/"
for(data_type in c("Normal", "Bimodal", "Peaked", "Left-skewed", "Right-skewed", "Discrete")) {

  plot_data <- normal_data |> 
    dplyr::filter(Type == data_type)
  
  my_qq <- plot_data |> 
  ggplot(aes(sample = Residual)) + 
    stat_smooth(aes(x = Residual, y = Residual), method=lm, 
                se=FALSE, formula=y~x-1, color = "grey", fullrange=TRUE,
                linetype = "dotted", linewidth = 0.7) +
    stat_qq(color = "black", alpha = 0.6, size = 1) +
  scale_x_continuous(limits = range(plot_data$Residual)) +
  scale_y_continuous(limits = range(plot_data$Residual)) +
    labs(x = "Theoretical quantiles",
         y = "Standardized residuals") + 
    theme_bw() + theme(panel.grid = element_blank(), 
                       text = element_text(face = "bold"))

my_histogram <- plot_data |>  
  ggplot(aes(x = Residual)) + 
  geom_histogram(aes(y = after_stat(density)), fill = "#56B4E9", color = "white") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "#009E73") +
  scale_x_continuous(breaks = seq(-4, 4, by = 2),
                     limits = range(normal_data$Residual), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .01))) +
  labs(y = "Density") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        text = element_text(face = "bold"))

# Saving plots
ggsave(plot = my_qq,
       filename = paste0(exam_folder, "qqs_", data_type, ".png"),
       width = 2000, height = 996, 
       units = "px")

# Saving plots
ggsave(plot = my_histogram,
       filename = paste0(exam_folder, "hist_", data_type, ".png"),
       width = 2000, height = 996, 
       units = "px")
}

# Outliers and Influential Points -----------------------------------------

# Function for plotting variance inflation factor (VIF) values 
vif_plot <- function(model_fit) {
  
  vifs <- car::vif(model_fit)
  
  if("GVIF^(1/(2*Df))" %in% colnames(vifs)) {
    vifs <- vifs[, "GVIF^(1/(2*Df))"]
  }
  
  preds <- names(vifs)
  
  vifGG <- tibble(Predictor = preds,
                  VIF = vifs) |> 
    dplyr::mutate(Predictor = fct_reorder(Predictor, -VIF)) |> 
    ggplot(aes(x = VIF, y = Predictor)) + 
    geom_segment(aes(xend = VIF, x = 0, 
                     yend = Predictor, y = Predictor)) +
    geom_vline(xintercept = 5, linetype = "dotted") +
    geom_vline(xintercept = 10, linetype = "dotted") +
    scale_x_continuous(limits = c(0, max(vifs)*1.1),
                       expand = expansion(mult = c(0, 0.10))) +
    geom_point(size = 3, color = "#56B4E9") +
    labs(title = "Variance inflation values")
  
  return(suppressWarnings(print(vifGG)))
}

# Set ggplot theme for visualizations
theme_set(ggthemes::theme_few())

set.seed(1994)

# Simulate data 
N <- 27
beta1 <- 2
beta0 <- -3
outliers_influentials <- tibble(X1 = c(rnorm(n = N-2, sd = 2), -6, 5)) |> 
  mutate(ID = 1:n(),
         Y = beta0 + X1*beta1 + c(rnorm(n = N-2, sd = 1.5), 6*(beta1), 0),
         Y = case_when(ID == 5 ~ -11, 
                       ID == N - 1 ~ -4, 
                       TRUE ~ Y),
         X2 = X1*beta1 + c(rnorm(n = N, sd = 2)),
         X3 = rnorm(n = N, sd = 2))

# Fitting linear model
mlr_model <- lm(Y ~ X1 + X2 + X3, data = outliers_influentials)

# Cook's Distance plot
cd_plot <- autoplot(mlr_model, which = 4, 
               label.repel = TRUE)@plots[[1]] +
  labs(x = "Observation",
       y = "Cook's D",
       title = "") +
  geom_hline(yintercept = 0, linetype = "solid",
             color = "black") +
  geom_point(shape = 1, color = "#05379b")

# Removing text labels
cd_plot[["layers"]][[2]] <- NULL

# Studentized residuals by leverage
sr_leverage_plot <- autoplot(mlr_model, which = 5, label.repel = TRUE)@plots[[1]] +
  geom_point(shape = 1, color = "#05379b") +
  labs(y = "RStudent",
       title = "") +
  geom_hline(yintercept = 2, linetype = "solid",
             color = "black") +
  geom_hline(yintercept = -2, linetype = "solid",
             color = "black") 
  
# Removing text labels and solid points
sr_leverage_plot[["layers"]][c(1, 2, 3, 5)] <- NULL

combinedGG <- (sr_leverage_plot / cd_plot)

# Variance inflation values
vif_plot(mlr_model)

# Saving for STA 216 Guide document
ggsave(combinedGG, filename = "guide_diagnostics.png",
       width = 15, height = 15, units = "cm")
