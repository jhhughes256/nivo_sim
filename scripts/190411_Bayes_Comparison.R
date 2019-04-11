# Comparing 1000 pop run of dosing_protocol study
# -----------------------------------------------------------------------------
# Was run for 1000 population for flat dosing and model based dosing. Checking
#   for errors and any problems before further runs are made.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(magrittr)
  library(ggplot2)
  library(scales)
  
# Source external scripts
  source("scripts/functions_utility.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Read in data
  trough_flat_df <- readr::read_rds("rds/flat_dosing_1000.rds")
  output_bayes_df <- readr::read_rds("rds/model_based_1000mod1.rds")
  
# Process bayes output
  trough_bayes_df <- output_bayes_df %>%
    dplyr::select(ID, bayes) %>%
    tidyr::unnest()
  
# Check for patients with fatal errors in Bayes estimation
  error_status <- output_bayes_df %>%
    dplyr::mutate(
      out = purrr::map(bayes, function(x) {
        any(is.na(dplyr::slice(x, -dim(x)[1])))
      } 
    )) %>%
    dplyr::pull(out) %>%
    unlist()
  
  sum(error_status)
  which(error_status)
  
# How many of these had a fatal error from the first fit
  trough_bayes_df %>%
    dplyr::filter(ID %in% which(error_status) & time == 28) %>%
    dplyr::pull(DV) %>% is.na() %>% sum()
  # all but one!

# What happened to that other one?
  trough_bayes_df %>%
    dplyr::filter(ID %in% which(error_status)[2])
  # likely had convergence from Bayes estimate for full model, but model that
  #   was fit resulted in NaN for DV
  
# Is there some strangeness with doses defaulting to 240?
  trough_bayes_df %>%
    dplyr::filter(time > 14) %>% 
    {dplyr::slice(., which(.$amt == 240))} %>% 
    {unique(.$ID)} %>% 
    length()
  # 151 patients have extra 240 mg doses other than the induction
  # this means they failed to optimise the dose
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# So now lets observe some plots
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Plot summary of patient concentrations
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = time, y = DV), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = DV), geom = "ribbon", data = trough_flat_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "red", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = DV), geom = "line", fun.y = "median", data = trough_bayes_df,
    colour = "blue", size = 1)
  p <- p + stat_summary(aes(x = time, y = DV), geom = "line", fun.y = "median", data = trough_flat_df,
    colour = "red", size = 1)
  p <- p + geom_hline(yintercept = 2.5, linetype = "dashed")
  p <- p + labs(x = "Time (days)", y = "Concentration (mg/mL)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  p
  
# Plot summary of patient tumour size
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = time, y = TUM), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = TUM), geom = "ribbon", data = trough_flat_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "red", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = TUM), geom = "line", fun.y = "median", data = trough_bayes_df,
    colour = "blue", size = 1)
  p <- p + stat_summary(aes(x = time, y = TUM), geom = "line", fun.y = "median", data = trough_flat_df,
    colour = "red", size = 1)
  p <- p + geom_hline(yintercept = 54.6, linetype = "dashed")
  p <- p + labs(x = "Time (days)", y = "Tumour Size (mm)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  p
  
# 
  p <- NULL
  p <- ggplot(data = trough_flat_df[trough_flat_df$time %in% c(14, 112, 350),])
  p <- p + geom_point(aes(x = DV, y = TUM), shape = 1)
  p <- p + scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = comma)
  p <- p + coord_cartesian(xlim = c(0.1, 1000))
  p <- p + facet_wrap(~time)
  p

# Plot summary of patient doses
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI80lo", fun.ymax = "CI80hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI70lo", fun.ymax = "CI70hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI60lo", fun.ymax = "CI60hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI50lo", fun.ymax = "CI50hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI40lo", fun.ymax = "CI40hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI30lo", fun.ymax = "CI30hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI20lo", fun.ymax = "CI20hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI10lo", fun.ymax = "CI10hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "line", fun.y = "median", data = trough_bayes_df,
    colour = "blue", size = 1)
  p <- p + geom_hline(yintercept = 240, linetype = "dashed")
  p <- p + labs(x = "Time (days)", y = "Amount Administered (mg)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  # p <- p + facet_wrap(~ ECOG)
  p
  
# Plot summary of patient doses per kg
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI80lo", fun.ymax = "CI80hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI70lo", fun.ymax = "CI70hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI60lo", fun.ymax = "CI60hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI50lo", fun.ymax = "CI50hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI40lo", fun.ymax = "CI40hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI30lo", fun.ymax = "CI30hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI20lo", fun.ymax = "CI20hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon", data = trough_bayes_df,
    fun.ymin = "CI10lo", fun.ymax = "CI10hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "line", fun.y = "median", data = trough_bayes_df,
    colour = "blue", size = 1)
  p <- p + stat_summary(aes(x = time, y = 240/BWT), geom = "line", fun.y = "median", data = trough_flat_df,
    colour = "red", size = 1)
  p <- p + labs(x = "Time (days)", y = "Amount Administered (mg/kg)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  # p <- p + facet_wrap(~ ECOG)
  p