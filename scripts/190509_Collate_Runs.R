# Set up environment
  library(magrittr)
  library(ggplot2)
  setwd("E:/Hughes/Git/nivo_sim/scripts/dosing_protocols")
  source("functions_utility.R")

# Read in data
  id_run <- paste0("id/id", (1:100), ".rds")
  trough_bayes_df <- purrr::map_dfr(id_run, readr::read_rds)
  trough_flat_df <- readr::read_rds("flat_dosing.rds") %>%
    dplyr::filter(ID != 80)

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

# Plot summary of patient doses
  p <- NULL
  p <- ggplot(data = trough_bayes_df)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI80lo", fun.ymax = "CI80hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI70lo", fun.ymax = "CI70hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI60lo", fun.ymax = "CI60hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI50lo", fun.ymax = "CI50hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI40lo", fun.ymax = "CI40hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI30lo", fun.ymax = "CI30hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI20lo", fun.ymax = "CI20hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "ribbon",
    fun.ymin = "CI10lo", fun.ymax = "CI10hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt), geom = "line", fun.y = "median",
    colour = "blue", size = 1)
  p <- p + geom_hline(yintercept = 240, linetype = "dashed")
  p <- p + labs(x = "Time (days)", y = "Amount Administered (mg)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  # p <- p + facet_wrap(~ ECOG)
  p

# Plot summary of patient doses per kg
  p <- NULL
  p <- ggplot(data = trough_bayes_df)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI80lo", fun.ymax = "CI80hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI70lo", fun.ymax = "CI70hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI60lo", fun.ymax = "CI60hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI50lo", fun.ymax = "CI50hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI40lo", fun.ymax = "CI40hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI30lo", fun.ymax = "CI30hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI20lo", fun.ymax = "CI20hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "ribbon",
    fun.ymin = "CI10lo", fun.ymax = "CI10hi", fill = "blue", alpha = 0.05)
  p <- p + stat_summary(aes(x = time, y = amt/BWT), geom = "line", fun.y = "median",
    colour = "blue", size = 1)
  p <- p + stat_summary(aes(x = time, y = 240/BWT), geom = "line", fun.y = "median", data = trough_flat_df,
    colour = "red", size = 1)
  p <- p + labs(x = "Time (days)", y = "Amount Administered (mg/kg)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  # p <- p + facet_wrap(~ ECOG)
  p
