# Set up environment
  library(magrittr)
  library(ggplot2)
  setwd("E:/Hughes/Git/nivo_sim/scripts/dosing_protocols")
  source("functions_utility.R")

# Read in data
  data1_df <- dplyr::bind_rows(
    readr::read_rds("model_based_mod1_1.rds"),
    readr::read_rds("model_based_mod1_2.rds"),
    readr::read_rds("model_based_mod1_3.rds"),
    readr::read_rds("model_based_mod1_4.rds"),
    readr::read_rds("model_based_mod1_5.rds"),
    readr::read_rds("model_based_mod1_6.rds"),
    readr::read_rds("model_based_mod1_7.rds"),
    readr::read_rds("model_based_mod1_8.rds"),
    # readr::read_rds("model_based_mod1_9.rds"),
    readr::read_rds("model_based_mod1_10.rds")
  )
  data2_df <- dplyr::bind_rows(
    readr::read_rds("model_based_mod2_1.rds"),
    readr::read_rds("model_based_mod2_2.rds"),
    readr::read_rds("model_based_mod2_3.rds"),
    readr::read_rds("model_based_mod2_4.rds"),
    readr::read_rds("model_based_mod2_5.rds"),
    readr::read_rds("model_based_mod2_6.rds"),
    readr::read_rds("model_based_mod2_7.rds"),
    readr::read_rds("model_based_mod2_8.rds"),
    # readr::read_rds("model_based_mod2_9.rds"),
    readr::read_rds("model_based_mod2_10.rds")
  )
  trough_flat_df <- readr::read_rds("flat_dosing.rds")

# Analyse different datasets
  trough_data1_df <- data1_df %>%
    dplyr::select(ID, bayes) %>%
    tidyr::unnest()
  trough_data2_df <- data2_df %>%
    dplyr::select(ID, bayes) %>%
    tidyr::unnest()

# Check for patients with fatal errors in Bayes estimation
  data1_errstat <- data1_df %>%
    dplyr::mutate(
      out = purrr::map(bayes, function(x) {
        any(is.na(dplyr::slice(x, -dim(x)[1])))
      }
    )) %>%
    dplyr::pull(out) %>%
    unlist()

  sum(data1_errstat)
  data1_errwhich <- which(data1_errstat)
  data1_errwhich <- data1_errwhich %>% {dplyr::if_else(. > 800, . + 100, . + 0)}

# How many of these had a fatal error from the first fit
  trough_data1_df %>%
    dplyr::filter(ID %in% data1_errwhich & time == 28) %>%
    dplyr::pull(DV) %>% is.na() %>% sum()
  # all of them

# Check for other dataset
  data2_errstat <- data2_df %>%
    dplyr::mutate(
      out = purrr::map(bayes, function(x) {
        any(is.na(dplyr::slice(x, -dim(x)[1])))
      }
    )) %>%
    dplyr::pull(out) %>%
    unlist()

  sum(data2_errstat)
  data2_errwhich <- which(data2_errstat)
  data2_errwhich <- data2_errwhich %>% {dplyr::if_else(. > 800, . + 100, . + 0)}

# How many of these had a fatal error from the first fit
  trough_data2_df %>%
    dplyr::filter(ID %in% data2_errwhich & time == 28) %>%
    dplyr::pull(DV) %>% is.na() %>% sum()
  # all of them!

# Are doses defaulting 240 much?
  trough_data2_df %>%
    dplyr::filter(time > 14) %>%
    {dplyr::slice(., which(.$amt == 240))} %>%
    {unique(.$ID)} %>%
    length()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# So now lets observe some plots
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Plot summary of patient concentrations
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = time, y = DV), geom = "ribbon", data = trough_data1_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = DV), geom = "ribbon", data = trough_flat_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "red", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = DV), geom = "line", fun.y = "median", data = trough_data1_df,
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
  p <- p + stat_summary(aes(x = time, y = TUM), geom = "ribbon", data = trough_data1_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "blue", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = TUM), geom = "ribbon", data = trough_flat_df,
    fun.ymin = "CI90lo", fun.ymax = "CI90hi", fill = "red", alpha = 0.1)
  p <- p + stat_summary(aes(x = time, y = TUM), geom = "line", fun.y = "median", data = trough_data1_df,
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
  p <- ggplot(data = trough_data1_df)
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
  p <- ggplot(data = trough_data2_df)
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
