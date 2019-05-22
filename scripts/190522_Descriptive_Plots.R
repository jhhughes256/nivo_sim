# Experimentation with different descriptive plots for manuscript
# -----------------------------------------------------------------------------
# For publication of this work, informative figures and results will need to
#   be presented. Aspects of the work that are of interest are:
# - How well do these TDM methods keep patients within the range?
# - How much money do these TDM methods save?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(magrittr)
  library(ggplot2)  # Graphical package
  
# Source functions utility
  source("scripts/functions_utility.R")
  
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Set colourblind palette
  cbPalette <- data.frame(
		grey = "#999999",
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)
  
# Read in data
  flat_df <- readr::read_rds("output/flat_dose.rds")
  prop_df <- readr::read_rds("output/proportional_tdm.rds")
  step_df <- readr::read_rds("output/stepwise_tdm_alt.rds")
  
# Convert data to same format
  flat_df <- flat_df %>%
    tibble::add_column(method = "flat")
  
  prop_df <- prop_df %>%
    dplyr::select(ID, bayes) %>%
    tidyr::unnest() %>%
    dplyr::select(-ID1) %>%
    tibble::add_column(method = "prop")
    
  step_df <- step_df %>%
    dplyr::select(ID, data) %>%
    tidyr::unnest() %>%
    dplyr::select(-ID1) %>%
    tibble::add_column(method = "step")
  
# Bind into plot dataset and calculate cumulative dose
  plot_df <- dplyr::bind_rows(flat_df, prop_df, step_df) %>%
    dplyr::group_by(ID, method) %>% tidyr::nest() %>%  # nest_by
    dplyr::mutate(data = purrr::map(data, function(df) {
      dplyr::mutate(df, 
        cumcost = cumsum(amt),  # *dose_cost),
        relcost = cumsum(amt)/c(cumsum(rep(240, 26)), 240*26)  # *dose_cost)
      )
    })) %>% tidyr::unnest() %>%
# Determine time in range by setting a flag
    dplyr::mutate(cflag = dplyr::if_else(DV >= 2.5 & DV <= 5, 1, 0))
  
  inrange_df <- plot_df %>%
    dplyr::group_by(method, time) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(df) {
      sum(df$cflag, na.rm = T)/dim(df)[1]
    })) %>% tidyr::unnest()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Evaluation Plots
# Trough Concentration-time plots
  palette1 <- with(cbPalette, c(red, blue, green))
  p <- NULL
  p <- ggplot(aes(x = time, y = DV, fill = method), data = plot_df)
  p <- p + ggtitle("Trough Concentrations over Time")
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI90lo", 
    fun.ymax = "CI90hi", alpha = 0.1)
  p <- p + stat_summary(aes(colour = method), geom = "line", fun.y = "median", 
    size = 1)
  p <- p + geom_hline(yintercept = 2.5, linetype = "dashed")
  p <- p + scale_colour_manual(values = palette1)
  p <- p + scale_fill_manual(values = palette1)
  p <- p + labs(x = "Time (days)", y = "Concentration (mg/mL)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  p
  
# Dose-time plots
  p <- NULL
  p <- ggplot(aes(x = time, y = amt, fill = method), 
    data = dplyr::filter(plot_df, time != 364))
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI90lo", 
    fun.ymax = "CI90hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI80lo", 
    fun.ymax = "CI80hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI70lo", 
    fun.ymax = "CI70hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI60lo", 
    fun.ymax = "CI60hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI50lo", 
    fun.ymax = "CI50hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI40lo", 
    fun.ymax = "CI40hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI30lo", 
    fun.ymax = "CI30hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI20lo", 
    fun.ymax = "CI20hi", alpha = 0.05)
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI10lo", 
    fun.ymax = "CI10hi", alpha = 0.05)
  p <- p + stat_summary(aes(colour = method), geom = "line", fun.y = "median", 
    size = 1)
  p <- p + facet_wrap(~method)
  p
  
# Time in Range plots
# Has issues with old induction script, needs resampling to work well!!
  p <- NULL
  p <- ggplot(data = inrange_df)
  p <- p + ggtitle("Proportion of Subjects in Therapeutic Range")
  p <- p + geom_line(aes(x = time, y = data*100, colour = method), size = 1)
  p <- p + scale_colour_manual(values = palette1)
  p <- p + labs(x = "Time (days)", y = "Subjects in Therapeutic Range (%)")
  p
  
# Cost Plots
  p <- NULL
  p <- ggplot(aes(x = time, y = cumcost, fill = method), 
    data = dplyr::filter(plot_df, time != 364))
  p <- p + ggtitle("Relative Cost per Patient")
  p <- p + stat_summary(geom = "ribbon", fun.ymin = "CI90lo", 
    fun.ymax = "CI90hi", alpha = 0.1)
  p <- p + stat_summary(aes(colour = method), geom = "line", fun.y = "median", 
    size = 1)
  p <- p + scale_colour_manual(values = palette1)
  p <- p + scale_fill_manual(values = palette1)
  p <- p + labs(x = "Time (days)", y = "Cost ($)")
  # p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
  p
  
# Cost table
  cost_tbl <- plot_df %>%
    dplyr::group_by(method) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(df) {
      data.frame(
        t112 = dplyr::filter(df, time == 112) %>%
          dplyr::pull(cumcost) %>% sum(),
        t224 = dplyr::filter(df, time == 224) %>%
          dplyr::pull(cumcost) %>% sum(),
        t336 = dplyr::filter(df, time == 336) %>%
          dplyr::pull(cumcost) %>% sum()
      )
    })) %>% tidyr::unnest()
  