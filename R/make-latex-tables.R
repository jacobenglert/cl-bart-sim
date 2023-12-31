# Program Name: make-latex-tables.R
# Author: Jacob Englert
# Date: 25JUL2023
# Description: Convert .csv files to nice LaTeX code to be copied to Overleaf.


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(knitr)
library(kableExtra)


# Tables ------------------------------------------------------------------

# Simulate dates
cart_sim_date <- '01AUG2023'
friedman_sim_date <- '31JUL2023'

# CART summary statistics - overall
read_csv(here::here('Tables','CART', cart_sim_date, 'cart-summary-stats.csv'),
         show_col_types = FALSE) |>
  filter(p == 10 & corstr == 'AR-1' & alpha_rho == 0.95) |>
  filter(subset == 'Overall') |>
  select(model, num_trees, bias, rpehe, coverage, width) |>
  kbl(col.names = c('Model','Trees','Bias','REPHE','Coverage','Width'),
      booktabs = TRUE, format = 'latex') |>
  collapse_rows(1, latex_hline = 'full', valign = 'middle')
  
# CART summary statistics - grouped
read_csv(here::here('Tables','CART', cart_sim_date, 'cart-summary-stats.csv'),
           show_col_types = FALSE) |>
    filter(p == 10 & corstr == 'AR-1' & alpha_rho == 0.95) |>
  filter(subset != 'Overall') |>
    select(subset, model, num_trees, bias, rpehe, coverage, width) |>
    kbl(col.names = c('Group','Model','Trees','Bias','REPHE','Coverage','Width'),
        booktabs = TRUE, format = 'latex') |>
    collapse_rows(1:2, latex_hline = 'full', valign = 'middle')

# Friedman summary statistics - overall
read_csv(here::here('Tables','friedman', friedman_sim_date,'friedman-summary-stats.csv'),
         show_col_types = FALSE) |>
  filter(p == 10 & alpha_rho == 0.95) |>
  select(model, num_trees, bias, rpehe, coverage, width) |>
  kbl(col.names = c('Model','Trees','Bias','REPHE','Coverage','Width'),
      booktabs = TRUE, format = 'latex') |>
  collapse_rows(1, latex_hline = 'full', valign = 'middle')
