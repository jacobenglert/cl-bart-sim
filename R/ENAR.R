# Figure for ENAR 2024 Student Paper Contest

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(patchwork)

# Load Results ------------------------------------------------------------

# Simulation dates
c_sim_date <- '01AUG2023'
f_sim_date <- '31JUL2023'

# Simulate variable importance plots
c_plt <- read_rds(here::here('Figures','CART', c_sim_date, 'cart-var-imp.rds')) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(size = 10)) +
  labs(title = 'CART Simulation')
f_plt <- read_rds(here::here('Figures','friedman', f_sim_date, 'friedman-var-imp.rds')) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(size = 10)) +
  labs(title = 'Friedman Simulation')

# Concatenate plots
c_plt + f_plt
ggsave(here::here('Figures','ENAR','sim-var-imp.png'), width = 6.5, height = 4.5)
