# Program Name: set-params.R
# Author: Jacob Englert
# Date: 21JUL2023
# Description: Generate parameter tables for CL-BART simulations.

# Load Packages -----------------------------------------------------------
library(tidyverse)

# CART Simulation ---------------------------------------------------------
rho <- 0.6
corstr <- c("'exchangeable'", "'independent'", "'AR-1'")
p <- c(10, 20)
alpha_rho <- c(0.95, 0.3)
beta_rho <- c(2, 3)
seed <- 1:100
num_trees <- c(1, 2, 5, 10)

params <- crossing(rho, corstr, p, nesting(alpha_rho, beta_rho), seed, num_trees) |>
  mutate(iter = 5000, thin = 5, warmup = 2500,
         sigma2_beta = 1, sigma2_beta_update_freq = 50, beta_acc_prob = 0.35,
         n_min = 5,
         alpha_sigma = 1, beta_sigma = 1,
         moves = "c('grow','prune','change')", move_probs = 'c(0.3,0.3,0.4)',
         ID = row_number()) |>
  select(ID, everything())

write_csv(params, here::here('Params','cart-params.csv'))

# Friedman Simulation -----------------------------------------------------
p <- c(10, 20)
alpha_rho <- c(0.95, 0.3)
beta_rho <- c(2, 3)
seed <- 1:100
num_trees <- c(5, 10, 25, 50)

params <- crossing(p, nesting(alpha_rho, beta_rho), seed, num_trees) |>
  mutate(iter = 6000, thin = 5, warmup = 3000,
         sigma2_beta = 1, sigma2_beta_update_freq = 50, beta_acc_prob = 0.35,
         n_min = 5,
         alpha_sigma = 1, beta_sigma = 1,
         moves = "c('grow','prune','change')", move_probs = 'c(0.3,0.3,0.4)',
         ID = row_number()) |>
  select(ID, everything())

write_csv(params, here::here('Params','friedman-params.csv'))
