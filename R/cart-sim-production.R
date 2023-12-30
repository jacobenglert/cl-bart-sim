# Program Name: cart-sim-production.R
# Author: Jacob Englert
# Date: 21JUL2023
# Description: Code for a CL-BART simulation study where the true log
# odds-ratio is determined by a single CART model - to be ran on Emory HPC.

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(survival)
library(clogitL1)
library(clbart2)
library(parallel)

# Get Parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c('index')
index <- as.numeric(args['index'])
params <- read_csv(here::here('Params','cart-params.csv'), show_col_types = F)[,-1]
for(i in 1:ncol(params)){
  assign(names(params[index,i]), eval(parse(text = params[index, i, drop = TRUE])))
}

# Set Simulation Parameters -----------------------------------------------
set.seed(seed)                              # Random Seed
n       <- 10000                            # Initial Population Size
start   <- as.Date("2020-01-01")            # Start of follow-up period
end     <- as.Date("2022-12-31")            # End of follow-up period
dates   <- seq.Date(start, end, by = 'day') # 3-year follow-up period
t       <- length(dates)                    # Length of follow-up period

# Create moderator correlation matrix
if(corstr == 'AR-1') w_cor <- sapply(1:p, \(x) rho^(abs(x - (1:p))))
if(corstr == 'exchangeable') w_cor <- (1 - rho) * diag(p) + matrix(rho, p, p)
if(corstr == 'independent') w_cor <- diag(p)

# Define function to calculate CART-based log odds-ratio
get_CART_logOR <- function(x){
  
  t1 <- 0.8 * I(x[,1] <= 0 & x[,2] <= 0) +
    1.1 * I(x[,1] <= 0 & x[,2] > 0) + 
    1.3 * I(x[,1] > 0 & x[,3] <= 0) +
    1.5 * I(x[,1] > 0 & x[,3] > 0)
  
  return(log(t1))
}

# Simulate shared exposure time series with sinusoidal time trend
z <- rnorm(t, sin(seq(0, 2 * pi * 3, length.out = t)))

# Create exposure data frame
expos_df <- data.frame(z = z, date = dates) |>
  mutate(year = year(date), month = month(date), dow = wday(date))

# Simulate binary exposure moderators
mod_df <- (mvtnorm::rmvnorm(n, sigma = w_cor) > 0) |> 
  as.data.frame() |>
  lapply(as.numeric) |>
  as.data.frame()
colnames(mod_df) <- paste0('W', 1:p)

# Follow each individual throughout time
data <- cross_join(data.frame(ID = 1:n, mod_df), expos_df)

# Compute the heterogeneous exposure effects on the log-odds scale
w_imp <- sample(colnames(mod_df), 3)
data$tau <- get_CART_logOR(data[w_imp])

# Convert to probability scale
expit <- function(x) exp(x) / (1 + exp(x))
alpha <- -8
data$p <- expit(alpha + data$tau * data$z)

# Simulate the binary outcome
y <- rbinom(n * t, 1, data$p)

# Identify which observations experienced the event
event_idx <- which(y == 1)

# Create case-crossover structure
data_cc <- data[event_idx,] |>  # subset cases
  mutate(ID = row_number()) |>  # assign unique ID (strata) to each case
  select(-z) |>
  left_join(expos_df, by = c('year', 'month', 'dow'), 
            relationship = 'many-to-many') |>
  mutate(y = ifelse(date.x == date.y, 1, 0)) |> # create case variable
  select(-date.x) |>
  rename(date = date.y)

# Remove large objects from memory
rm(data, mod_df, expos_df)

# Extract individual elements
w_cc <- select(data_cc, starts_with('W'))
y_cc <- data_cc$y
z_cc <- data_cc$z
tau_cc <- data_cc$tau
strata_cc <- data_cc$ID


# Fit "Oracle" Conditional Logistic Regression Model ----------------------
w_cc_imp <- data_cc[w_imp]
w_cc_true <- cbind(w_cc_imp[,1] <= 0 & w_cc_imp[,2] <= 0,
                   w_cc_imp[,1] <= 0 & w_cc_imp[,2] > 0,
                   w_cc_imp[,1] > 0 & w_cc_imp[,3] <= 0,
                   w_cc_imp[,1] > 0 & w_cc_imp[,3] > 0) |>
  as.numeric() |>
  matrix(ncol = 4)

fit_clogit <- clogit(y_cc ~ w_cc_true : z_cc + strata(strata_cc))

# Fit the LASSO model -----------------------------------------------------
w_cc_2 <- as.data.frame(model.matrix(~ .^2, data = w_cc))[,-1] * z_cc
fit_lasso <- clogitL1(x = w_cc_2, y = y_cc, strata = strata_cc)
fit_lasso_cv <- cv.clogitL1(fit_lasso, numFolds = 5)
lambda_idx <- which(fit_lasso_cv$lambda == fit_lasso_cv$minCV_lambda)
beta_lasso_cv <- fit_lasso_cv$beta[lambda_idx,]
tau_lasso <- as.numeric(as.matrix(w_cc_2) %*% beta_lasso_cv)[c(T,F,F,F)]


# Fit CL-BART model -------------------------------------------------------

fit_clbarts <- mclapply(1:5, function(s){
  clbart(w = w_cc, y = y_cc, z = z_cc, stratum = strata_cc,
         num_trees = num_trees, seed = s,
         iter = iter, thin = thin, warmup = warmup,
         sigma2_beta = sigma2_beta,
         sigma2_beta_update_freq = sigma2_beta_update_freq,
         beta_acc_prob = beta_acc_prob,
         n_min = n_min,
         moves = moves, move_probs = move_probs,
         alpha_rho = alpha_rho, beta_rho = beta_rho,
         alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)
}, mc.cores = 5)

fit_clbart <- fit_clbarts[[which.min(lapply(fit_clbarts, '[[', 'WAIC'))]]


# Summarize Results -------------------------------------------------------

# Find first index for each strata
firsts <- data_cc |> 
  mutate(idx = row_number()) |> 
  filter(row_number() == 1, .by = ID) |>
  select(idx) |>
  as_vector()

# true log odds-ratios
tau_true <- tau_cc[firsts]

# clogit estimates and intervals
tau_clogit <- as.numeric(w_cc_true %*% coef(fit_clogit))[firsts]
tau_clogit_lower <- as.numeric(w_cc_true %*% confint(fit_clogit)[,1])[firsts]
tau_clogit_upper <- as.numeric(w_cc_true %*% confint(fit_clogit)[,2])[firsts]

# clbart estimates and intervals
w_unique <- select(data_cc[firsts,], starts_with('W'))
tau_clbart_post <- lapply(fit_clbart$forests, 
                            \(f) predict_forest(f, w_unique)) |>
  do.call(what = rbind)
tau_clbart <- colMeans(tau_clbart_post)
tau_clbart_lower <- apply(tau_clbart_post, 2, \(t) quantile(t, 0.025))
tau_clbart_upper <- apply(tau_clbart_post, 2, \(t) quantile(t, 0.975))

# Function to generate simulation statistics
subset_stats <- function(x){
  
  if(length(x) == length(tau_true)) subset <- 'overall'
  else{
    tmp <- w_unique[x, w_imp] |> 
      Filter(f = \(x) var(x) == 0) |> 
      unique() 
    colnames(tmp) <- paste0('W', 1:length(w_imp))[names(w_unique[x, w_imp]) %in% names(tmp)]
    subset <- paste(paste(names(tmp), '=', tmp), collapse = ' & ')
  }
  
  est_clogit      <- mean(tau_clogit[x])
  bias_clogit     <- mean(tau_clogit[x] - tau_true[x])
  rpehe_clogit    <- sqrt(mean((tau_clogit[x] - tau_true[x])^2))
  coverage_clogit <- mean(between(tau_true[x], tau_clogit_lower[x], tau_clogit_upper[x]))
  width_clogit    <- mean(tau_clogit_upper[x] - tau_clogit_lower[x])
  
  est_lasso <- mean(tau_lasso[x])
  bias_lasso <- mean(tau_lasso[x] - tau_true[x])
  rpehe_lasso <- sqrt(mean((tau_lasso[x] - tau_true[x])^2))
  
  est_clbart      <- mean(tau_clbart[x])
  bias_clbart     <- mean(tau_clbart[x] - tau_true[x])
  rpehe_clbart    <- sqrt(mean((tau_clbart[x] - tau_true[x])^2))
  coverage_clbart <- mean(between(tau_true[x], tau_clbart_lower[x], tau_clbart_upper[x]))
  width_clbart    <- mean(tau_clbart_upper[x] - tau_clbart_lower[x])
  
  tau_bar_mean  <- mean(rowMeans(tau_clbart_post[,x]))
  tau_bar_lower <- quantile(rowMeans(tau_clbart_post[,x]), 0.025)
  tau_bar_upper <- quantile(rowMeans(tau_clbart_post[,x]), 0.975)
  tau_bar_bias  <- mean(tau_bar_mean - tau_true[x])
  tau_bar_rpehe <- sqrt(mean((tau_bar_mean - tau_true[x])^2))
  tau_bar_coverage <- mean(between(tau_true[x], tau_bar_lower, tau_bar_upper))
  tau_bar_width <- tau_bar_upper - tau_bar_lower
  
  results <- data.frame(model = c('clogit', 'LASSO', 'CL-BART', 'CL-BART-avg'),
                        subset = subset,
                        seed = seed,
                        num_trees = num_trees,
                        rho = rho,
                        corstr = corstr,
                        p = p,
                        alpha_rho = alpha_rho,
                        beta_rho = beta_rho,
                        alpha_sigma = alpha_sigma,
                        beta_sigma = beta_sigma,
                        n = length(tau_true[x]),
                        est = c(est_clogit, est_lasso, est_clbart, tau_bar_mean),
                        bias = c(bias_clogit, bias_lasso, bias_clbart, tau_bar_bias),
                        rpehe = c(rpehe_clogit, rpehe_lasso, rpehe_clbart, tau_bar_rpehe),
                        coverage = c(coverage_clogit, NA, coverage_clbart, tau_bar_coverage),
                        width = c(width_clogit, NA, width_clbart, tau_bar_width))
  
  if(subset == 'overall'){
    results2 <- data.frame(model = 'CL-BART', waic = fit_clbart$WAIC) |>
      bind_cols(t(c(
        colMeans(fit_clbart$split_props, na.rm = TRUE)[w_imp],
        colMeans(fit_clbart$split_props, na.rm = TRUE)[!(colnames(w_cc) %in% w_imp)]) |> 
          setNames(paste0('W', 1:p))))
    results <- left_join(results, results2, by = 'model')
  }
  
  return(results)
}

overall_results <- subset_stats(x = 1:length(tau_true))
subset_results <- sapply(unique(tau_true), 
                         \(x) subset_stats(which(tau_true == x)),
                         simplify = FALSE) |>
  do.call(what = rbind)

results <- bind_rows(overall_results, subset_results)


# Output Results ----------------------------------------------------------
message(paste0("SLURM_ARRAY_TASK_ID is : ", index))
if(!dir.exists(here::here('Results','CART','temp'))) dir.create(here::here('Results','CART','temp'))
write_csv(results, here::here('Results','CART','temp', paste0(sprintf("%05d", index), ".csv")))
