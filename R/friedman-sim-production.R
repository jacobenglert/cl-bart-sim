# Program Name: friedman_sim.R
# Author: Jacob Englert
# Date: 24JUL2023
# Description: A CL-BART simulation study where the true log odds-ratio 
# is determined by a modified version of the Friedman function.

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(survival)
# library(clogitL1)
library(clbart2)
library(parallel)

# Get Parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c('index')
index <- as.numeric(args['index'])
params <- read_csv(paste0('Params/friedman_params.csv'), show_col_types = F)[,-1]
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
w_cor <- diag(p)

# Define function to calculate Friedman function log odds-ratio
get_complex_logOR <- function(x){
  
  t1 <- (0.5 * sin(pi*x[,1]*x[,2]) + 2*(x[,3]-.5)^2 + .5*x[,4] + .5*x[,5])
  
  return(log(t1))
}

# Simulate shared exposure time series with sinusoidal time trend
z <- rnorm(t, sin(seq(0, 2 * pi * 3, length.out = t)))

# Create exposure data frame
expos_df <- data.frame(z = z, date = dates) |>
  mutate(year = year(date), month = month(date), dow = wday(date))

# Simulate continuous exposure moderators
mod_df <- replicate(p, runif(n)) |> as.data.frame() |> round(2)
colnames(mod_df) <- paste0('W', 1:p)

# Follow each individual throughout time
data <- cross_join(data.frame(ID = 1:n, mod_df), expos_df)

# Compute the heterogeneous exposure effects on the log-odds scale
w_imp <- sample(colnames(mod_df), 5)
data$tau <- get_complex_logOR(data[w_imp])

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

rm(data, mod_df, expos_df)

# Extract individual elements
w_cc <- select(data_cc, starts_with('W'))
y_cc <- data_cc$y
z_cc <- data_cc$z
tau_cc <- data_cc$tau
strata_cc <- data_cc$ID


# Fit "Oracle" Conditional Logistic Regression Model ----------------------
w_cc_imp <- data_cc[w_imp]
w_cc_true <- cbind(sin(pi * w_cc_imp[,1] * w_cc_imp[,2]),
                   (w_cc_imp[,3] - .5)^2,
                   w_cc_imp[,4],
                   w_cc_imp[,5])

fit_clogit <- clogit(y_cc ~ w_cc_true : z_cc + strata(strata_cc))


# # Fit the LASSO model -----------------------------------------------------
# w_cc_2 <- as.data.frame(model.matrix(~ .^2, data = w_cc))[,-1] * z_cc
# fit_lasso <- clogitL1(x = w_cc_2, y = y_cc, strata = strata_cc)
# fit_lasso_cv <- cv.clogitL1(fit_lasso, numFolds = 5)
# lambda_idx <- which(fit_lasso_cv$lambda == fit_lasso_cv$minCV_lambda)
# beta_lasso_cv <- fit_lasso_cv$beta[lambda_idx,]
# tau_lasso <- as.numeric(as.matrix(w_cc_2) %*% beta_lasso_cv)[c(T,F,F,F)]


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

# cl <- makeCluster(5)
# 
# clusterExport(cl, varlist = c("w_cc", "y_cc", "z_cc", "strata_cc",
#                               "num_trees", "iter", "thin", "warmup",
#                               "sigma2_beta", "sigma2_beta_update_freq",
#                               "beta_acc_prob", "n_min", "moves", "move_probs",
#                               "alpha_rho", "beta_rho", "alpha_sigma", "beta_sigma"))
# 
# 
# fit_clbarts <- clusterApply(cl, 1:5, \(s){ 
#   clbart2::clbart(seed = s, 
#                   w = w_cc, y = y_cc, z = z_cc, stratum = strata_cc, 
#                   num_trees = num_trees,
#                   iter = iter, thin = thin, warmup = warmup, 
#                   sigma2_beta = sigma2_beta, 
#                   sigma2_beta_update_freq = sigma2_beta_update_freq,
#                   beta_acc_prob = beta_acc_prob, 
#                   n_min = n_min, 
#                   moves = moves, move_probs = move_probs, 
#                   alpha_rho = alpha_rho, beta_rho = beta_rho, 
#                   alpha_sigma = alpha_sigma, beta_sigma = beta_sigma)})
# 
# stopCluster(cl)

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


# Calculate simulation statistics
est_clogit      <- mean(tau_clogit)
bias_clogit     <- mean(tau_clogit - tau_true)
rpehe_clogit    <- sqrt(mean((tau_clogit - tau_true)^2))
coverage_clogit <- mean(between(tau_true, tau_clogit_lower, tau_clogit_upper))
width_clogit    <- mean(tau_clogit_upper - tau_clogit_lower)

# est_lasso <- mean(tau_lasso)
# bias_lasso <- mean(tau_lasso - tau_true)
# rpehe_lasso <- sqrt(mean((tau_lasso - tau_true)^2))

est_clbart      <- mean(tau_clbart)
bias_clbart     <- mean(tau_clbart - tau_true)
rpehe_clbart    <- sqrt(mean((tau_clbart - tau_true)^2))
coverage_clbart <- mean(between(tau_true, tau_clbart_lower, tau_clbart_upper))
width_clbart    <- mean(tau_clbart_upper - tau_clbart_lower)

tau_bar_mean  <- mean(rowMeans(tau_clbart_post))
tau_bar_lower <- quantile(rowMeans(tau_clbart_post), 0.025)
tau_bar_upper <- quantile(rowMeans(tau_clbart_post), 0.975)
tau_bar_bias  <- mean(tau_bar_mean - tau_true)
tau_bar_rpehe <- sqrt(mean((tau_bar_mean - tau_true)^2))
tau_bar_coverage <- mean(between(tau_true, tau_bar_lower, tau_bar_upper))
tau_bar_width <- tau_bar_upper - tau_bar_lower

results <- data.frame(model = c('clogit', 'CL-BART', 'CL-BART-avg'), # 'LASSO'
                      seed = seed,
                      num_trees = num_trees,
                      p = p,
                      alpha_rho = alpha_rho,
                      beta_rho = beta_rho,
                      alpha_sigma = alpha_sigma,
                      beta_sigma = beta_sigma,
                      n = length(tau_true),
                      est = c(est_clogit, est_clbart, tau_bar_mean), # est_lasso
                      bias = c(bias_clogit, bias_clbart, tau_bar_bias), # bias_lasso
                      rpehe = c(rpehe_clogit, rpehe_clbart, tau_bar_rpehe), # rpehe_lasso
                      coverage = c(coverage_clogit, coverage_clbart, tau_bar_coverage),
                      width = c(width_clogit, width_clbart, tau_bar_width))

results2 <- data.frame(model = 'CL-BART', waic = fit_clbart$WAIC) |>
  bind_cols(t(c(colMeans(fit_clbart$split_props, na.rm = TRUE)[w_imp],
                colMeans(fit_clbart$split_props, na.rm = TRUE)[!(colnames(w_cc) %in% w_imp)])|> 
                setNames(paste0('W', 1:p))))

results <- left_join(results, results2, by = 'model')


# Output Results ----------------------------------------------------------
message(paste0("SLURM_ARRAY_TASK_ID is : ", index))
if(!dir.exists("Results/friedman/temp")) dir.create("Results/friedman/temp")
write_csv(results, paste0("Results/friedman/temp/", sprintf("%05d", index), ".csv"))
