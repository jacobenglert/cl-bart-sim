# Program Name: friedman_example.R
# Author: Jacob Englert
# Date: 24JUL2023
# Description: A CL-BART example where the true log odds-ratio is determined
# by a modified version of the Friedman function.

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(survival)
library(clbart2)

# Set Simulation Parameters -----------------------------------------------
seed    <- 2187
set.seed(seed)                              # Random Seed
n       <- 10000                            # Initial Population Size
start   <- as.Date("2020-01-01")            # Start of follow-up period
end     <- as.Date("2022-12-31")            # End of follow-up period
dates   <- seq.Date(start, end, by = 'day') # 3-year follow-up period
t       <- length(dates)                    # Length of follow-up period
p       <- 10                               # Number of effect moderators

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

# Visualize probabilities over time
data |>
  filter(ID %in% sample(data$ID, 10)) |>
  ggplot(aes(x = date, y = p, color = factor(ID))) +
  geom_line(alpha = 0.4, show.legend = FALSE) +
  geom_smooth(se = FALSE, show.legend = FALSE) +
  theme_bw() +
  labs(x = 'Date', 
       y = 'Probability of ER Admission',
       title = 'Probability of ER Admission During the Study Period')

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

rm(data)

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
summary(fit_clogit)

# # Fit the LASSO model -----------------------------------------------------
# w_cc_2 <- as.data.frame(model.matrix(~ .^2, data = w_cc))[,-1] * z_cc
# fit_lasso <- clogitL1(x = w_cc_2, y = y_cc, strata = strata_cc)
# fit_lasso_cv <- cv.clogitL1(fit_lasso, numFolds = 5)
# lambda_idx <- which(fit_lasso_cv$lambda == fit_lasso_cv$minCV_lambda)
# beta_lasso_cv <- fit_lasso_cv$beta[lambda_idx,]
# tau_lasso <- as.numeric(as.matrix(w_cc_2) %*% beta_lasso_cv)[c(T,F,F,F)]


# Fit CL-BART model -------------------------------------------------------
fit_clbart <- clbart(w = w_cc, y = y_cc, z = z_cc, stratum = strata_cc, 
                     num_trees = 1, seed = 1, iter = 6000, thin = 5, warmup = 3000, 
                     alpha_sigma = 1, beta_sigma = 1)

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
tau_clbart_post <- lapply(fit_clbart$forests, \(f) predict_forest(f, w_unique)) |>
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
                      num_trees = fit_clbart$num_trees,
                      p = p,
                      alpha_rho = fit_clbart$alpha_rho,
                      beta_rho = fit_clbart$beta_rho,
                      alpha_sigma = 1,
                      beta_sigma = 1,
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




# Model Fit Exploration ---------------------------------------------------


# Visualize CL-BART Model Diagnostics
fit_clbart[c('tree_acc_rate','time','logLik')] |>
  bind_cols() |>
  mutate(Iteration = row_number()) |>
  pivot_longer(cols = -Iteration, names_to = 'Parameter', values_to = 'Value') |>
  ggplot(aes(x = Iteration, y = Value)) +
  geom_line() +
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw()

# Visualize Tree Complexity
fit_clbart[c('avg_tree_depth','avg_num_nodes','avg_num_leaves')] |>
  bind_cols() |>
  mutate(Iteration = row_number()) |>
  pivot_longer(cols = -Iteration, names_to = 'Parameter', values_to = 'Value') |>
  ggplot(aes(x = Iteration, y = Value, color = Parameter)) +
  geom_line() +
  theme_bw()

# Visualize Tree Splitting Variables
fit_clbart$split_props |>
  data.frame() |> 
  mutate(Iteration = row_number()) |>
  pivot_longer(cols = -Iteration, names_to = 'Covariate', values_to = 'Value') |>
  summarise(post_mean = mean(Value, na.rm = TRUE),
            post_l95 = quantile(Value, 0.025, na.rm = TRUE),
            post_u95 = quantile(Value, 0.975, na.rm = TRUE),
            .by = Covariate) |>
  ggplot(aes(x = Covariate, y = post_mean)) +
  geom_col(fill = 'white', color = 'black') +
  geom_errorbar(aes(ymin = post_l95, ymax = post_u95), width = 0.5) +
  theme_bw()

fit_clbart$split_props |>
  data.frame() |>
  mutate(Iteration = row_number()) |>
  pivot_longer(cols = -Iteration, names_to = 'Covariate', values_to = 'Value') |>
  ggplot(aes(x = Iteration, y = Covariate, fill = Value)) +
  geom_tile() +
  theme_bw()

# Calculate CLORs (Conditional Log Odds Ratios)
clor <- tau_clbart_post |>
  as.data.frame() |>
  mutate(Iteration = row_number())

# Visualize 10 CLORs
clor |> 
  pivot_longer(cols = -Iteration, names_to = 'ID', values_to = 'CLOR') |>
  group_by(ID) |>
  group_split() |>
  sample(20) |>
  bind_rows() |>
  mutate(color = mean(CLOR), .by = ID) |>
  ggplot(aes(x = Iteration, y = exp(CLOR), color = color)) +
  geom_line(show.legend = FALSE) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_wrap(~ID) +
  theme_bw()

# Visualize ACLOR (Average Conditional Log Odds Ratio)
clor |> 
  pivot_longer(cols = -Iteration, names_to = 'ID', values_to = 'CLOR') |>
  summarise(ACLOR = mean(CLOR), .by = Iteration) |>
  ggplot(aes(x = Iteration, y = exp(ACLOR))) +
  geom_line() +
  theme_bw()

# Calculate Marginal Partial ACLORs
library(parallel) # to speed up computation
pACLOR <- mclapply(colnames(mod_df), \(x) comp_pd(model = fit_clbart, model_data = w_unique, vars = x, resolution = 10), mc.cores = 4)
pACLOR_clean <- pACLOR |>
  bind_rows() |>
  pivot_longer(cols = all_of(colnames(mod_df)), names_to = 'Variable', values_to = 'Value') |>
  filter(!is.na(Value))

# Visualize Marginal pACLORs
pACLOR_clean |>
  summarise(pACLOR = mean(pd),
            lower = quantile(pd, .025),
            upper = quantile(pd, .975),
            .by = c(Variable, Value)) |>
  ggplot(aes(x = Value, y = pACLOR)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  # geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~Variable, scales = 'free') +
  theme_bw()


# Partial Dependence Function
comp_pd <- function(model = NULL, model_data = NULL, new_data = NULL, vars = NULL, resolution = 5){
  
  if(is.null(new_data)){
    pd_groups <- sapply(vars, \(v) seq(min(model_data[v]), max(model_data[v]), length.out = resolution), simplify = FALSE) |>
      expand.grid()
    
    #pd_groups <- seq(min(model_data[vars]), max(model_data[vars]), length.out = 5)
    new_data <- model_data[!(names(model_data) %in% vars)] |>
      dplyr::group_by(across(everything())) |>
      dplyr::summarise(n = n(), .groups = 'drop') |>
      merge(pd_groups)
  }
  
  pd <- model$forests |>
    lapply(\(f){
      cbind(new_data, preds = predict_forest(f, new_data = new_data[!(names(new_data) == 'n')])) |>
        dplyr::group_by(across(all_of(vars))) |>
        dplyr::summarise(pd = weighted.mean(preds, n), .groups = 'drop')
    }) |>
    do.call(what = rbind)
  
  pd$Iteration <- rep(1:(nrow(pd) / nrow(pd_groups)), each = nrow(pd_groups))
  
  return(pd)
  
}



