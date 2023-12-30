# Program Name: friedman_summary.R
# Author: Jacob Englert
# Date: 25JUL2023
# Description: Analyze the results of the Friedman simulation study.

# Load Packages -----------------------------------------------------------
library(tidyverse)

# Load Results ------------------------------------------------------------
sim_date <- "31JUL2023"

# Read in simulation results
results <- read_csv(paste0('Results/friedman/', sim_date, '.csv'),
                    show_col_types = FALSE, 
                    col_types = cols(waic = col_double(), 
                                     W1 = col_double(), W2 = col_double(), W3 = col_double(),
                                     W4 = col_double(), W5 = col_double(), W6 = col_double(),
                                     W7 = col_double(), W8 = col_double(), W9 = col_double(),
                                     W10 = col_double(), W11 = col_double(), W12 = col_double(),
                                     W13 = col_double(), W14 = col_double(), W15 = col_double(),
                                     W16 = col_double(), W17 = col_double(), W18 = col_double(),
                                     W19 = col_double(), W20 = col_double()))

# 0.95, 50 tree, seed 18 & 69, p = 10

# Create Monte Carlo 95% Intervals ----------------------------------------

summaries <- results |>
  pivot_longer(cols = n:W20, names_to = 'Metric', values_to = 'Value') |>
  summarise(mean = mean(Value, na.rm = TRUE),
            l95 = quantile(Value, 0.025, na.rm = TRUE),
            u95 = quantile(Value, 0.975, na.rm = TRUE),
            .by = c(model, num_trees, p, 
                    alpha_rho, beta_rho, Metric)) |>
  mutate(ref = case_when(Metric == 'bias' ~ 0,
                         Metric == 'coverage' ~ .95,
                         Metric == 'rpehe' ~ 0)) |>
  filter(!is.na(mean))

# Overall Results ---------------------------------------------------------

fmt_dec <- function(num, ndec = 2){
  new_num <- round(num, ndec)
  if(new_num < 0) sprintf(paste0("%.", ndec, "f"), new_num)
  else sprintf(paste0("%.", ndec, "f"), abs(new_num))
}

# First Table
t1 <- summaries |>
  filter(model == 'CL-BART' | (model == 'clogit' & num_trees == 5)) |>
  mutate(num_trees = ifelse(model == 'clogit', NA, num_trees),
         model = ifelse(model == 'clogit', 'Oracle CLR', model)) |>
  mutate(model = factor(model, levels = c('Oracle CLR','CL-BART'))) |>
  filter(Metric %in% c('bias','rpehe','coverage','width')) |>
  rowwise() |>
  mutate(Value = paste0(fmt_dec(mean), ' (', 
                        fmt_dec(l95), ', ', 
                        fmt_dec(u95), ')')) |>
  pivot_wider(id_cols = c(p, alpha_rho, beta_rho, model, num_trees),
              names_from = Metric, values_from = Value) |>
  arrange(p, alpha_rho, beta_rho, model, num_trees)
write_csv(t1, paste0('Tables/friedman/', sim_date, '/t1.csv'))

# Table
t1_overall <- summaries |>
  filter(model == 'CL-BART' | (model == 'clogit' & num_trees == 5)) |>
  mutate(num_trees = ifelse(model == 'clogit', NA, num_trees)) |>
  filter(Metric %in% c('bias','rpehe','coverage','width')) |>
  mutate(Value = paste0(round(mean, 2), ' (', round(l95, 2), ', ', round(u95, 2), ')')) |>
  pivot_wider(id_cols = c(p, alpha_rho, beta_rho, num_trees, model),
              names_from = Metric, values_from = Value) |>
  arrange(p, alpha_rho, beta_rho, desc(model), num_trees)
write_csv(t1_overall, 'Tables/friedman/t1.csv')

# Figure (option 1, with clogit)
f1 <- summaries |>
  filter(model %in% c('CL-BART','clogit')) |> #,'CL-BART-avg')) |>
  filter(Metric %in% c('bias','coverage','rpehe')) |>
  ggplot(aes(x = num_trees, y = mean, color = model, lty = interaction(alpha_rho, beta_rho))) +
  geom_pointrange(aes(ymin = l95, ymax = u95), size = .3, position = position_dodge(.5)) +
  geom_line(position = position_dodge(.5)) +
  geom_hline(aes(yintercept = ref), lty = 2, alpha = 0.4) +
  scale_x_continuous(breaks = unique(results$num_trees)) +
  facet_grid(Metric~paste0('# Predictors: ', p), scales = 'free') + 
  theme_bw() +
  theme(legend.position = 'top') +
  labs(title = 'Friedman Simulation Results',
       x = 'Number of Trees in Ensemble',
       y = 'Estimate w/ 95% Monte Carlo Confidence Interval',
       color = 'Model')
f1
ggsave(paste0('Figures/friedman/', sim_date, '/f1.png'), width = 6, height = 5)

# summaries |>
#   filter(model == 'CL-BART' | (model == 'clogit' & num_trees == 5)) |>
#   filter(Metric %in% c('bias','coverage','rpehe')) |>
#   pivot_wider(id_cols = c(p, alpha_rho, beta_rho, num_trees, Metric),
#               names_from = model, values_from = c(mean, l95, u95)) |>
#   ggplot(aes(x = num_trees, y = `mean_CL-BART`, color = interaction(alpha_rho, beta_rho))) +
#   geom_pointrange(aes(ymin = `l95_CL-BART`, ymax = `u95_CL-BART`), size = .3, position = position_dodge(.5)) +
#   geom_line(position = position_dodge(.5)) +
#   geom_hline(aes(yintercept = mean_clogit), lty = 2, alpha = 0.5) +
#   geom_hline(aes(yintercept = l95_clogit), lty = 2, alpha = 0.5) +
#   geom_hline(aes(yintercept = u95_clogit), lty = 2, alpha = 0.5) +
#   facet_grid(Metric~paste0('# Predictors: ', p), scales = 'free') +
#   theme_bw() +
#   scale_x_continuous(breaks = unique(results$num_trees)) +
#   labs(title = 'Simulation Results - Overall',
#        subtitle = paste0('Predictors sampled from a set of predictors with an ', unique(results$corstr), ' correlation structure.'),
#        x = 'Number of Trees in Ensemble',
#        y = 'Estimate w/ 95% Monte Carlo Confidence Interval',
#        color = 'Regularization')



# In-Depth Analysis of CL-BART Model Fits ---------------------------------

# Just Summaries
f2 <- summaries |>
  filter(model == 'CL-BART') |>
  filter(Metric %in% c('bias','coverage','rpehe','waic')) |>
  ggplot(aes(x = factor(num_trees), y = mean, lty = interaction(alpha_rho, beta_rho))) +
  geom_pointrange(aes(ymin = l95, ymax = u95), size = .3, position = position_dodge(.5)) +
  #geom_line(position = position_dodge(.5)) +
  #geom_hline(aes(yintercept = ref), lty = 2, alpha = 0.4) +
  #scale_x_continuous(breaks = unique(results$num_trees)) +
  facet_grid(Metric~paste0('# Predictors: ', p), scales = 'free') + 
  theme_bw() +
  theme(legend.position = 'top') +
  labs(title = 'Friedman Simulation Results',
       x = 'Number of Trees in Ensemble',
       y = 'Estimate w/ 95% Monte Carlo Confidence Interval')
f2

# Analyze the relative performance change when adding more trees
clbart_rel <- results |>
  filter(model == 'CL-BART') |>
  pivot_longer(cols = bias:waic, names_to = 'Metric', values_to = 'Value') |>
  mutate(ref_val = min(abs(Value)), 
         rel_val = abs(Value) / ref_val,
         .by = c(Metric, p, alpha_rho, beta_rho, seed)) |>
  summarise(mean_rel_val = mean(rel_val), 
            .by = c(Metric, p, alpha_rho, beta_rho, num_trees)) |>
  pivot_wider(id_cols = c(p, alpha_rho, beta_rho, num_trees), 
              names_from = Metric, values_from = mean_rel_val)
clbart_rel

clbart_rel |>
  pivot_longer(cols = bias:waic, names_to = 'Metric', values_to = 'Value') |>
  filter(Metric != 'width') |>
  ggplot(aes(x = factor(num_trees), y = Value, group = interaction(alpha_rho, beta_rho),
             color = interaction(alpha_rho, beta_rho))) +
  geom_point() +
  geom_line() +
  facet_grid(Metric~paste0('# Predictors: ', p), scales = 'free') + 
  theme_bw() +
  theme(legend.position = 'top')

# WAIC plot
fried_waic <- clbart_rel |>
  filter(p == 10 & alpha_rho == 0.95) |>
  pivot_longer(cols = bias:waic, names_to = 'Metric', values_to = 'Value') |>
  filter(Metric == 'waic') |>
  ggplot(aes(x = factor(num_trees), y = Value, group = interaction(alpha_rho, beta_rho))) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  theme_bw() +
  #  facet_wrap(~Metric, scales = 'free') +
  theme(legend.position = 'top') +
  labs(x = '# Trees',
       y = 'Average Relative WAIC')
ggsave(paste0('Figures/friedman/', sim_date, '/friedman_waic.png'), width = 4, height = 4)

# Variable Importance
results |>
  filter(model == 'CL-BART') |>
  pivot_longer(cols = W1:W20, names_to = 'Moderator', values_to = 'Value') |>
  mutate(Value = coalesce(Value, 0),
         Moderator = factor(Moderator, levels = paste0('W',1:20))) |>
  filter(!(p == 10 & Moderator %in% paste0('W',11:20))) |>
  summarise(mean = mean(Value),
            .by = c(Moderator, p, num_trees, alpha_rho, beta_rho)) |>
  ggplot(aes(x = Moderator, y = mean, color = factor(num_trees), group = num_trees)) +
  geom_point() +
  geom_line(lty = 5) +
  facet_grid(interaction(alpha_rho, beta_rho)~p, scales = 'free') +
  theme_bw()

# variable importance plot
fried_imp <- results |>
  filter(model == 'CL-BART') |>
  filter(alpha_rho == 0.95 & p == 10) |>
  pivot_longer(cols = W1:W20, names_to = 'Moderator', values_to = 'Value') |>
  mutate(Value = coalesce(Value, 0),
         Moderator = factor(Moderator, levels = paste0('W',1:20))) |>
  filter(!(p == 10 & Moderator %in% paste0('W',11:20))) |>
  summarise(mean = mean(Value),
            .by = c(Moderator, num_trees, alpha_rho, beta_rho)) |>
  ggplot(aes(x = Moderator, y = mean, color = factor(num_trees), group = num_trees)) +
  geom_point() +
  geom_line(lty = 5) +
  theme_bw() +
  theme(legend.position = c(.95,.95),
        legend.box.background = element_rect(color = "black", linewidth = 0.5),
        legend.justification = c(1,1)) +
  labs(x = 'Effect Moderator',
       y = 'Average Variable Split Proportion',
       color = '# Trees')
ggsave(paste0('Figures/friedman/', sim_date, '/friedman_imp.png'), width = 6.5, height = 5)
  
library(patchwork)
fried_waic + fried_imp + plot_annotation(tag_levels = 'A')
ggsave(paste0('Figures/friedman/', sim_date, '/friedman_waic_imp.png'), width = 6.5, height = 3.5)
