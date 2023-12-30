# Program Name: CART_summary.R
# Author: Jacob Englert
# Date: 24JUL2023
# Description: Analyze the results of the CART simulation study.

# Load Packages -----------------------------------------------------------
library(tidyverse)

# Load Results ------------------------------------------------------------
sim_date <- "01AUG2023"

# Read in simulation results
results0 <- read_csv(paste0('Results/CART/', sim_date, '.csv'),
                     show_col_types = FALSE, 
                     col_types = cols(waic = col_double(), 
                                      W1 = col_double(), W2 = col_double(), W3 = col_double(),
                                      W4 = col_double(), W5 = col_double(), W6 = col_double(),
                                      W7 = col_double(), W8 = col_double(), W9 = col_double(),
                                      W10 = col_double(), W11 = col_double(), W12 = col_double(),
                                      W13 = col_double(), W14 = col_double(), W15 = col_double(),
                                      W16 = col_double(), W17 = col_double(), W18 = col_double(),
                                      W19 = col_double(), W20 = col_double()))

# Adjust the overall results for CL-BART-avg
clbart_avg_overall <- results0 |>
  filter(model == 'CL-BART-avg') |>
  filter(subset != 'overall') |>
  summarise(coverage = sum(n * coverage) / sum(n),
            bias = sum(n * bias) / sum(n),
            rpehe = sqrt(sum(n * rpehe^2) / sum(n)),
            est = sum(n * est) / sum(n),
            width = sum(n * width) / sum(n),
            .by = c(seed, p, num_trees, corstr, alpha_rho, beta_rho)) |>
  mutate(model = 'CL-BART-avg',
         subset = 'overall')

# Replace corresponding rows
results <- results0 |>
  filter(!(model == 'CL-BART-avg' & subset == 'overall')) |>
  bind_rows(clbart_avg_overall)

# Create Monte Carlo 95% Intervals ----------------------------------------

summaries <- results |>
  pivot_longer(cols = n:W20, names_to = 'Metric', values_to = 'Value') |>
  summarise(mean = mean(Value, na.rm = TRUE),
            l95 = quantile(Value, 0.025, na.rm = TRUE),
            u95 = quantile(Value, 0.975, na.rm = TRUE),
            .by = c(model, subset, num_trees, corstr, p, 
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
  filter(model == 'CL-BART' | (model == 'clogit' & num_trees == 1) | model == 'CL-BART-avg') |>
  mutate(num_trees = ifelse(model == 'clogit', NA, num_trees),
         subset = factor(subset, levels = c('overall', 'W1 = 0 & W2 = 0', 'W1 = 0 & W2 = 1', 'W1 = 1 & W3 = 0', 'W1 = 1 & W3 = 1'),
                         labels = c('Overall','A','B','C','D')),
         model = ifelse(model == 'clogit', 'Oracle CLR', model)) |>
  mutate(model = factor(model, levels = c('Oracle CLR','CL-BART','CL-BART-avg'))) |>
  filter(Metric %in% c('bias','rpehe','coverage','width')) |>
  rowwise() |>
  mutate(Value = paste0(fmt_dec(mean), ' (', 
                        fmt_dec(l95), ', ', 
                        fmt_dec(u95), ')')) |>
  pivot_wider(id_cols = c(p, corstr, subset, alpha_rho, beta_rho, model, num_trees),
              names_from = Metric, values_from = Value) |>
  arrange(p, corstr, subset, alpha_rho, beta_rho, model, num_trees)
write_csv(t1, paste0('Tables/CART/', sim_date, '/t1.csv'))



# Figure (option 1, with clogit)
for(c in unique(summaries$corstr)){
  for(a in unique(summaries$alpha_rho)){
    f1 <- summaries |>
      filter(corstr == c & alpha_rho == a) |>
      filter(model %in% c('CL-BART','clogit','CL-BART-avg')) |>
      filter(subset == 'overall') |>
      filter(Metric %in% c('bias','coverage','rpehe','width')) |>
      ggplot(aes(x = num_trees, y = mean, color = model)) +
      geom_pointrange(aes(ymin = l95, ymax = u95), size = .3, position = position_dodge(.5)) +
      geom_line(position = position_dodge(.5)) +
      geom_hline(aes(yintercept = ref), lty = 2, alpha = 0.4) +
      scale_x_continuous(breaks = unique(results$num_trees)) +
      facet_grid(Metric~paste0('# Predictors: ', p), scales = 'free') + 
      theme_bw() +
      theme(legend.position = 'top') +
      labs(title = 'CART Simulation Results - Overall',
           subtitle = paste0('Predictors have a ', c, ' correlation structure. Regularization Priors = ', a),
           x = 'Number of Trees in Ensemble',
           y = 'Estimate w/ 95% Monte Carlo Confidence Interval',
           color = 'Model')
    ggsave(paste0('Figures/CART/', sim_date, '/f1_', substr(c, 1, 4), '_', a, '.png'), 
           plot = f1, width = 7, height = 6)
  }
}


# summaries |>
#   filter(corstr == 'AR-1') |>
#   filter(model == 'CL-BART' | (model == 'clogit' & num_trees == 1) | model == 'CL-BART-avg') |>
#   filter(subset == 'overall') |>
#   filter(Metric %in% c('bias','coverage','rpehe')) |>
#   pivot_wider(id_cols = c(p, corstr, alpha_rho, beta_rho, num_trees, Metric),
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




# Subgroup Results --------------------------------------------------------

# Figure (option 1, with clogit)
for(c in unique(summaries$corstr)){
  for(a in unique(summaries$alpha_rho)){
    f2 <- summaries |>
      filter(corstr == c & alpha_rho == a) |>
      filter(model %in% c('CL-BART','clogit','CL-BART-avg')) |>
      filter(subset != 'overall') |>
      filter(Metric %in% c('est','bias','coverage','rpehe','width')) |>
      ggplot(aes(x = factor(num_trees), y = mean, shape = model, color = subset)) +
      geom_pointrange(aes(ymin = l95, ymax = u95), size = .3, position = position_dodge(.5)) +
      #geom_line(position = position_dodge(.5)) +
      geom_hline(aes(yintercept = ref), lty = 2, alpha = 0.4) +
      #scale_x_continuous(breaks = unique(results$num_trees)) +
      facet_grid(Metric~paste0('# Predictors: ', p), scales = 'free') + 
      theme_bw() +
      theme(legend.position = 'top') +
      labs(title = 'CART Simulation Results - Subgroups',
           subtitle = paste0('Predictors have a ', c, ' correlation structure. Regularization Priors = ', a),
           x = 'Number of Trees in Ensemble',
           y = 'Estimate w/ 95% Monte Carlo Confidence Interval',
           shape = 'Model',
           color = 'Subgroup')
    ggsave(paste0('Figures/CART/', sim_date, '/f2_', substr(c, 1, 4), '_', a, '.png'), 
           plot = f2, width = 9, height = 5)
  }
}


# In-Depth Analysis of CL-BART Model Fits ---------------------------------

# Analyze the relative performance change when adding more trees
clbart_rel <- results |>
  filter(subset == 'overall') |>
  filter(model == 'CL-BART') |>
  pivot_longer(cols = bias:waic, names_to = 'Metric', values_to = 'Value') |>
  mutate(ref_val = min(abs(Value)), 
         rel_val = abs(Value) / ref_val,
         .by = c(Metric, p, corstr, alpha_rho, beta_rho, seed)) |>
  summarise(mean_rel_val = mean(rel_val), 
            .by = c(Metric, p, corstr, alpha_rho, beta_rho, num_trees)) |>
  pivot_wider(id_cols = c(corstr, p, alpha_rho, beta_rho, num_trees), 
              names_from = Metric, values_from = mean_rel_val)
clbart_rel

# Visualize dip in WAIC
clbart_rel |>
  filter(corstr == 'AR-1') |>
  pivot_longer(cols = bias:waic, names_to = 'Metric', values_to = 'Value') |>
  filter(Metric != 'width') |>
  ggplot(aes(x = factor(num_trees), y = Value, group = interaction(alpha_rho, beta_rho, corstr),
             color = interaction(alpha_rho, beta_rho), lty = corstr)) +
  geom_point() +
  geom_line() +
  facet_grid(Metric~paste0('# Predictors: ', p), scales = 'free') + 
  theme_bw() +
  theme(legend.position = 'top')


# WAIC plot
cart_waic <- clbart_rel |>
  filter(p == 10 & corstr == 'AR-1' & alpha_rho == 0.95) |>
  pivot_longer(cols = bias:waic, names_to = 'Metric', values_to = 'Value') |>
  filter(Metric == 'waic') |>
  ggplot(aes(x = factor(num_trees), y = Value, group = interaction(alpha_rho, beta_rho, corstr), lty = corstr)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  theme_bw() +
#  facet_wrap(~Metric, scales = 'free') +
  theme(legend.position = 'top') +
  labs(x = '# Trees',
       y = 'Average Relative WAIC')
ggsave(cart_waic, paste0('Figures/CART/', sim_date, '/waic.png'), width = 4, height = 4)

# clbart_avg_rel <- results |>
#   filter(subset != 'overall') |>
#   filter(model == 'CL-BART-avg') |>
#   pivot_longer(cols = coverage:width, names_to = 'Metric', values_to = 'Value') |>
#   mutate(ref_val = min(abs(Value)), 
#          rel_val = abs(Value) / ref_val,
#          .by = c(subset, Metric, p, corstr, alpha_rho, beta_rho, seed)) |>
#   summarise(mean_rel_val = mean(rel_val), 
#             .by = c(subset, Metric, p, corstr, alpha_rho, beta_rho, num_trees)) |>
#   pivot_wider(id_cols = c(corstr, p, alpha_rho, beta_rho, subset, num_trees), 
#               names_from = Metric, values_from = mean_rel_val) |>
#   arrange(corstr, p, alpha_rho, beta_rho, subset, num_trees)
# clbart_avg_rel
  

# Distribution of Predictions ---------------------------------------------
results |>
  filter(model %in% c('clogit','CL-BART')) |>
  ggplot(aes(x = est, fill = model)) +
  geom_histogram(alpha = 0.4, position = 'identity') +
  geom_vline(xintercept = log(c(.8,1.1,1.3,1.5))) +
  facet_grid(num_trees~alpha_rho, scales = 'free')


# Variable Importance
results |>
  filter(model == 'CL-BART') |>
  filter(subset == 'overall') |>
  filter(corstr == 'AR-1') |>
  pivot_longer(cols = W1:W20, names_to = 'Moderator', values_to = 'Value') |>
  mutate(Value = coalesce(Value, 0),
         Moderator = factor(Moderator, levels = paste0('W',1:20))) |>
  filter(!(p == 10 & Moderator %in% paste0('W',11:20))) |>
  summarise(mean = mean(Value),
            .by = c(Moderator, corstr, p, num_trees, alpha_rho, beta_rho)) |>
  ggplot(aes(x = Moderator, y = mean, color = factor(num_trees), group = num_trees)) +
  geom_point() +
  geom_line(lty = 5) +
  facet_grid(interaction(alpha_rho, beta_rho)~p, scales = 'free') +
  theme_bw()


# variable importance plot
cart_imp <- results |>
  filter(model == 'CL-BART') |>
  filter(subset == 'overall') |>
  filter(corstr == 'AR-1' & alpha_rho == 0.95 & p == 10) |>
  pivot_longer(cols = W1:W20, names_to = 'Moderator', values_to = 'Value') |>
  mutate(Value = coalesce(Value, 0),
         Moderator = factor(Moderator, levels = paste0('W',1:20))) |>
  filter(!(p == 10 & Moderator %in% paste0('W',11:20))) |>
  summarise(mean = mean(Value),
            .by = c(Moderator, corstr, p, num_trees, alpha_rho, beta_rho)) |>
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
ggsave(cart_imp, paste0('Figures/CART/', sim_date, '/CART_imp.png'), width = 6.5, height = 5)


library(patchwork)
cart_waic + cart_imp + plot_annotation(tag_levels = 'A')
ggsave(paste0('Figures/CART/', sim_date, '/CART_waic_imp.png'), width = 6.5, height = 3.5)



t2 <- results |>
  filter(model == 'CL-BART') |>
  filter(subset == 'overall') |>
  pivot_longer(cols = W1:W20, names_to = 'Moderator', values_to = 'Value') |>
  mutate(Value = coalesce(Value, 0),
         Moderator = factor(Moderator, levels = paste0('W',1:20))) |>
  filter(!(p == 10 & Moderator %in% paste0('W',11:20))) |>
  summarise(mean = mean(Value),
            .by = c(Moderator, subset, corstr, p, num_trees, alpha_rho, beta_rho)) |>
  pivot_wider(id_cols = subset:beta_rho, names_from = Moderator, values_from = mean)



# Helpful overall summaries -----------------------------------------------

refs <- summaries |>
  filter(p == 10) |>
  filter(subset == 'overall') |>
  filter(model == 'clogit' & num_trees == 1) |>
  select(corstr, p, alpha_rho, beta_rho, Metric, refmid = mean, reflow = l95, refhigh = u95)
  
summaries |>
  filter(p == 10) |>
  filter(subset == 'overall') |>
  filter(model %in% c('CL-BART','CL-BART-avg')) |>
  filter(Metric %in% c('bias','coverage','rpehe')) |>
  left_join(refs) |>
  ggplot(aes(x = factor(num_trees), y = mean, color = model, 
             shape = interaction(alpha_rho, beta_rho), 
             lty = interaction(alpha_rho, beta_rho))) +
  geom_pointrange(aes(ymin = l95, ymax = u95), size = .3, position = position_dodge(.5)) +
  geom_line(position = position_dodge(.5)) +
  geom_hline(aes(yintercept = refmid), lty = 2, alpha = 0.5) +
  geom_hline(aes(yintercept = reflow), lty = 2, alpha = 0.5) +
  geom_hline(aes(yintercept = refhigh), lty = 2, alpha = 0.5) +
  facet_grid(Metric~corstr, scales = 'free') +
  theme_bw() +
#  scale_x_continuous(breaks = unique(results$num_trees)) +
  labs(title = 'Simulation Results - Overall',
       x = 'Number of Trees in Ensemble',
       y = 'Estimate w/ 95% Monte Carlo Confidence Interval',
       color = 'Regularization')
