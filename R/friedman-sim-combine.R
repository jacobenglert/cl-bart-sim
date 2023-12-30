#!/usr/bin/env Rscript

library(tidyverse)

today <- toupper(format(Sys.Date(), "%d%b%Y"))

results <- list.files(here::here('Results','friedman','temp'), full.names = TRUE) |>
lapply(read.csv) |>
bind_rows() |>
write_csv(here::here('Results','friedman', paste0(today, '.csv')))
