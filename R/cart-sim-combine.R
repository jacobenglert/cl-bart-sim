#!/usr/bin/env Rscript

library(tidyverse)

today <- toupper(format(Sys.Date(), "%d%b%Y"))

results <- list.files(here::here('Results','CART','temp'), full.names = TRUE) |>
lapply(read.csv) |>
bind_rows() |>
write_csv(here::here('Results','CART', paste0(today, '.csv')))
