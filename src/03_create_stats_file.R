library(here)
library(tidyverse)
library(data.table)
fib <- readRDS(here("results/02/Fibroblast/Fibroblast_stats.rds"))[["n_matching"]] %>% unlist() %>% data.frame() %>% rownames_to_column("rn")
hg2 <- readRDS(here("results/02/HG002/HG002_stats.rds"))[["n_matching"]] %>% unlist() %>% data.frame() %>% rownames_to_column("rn")
hg4 <- readRDS(here("results/02/HG004/HG004_stats.rds"))[["n_matching"]] %>% unlist() %>% data.frame() %>% rownames_to_column("rn")
hg5 <- readRDS(here("results/02/HG005/HG005_stats.rds"))[["n_matching"]] %>% unlist() %>% data.frame() %>% rownames_to_column("rn")

stats <- left_join(fib, hg2, by = "rn")
stats <- left_join(stats, hg4, by = "rn")
stats <- left_join(stats, hg5, by = "rn")
colnames(stats) <- c("rn", "Fibroblast", "HG002", "HG004", "HG005")

fwrite(stats, here("results/03/n_matching_stats.csv"), bom = TRUE)
