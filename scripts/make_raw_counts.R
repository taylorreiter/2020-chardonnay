library(dplyr)
library(purrr)
library(readr)
library(tidyr)

files <- unlist(snakemake@input)
#files <- list.files("outputs/htseq", pattern = ".txt$", full.names = T)

raw_counts <- files %>%
  purrr::set_names() %>% 
  map_dfr(function(x) read_tsv(x, col_names = c("gene", "count")), .id = "sample") %>%
  mutate(sample = gsub("outputs\\/htseq\\/", "", sample)) %>%
  mutate(sample = gsub("_UMI.*", "", sample)) %>%
  pivot_wider(id_cols = gene, names_from = sample, values_from = count)

write_tsv(raw_counts, unlist(snakemake@output))
