library(dplyr)
library(ggplot2)
library(vroom)
library(stringr)
library(here)
library(ggpubr)
library(knitr)
library(GenomicRanges)
library(magrittr)
library(tidyr)

qtl_set <- here("pigQTL/pigQTLdb.gff") %>% vroom(delim = "\t", comment = "#", col_names = FALSE) %>% 
  dplyr::rename(seqnames = X1, db = X2, qtl_class = X3, start = X4, end = X5, qtl_details = X9) %>% 
  filter(!is.na(start) | !is.na(end)) %>% 
  mutate(
    qtl_name = qtl_details %>% 
      str_split("Name=") %>% 
      sapply(function(x) x[2]) %>% 
      str_split(";") %>% 
      sapply(function(x) x[1]) %>% 
      str_remove_all("\"") %>% 
      str_trim()
  ) %>% 
  dplyr::select(seqnames, start, end, qtl_class, qtl_name, qtl_details) %>% 
  mutate(seqnames = str_remove(seqnames, "Chr.")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
  sort()

# export for later use
saveRDS(qtl_set,  here("pigQTL/qtl_set.rds"))
