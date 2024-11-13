gc()
rm(list = ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)

MINCELLS  <- 50
MINGENES  <- 5

dataset.name <- "GSE127465"

raw.data.dir <- "/media/hieunguyen/GSHD_HN01/raw_data/UKK_Lung/raw_data"
savedir <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/UKK_Lung/preprocessed_raw_data", dataset.name)

all.files <- Sys.glob(file.path(raw.data.dir, dataset.name, dataset.name, "*.tsv"))
print(sprintf("Number of samples in this dataset %s: %s", dataset.name, length(all.files)))

path_to_expr <- all.files[[1]]
input.data <- read.table(path_to_expr, sep = "\t", header = TRUE) %>%
  column_to_rownames("barcode") %>% t()

sample.id <- str_split(basename(path_to_expr), "[.]")[[1]][[1]]

PROJECT <- sprintf("%s_%s", dataset.name, sample.id)
s.obj <- CreateSeuratObject(counts = input.data, 
                            min.cells = MINCELLS, 
                            min.features = MINGENES, 
                            project = PROJECT)
saveRDS(object = s.obj, file.path(savedir, sprintf("%s.rds", sample.id)))