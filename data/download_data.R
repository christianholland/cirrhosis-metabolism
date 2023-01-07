# Copyright (c) [2022] [Christian H. Holland]
# cholland2408@gmail.com

# This script downloads data required for analyses

library(here)
library(purrr)


# Download data from Zenodo (doi: 10.5281/zenodo.7242764) -----------------
# set timeout to 10 mins i.e. 600 seconds
options(timeout = 600)

base_url <- "https://zenodo.org/record/7242764/files"

# download data folder and extract certain files
download.file(file.path(base_url, "data.zip?download=1"), here("data/data.zip"))
data_files <- c(
  "data/mouse-chronic-ccl4/meta_data.rds",
  "data/mouse-chronic-ccl4/count_matrix.rds",
  "data/meta-mouse-vs-human/contrast_annotation.rds",
  "data/annotation/gene_id_annotation.rds",
  "data/human-hoang-nafld/meta_data.rds"
)

unzip(here("data/data.zip"), files = data_files, exdir = here())

# download output folder and extract certain files
download.file(
  file.path(base_url, "output.zip?download=1"),
  here("data/output.zip")
)
output_files <- c(
  "output/mouse-chronic-ccl4/limma_result.rds",
  "output/mouse-chronic-ccl4/stem_result.rds",
  "output/mouse-chronic-ccl4/normalized_expression.rds",
  "output/meta-mouse-vs-human/limma_result.rds",
  "output/meta-mouse-vs-human/meta_data.rds",
  "output/human-hoang-nafld/normalized_expression.rds"
)

unzip(here("data/output.zip"), files = output_files, exdir = here("data"))

# move extracted files to
walk(
  list.files(here("data/output"), pattern = "*.rds", recursive = TRUE),
  function(file) {
    file.rename(here("data", "output", file), here("data", file))
  }
)

# remove folders that were created by unzipping
unlink(here("data/output/"), recursive = TRUE)


# Download pericentral genes ----------------------------------------------

pericentral_genes <- read_csv("https://raw.githubusercontent.com/saezlab/LiverPeriportalization/master/paper/material/consensus_pericentral_geneset.csv")
saveRDS(pericentral_genes, here("data/annotation/pericentral_genes.rds"))
