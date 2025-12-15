#!/usr/bin/env Rscript

# -----------------------------------------------------------
# Script to load configuration from a YAML file and set up the R environment

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(purrr)
  library(ggplot2)
}) # Load necessary libraries

args <- commandArgs(trailingOnly = TRUE)

gff_files  <- character()
yaml_files <- character()
gene_query <- NULL
window     <- 3
top_n      <- 8
plot_file  <- NULL

for (i in seq_along(args)) {
    if (args[i] == "--gff") {
        gff_files <- c(gff_files, args[i + 1])
    }
}