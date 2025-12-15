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

# -----------------------------------------------------------
# Parse command line arguments
for (i in seq_along(args)) {
    if (args[i] == "--gff"  || args[i] == "-g") {
        gff_files <- c(gff_files, args[i + 1])
    } else if (args[i] == "--yaml" || args[i] == "-y") {
        yaml_files <- c(yaml_files, args[i + 1])
    } else if (args[i] == "--gene" || args[i] == "-g") {
        gene_query <- args[i + 1]
    } else if (args[i] == "--window" || args[i] == "-w") {
    window <- as.integer(args[i + 1])
    } else if (args[i] == "--top" || args[i] == "-t") {
    top_n <- as.integer(args[i + 1])
    }
}