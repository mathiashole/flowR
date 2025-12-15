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
    } else if (args[i] == "--plot" || args[i] == "-p") {
        plot_file <- args[i + 1]
    }
}

# -----------------------------------------------------------
# Validate required arguments
if (length(gff_files) == 0 || length(yaml_files) == 0 || is.null(gene_query)) {
    stop("Usage: Rscript flowR.R --gff genome.gff --yaml genome.yaml --gene DGF-1 [--window 3] [--top 8] [--plot out.png]") # Check for required arguments
}

if (length(gff_files) != length(yaml_files)) {
    stop("Number of GFF files must match number of YAML files") # Ensure matching number of files
}