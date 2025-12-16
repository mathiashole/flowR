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

# -----------------------------------------------------------
# auxiliary functions

extract_description <- function(attr, keys) { # Extract description from GFF attributes
    for (k in keys) {
        m <- sub(paste0(".*", k, "=([^;]*).*"), "\\1", attr) # Try to match each key
        if (m != attr) return(m) # Return the matched description
    }
    NA # return NA if no match found
}

normalize_description <- function(desc, norm) {
    if (is.na(desc)) return(NA)

    if (isTRUE(norm$decode_url)) {
        desc <- gsub("%2C", ",", desc) # Decode URL-encoded characters
        desc <- gsub("%20", " ", desc) # Decode URL-encoded characters
    }
    if (isTRUE(norm$remove_parentheses)) {
        desc <- gsub("\\(.*?\\)", "", desc) # Remove text within parentheses
    }
    if (isTRUE(norm$trim)) {
        desc <- trimws(desc) # Trim leading and trailing whitespace
    }
    if (isTRUE(norm$collapse_spaces)) {
        desc <- gsub("\\s+", " ", desc) # Collapse multiple spaces
    }
    if (isTRUE(norm$to_lower)) {
        desc <- tolower(desc) # Convert to lowercase
    }
    if (isTRUE(norm$capitalize) && nchar(desc) > 0) {
        desc <- paste0(toupper(substr(desc, 1, 1)), substr(desc, 2, nchar(desc))) # Capitalize first letter
    }
    desc # return normalized description
}

# -----------------------------------------------------------
# Biological function: contextual 5' and 3' gene orientation

get_oriented_context <- function(anchor_row, genes, n = 3) {

same_chr <- genes %>% # Select genes on the same chromosome as the anchor gene
    filter(chr == anchor_row$chr) %>%
    arrange(start) # Arrange by start position

before <- same_chr %>% filter(end < anchor_row$start) # Genes before the anchor gene
after  <- same_chr %>% filter(start > anchor_row$end) # Genes after the anchor gene

}
