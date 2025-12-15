#!/usr/bin/env Rscript

# -----------------------------------------------------------
# Script to load configuration from a YAML file and set up the R environment

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(purrr)
  library(ggplot2)
}) # Load necessary libraries
