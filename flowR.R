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
fill_file <- NULL

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
    } else if (args[i] == "--fill_file" || args[i] == "-f") {
        fill_file <- args[i + 1]
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

if (anchor_row$strand == "+") {
    upstream   <- tail(before, n) # Upstream genes (5' direction)
    downstream <- head(after, n) # Downstream genes (3' direction)
} else {
    upstream   <- head(after, n) # Upstream genes (5' direction)
    downstream <- tail(before, n) # Downstream genes (3' direction)
}

upstream <- upstream %>% # Assign relative positions to upstream genes
    arrange(start) %>% # Arrange by start position
    mutate(relative_position = -rev(seq_len(nrow(.)))) # Negative positions for upstream

downstream <- downstream %>% # Assign relative positions to downstream genes
    arrange(start) %>% # Arrange by start position
    mutate(relative_position = seq_len(nrow(.))) # Positive positions for downstream

center <- anchor_row %>%
    mutate(relative_position = 0) # Center gene with position 0

bind_rows(upstream, center, downstream) # Combine upstream, center, and downstream genes
}

# -----------------------------------------------------------
# Process re-annotated genes
if (!is.null(fill_file)) {
    override_raw <- read.delim(fill_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(override_raw)[1:3] <- c("chr", "v_tmp1", "v_tmp2")

    override <- override_raw %>%
        mutate(
        strand = ifelse(v_tmp1 <= v_tmp2, "+", "-"), # Detecte strand based on coordinates before normalization
        start = pmin(as.numeric(v_tmp1), as.numeric(v_tmp2)), # force start to be the smaller value is important for filtering
        end   = pmax(as.numeric(v_tmp1), as.numeric(v_tmp2)), # force end to be the larger value is important for filtering
        protein = if(ncol(override_raw) >= 4) as.character(override_raw[[4]]) else gene_query # unique protein name if not provided
    ) %>%
    select(chr, start, end, strand, protein) # Select relevant columns
    
    message("-> Normalized fill file: ", nrow(override), " gene charged.")
}

# -----------------------------------------------------------
# Main processing genomes
all_contexts <- list()
all_genes    <- list()

for (g in seq_along(gff_files)) {
    message("-> Processing: ", gff_files[g])
    gff <- read.delim(gff_files[g], header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
    cfg <- read_yaml(yaml_files[g])

    genes <- gff %>%
        filter(V3 == cfg$features$gene_type) %>% # Filter for gene features
        transmute(
        chr    = V1,
        start  = as.numeric(V4),
        end    = as.numeric(V5),
        strand = V7,
        attr   = V9
        ) %>% # Extract relevant columns
        arrange(chr, start) # Arrange by chromosome and start position
        
    all_genes[[g]] <- genes # Store all genes for potential further analysis

# Extract protein descriptions
    genes$protein <- vapply(
        genes$attr,
        extract_description,
        character(1),
        keys = cfg$features$attribute_keys
    )
# Normalize protein descriptions
    genes$protein <- vapply(
        genes$protein,
        normalize_description,
        character(1),
        norm = cfg$normalization
    )
# Apply replacements if specified
    if (!is.null(cfg$replacements)) { # if replacements are defined in the config
        genes$protein <- dplyr::recode(
        genes$protein,
        !!!cfg$replacements # Unquote-splice the replacements list
        )
    }
# Get contexts around anchor genes
    if (!is.null(fill_file)) { # Use fill file if provided
        anchors <- override # Filter for anchor gene from fill file
    } else {
    anchors <- genes %>%
        filter(grepl(gene_query, attr, ignore.case = TRUE)) # Filter for anchor gene
    }

    if (nrow(anchors) == 0) next # Skip if no anchor gene found
# Get contexts for each anchor gene
    ctx <- map_dfr(
        split(anchors, seq_len(nrow(anchors))),
        ~ get_oriented_context(.x, genes, n = window)
    ) # Combine contexts into a single data frame

# Get contexts for each anchor gene
    ctx$n_anchors <- nrow(anchors)
    ctx$genome <- cfg$genome
    ctx <- ctx %>% filter(relative_position != 0) # Remove anchor gene

    all_contexts[[g]] <- ctx
}

context_data <- bind_rows(all_contexts) # Combine all contexts into a single data frame
stopifnot("relative_position" %in% colnames(context_data)) # Ensure relative_position column exists

all_genes_df <- bind_rows(all_genes) # Combine all genes into a single data frame
universe <- all_genes_df$protein # Define universe of proteins
universe <- universe[!is.na(universe)] # Remove NA values

# -----------------------------------------------------------
# Fisher's Exact Test for enrichment analysis

# -----------------------------------------------------------
# Select top N frequencies proteins

top_proteins <- context_data %>% # Identify top N most frequent proteins
    count(protein, sort = TRUE) %>%
    slice_head(n = top_n) %>%
    pull(protein)

plot_data <- context_data %>% # prepare data for plotting
    filter(protein %in% top_proteins) %>%
    count(protein, relative_position)
# continue visualization even if some positions have zero counts
plot_data <- bind_rows(
    plot_data,
    expand.grid(
        protein = top_proteins,
        relative_position = 0,
        n = 0
    )
)
# save data
write.table(
    plot_data,
    file = sub("\\.png$", "_context_counts.tsv", plot_file), # Output file name
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# -----------------------------------------------------------
# Visualization
# draw plot only if plot_file is provided
if (!is.null(plot_file)) {

# n_anchors <- nrow(anchors) # Number of anchor genes found
n_anchors <- sum(unique(context_data$n_anchors)) # Total number of anchor genes across genomes

gg <- ggplot(plot_data, aes(x = relative_position, y = n, color = protein, group = protein)) + # Create ggplot object
    geom_line(linewidth = 1) + # Add lines for each protein
    geom_point(size = 3, alpha = 0.8) + # Add points for each protein
    annotate("text", x = Inf, y = Inf, label = paste0("n = ", n_anchors), hjust = 1.1, vjust = 1.5, size = 4) + # Annotate number of genes
    scale_x_continuous(breaks = -window:window, labels = c(paste0("5′ -", window:1), gene_query, paste0("3′ +", 1:window))) + # Customize x-axis labels
    labs(
        x = paste0("5′ <- ", gene_query, " -> 3′"),
        y = "Number of genes",
        color = "Protein",
    ) + # Add labels
    theme_classic() +
    theme(legend.position = "bottom", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) # Customize theme

ggsave(plot_file, gg, width = 12, height = 8, dpi = 800) # Save plot to file
} # debug make plot every time