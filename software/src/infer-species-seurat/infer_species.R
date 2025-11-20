#!/usr/bin/env Rscript
# Infer species from Seurat RDS file based on gene identifiers

library(SeuratObject)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to parse arguments
parse_args <- function(args) {
  parsed <- list(
    format = NULL,
    input_file = NULL
  )
  
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--format" || args[i] == "-f") {
      parsed$format <- args[i + 1]
      i <- i + 2
    } else if (!startsWith(args[i], "--")) {
      # First non-flag argument is the input file (positional)
      if (is.null(parsed$input_file)) {
        parsed$input_file <- args[i]
      }
      i <- i + 1
    } else {
      i <- i + 1
    }
  }
  
  return(parsed)
}

# Detect gene identifier type
detect_gene_format <- function(gene_names) {
  # Sample a subset of genes for analysis
  sample_size <- min(1000, length(gene_names))
  sample_genes <- sample(gene_names, sample_size)
  
  # Check for Ensembl IDs
  ensembl_human <- sum(grepl("^ENSG\\d{11}$", sample_genes))
  ensembl_mouse <- sum(grepl("^ENSMUSG\\d{11}$", sample_genes))
  ensembl_rat <- sum(grepl("^ENSRNOG\\d{11}$", sample_genes))
  ensembl_zebrafish <- sum(grepl("^ENSDARG\\d{11}$", sample_genes))
  ensembl_chicken <- sum(grepl("^ENSGALG\\d{11}$", sample_genes))
  ensembl_cow <- sum(grepl("^ENSBTAG\\d{11}$", sample_genes))
  ensembl_pig <- sum(grepl("^ENSSSCG\\d{11}$", sample_genes))
  
  total_ensembl <- ensembl_human + ensembl_mouse + ensembl_rat + ensembl_zebrafish + 
                   ensembl_chicken + ensembl_cow + ensembl_pig
  
  if (total_ensembl > length(sample_genes) * 0.5) {
    # Mostly Ensembl IDs
    if (ensembl_human > ensembl_mouse && ensembl_human > ensembl_rat) {
      return(list(type = "ensembl", species = "homo-sapiens"))
    } else if (ensembl_mouse > ensembl_rat) {
      return(list(type = "ensembl", species = "mus-musculus"))
    } else if (ensembl_rat > 0) {
      return(list(type = "ensembl", species = "rattus-norvegicus"))
    } else if (ensembl_zebrafish > 0) {
      return(list(type = "ensembl", species = "danio-rerio"))
    } else if (ensembl_chicken > 0) {
      return(list(type = "ensembl", species = "gallus-gallus"))
    } else if (ensembl_cow > 0) {
      return(list(type = "ensembl", species = "bos-taurus"))
    } else if (ensembl_pig > 0) {
      return(list(type = "ensembl", species = "sus-scrofa"))
    } else {
      return(list(type = "ensembl", species = "homo-sapiens"))  # Default
    }
  }
  
  # Check for Entrez IDs (numeric)
  numeric_genes <- sum(grepl("^\\d+$", sample_genes))
  if (numeric_genes > length(sample_genes) * 0.5) {
    # Mostly numeric - likely Entrez IDs
    # Try to infer species from ranges (simplified)
    gene_nums <- as.numeric(sample_genes[grepl("^\\d+$", sample_genes)])
    gene_nums <- gene_nums[!is.na(gene_nums)]
    
    if (length(gene_nums) > 0) {
      median_id <- median(gene_nums)
      if (median_id < 100000) {
        return(list(type = "entrez", species = "homo-sapiens"))
      } else if (median_id < 200000) {
        return(list(type = "entrez", species = "mus-musculus"))
      } else if (median_id < 300000) {
        return(list(type = "entrez", species = "rattus-norvegicus"))
      } else {
        return(list(type = "entrez", species = "homo-sapiens"))  # Default
      }
    }
    return(list(type = "entrez", species = "homo-sapiens"))
  }
  
  # Check gene symbol patterns
  # Human: mostly uppercase
  uppercase_genes <- sum(grepl("^[A-Z]{2,}$", sample_genes))
  # Mouse: mixed case starting with uppercase
  mixed_case_genes <- sum(grepl("^[A-Z][a-z]+", sample_genes))
  
  if (uppercase_genes > length(sample_genes) * 0.3) {
    # Likely human gene symbols
    return(list(type = "symbol", species = "homo-sapiens"))
  } else if (mixed_case_genes > length(sample_genes) * 0.3) {
    # Likely mouse gene symbols
    return(list(type = "symbol", species = "mus-musculus"))
  }
  
  # Default to human if uncertain
  return(list(type = "symbol", species = "homo-sapiens"))
}

args_parsed <- parse_args(args)

# Validate required arguments
if (is.null(args_parsed$input_file)) {
  stop("Error: input file is required as positional argument")
}

# Output filenames
output_file <- "species.txt"
format_output_file <- "gene_format.txt"

tryCatch({
  # Read the Seurat object
  cat(paste("Reading Seurat RDS file:", args_parsed$input_file, "\n"), file = stderr())
  seurat_obj <- readRDS(args_parsed$input_file)
  
  # Check if it's a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input file is not a Seurat object")
  }
  
  # Get gene names from counts matrix
  counts <- GetAssayData(seurat_obj, slot = "counts")
  gene_names <- rownames(counts)
  
  cat(paste("Analyzing", length(gene_names), "genes...\n"), file = stderr())
  
  # Detect gene format and infer species
  result <- detect_gene_format(gene_names)
  
  # Map identifier type to format name
  format_mapping <- list(
    "symbol" = "gene symbol",
    "ensembl" = "Ensembl Id",
    "entrez" = "Entrez Id"
  )
  
  gene_format <- format_mapping[[result$type]]
  if (is.null(gene_format)) {
    gene_format <- "unknown"
  }
  
  # Calculate data size: rows = genes, columns = cells
  num_genes <- nrow(counts)
  num_cells <- ncol(counts)
  
  # Write species output
  cat(paste("Inferred species:", result$species, "\n"), file = stderr())
  cat(paste("Gene format:", gene_format, "\n"), file = stderr())
  cat(paste("Data size:", num_genes, "genes (rows) ×", num_cells, "cells (columns)\n"), file = stderr())
  
  # Write without trailing newline
  cat(result$species, file = output_file, sep = "")
  cat(gene_format, file = format_output_file, sep = "")
  
  # Write data size output (rows columns)
  data_size_output_file <- "data_size.txt"
  cat(paste(num_genes, num_cells), file = data_size_output_file, sep = "")
  
  cat("Species inference completed successfully\n", file = stderr())
  
}, error = function(e) {
  cat(paste("Error:", conditionMessage(e), "\n"), file = stderr())
  quit(status = 1)
})

