#!/usr/bin/env Rscript
# Extract counts from Seurat RDS file and convert to long-format Parquet

library(SeuratObject)
library(Matrix)
library(dplyr)
library(data.table)
library(arrow)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to parse arguments
parse_args <- function(args) {
  parsed <- list(
    format = NULL,
    rds = NULL,
    sample_name = NULL,
    sample_column_name = NULL,
    sample_id = NULL,
    output = NULL
  )
  
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--format") {
      parsed$format <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--rds") {
      parsed$rds <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--sample-name") {
      parsed$sample_name <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--sample-column-name") {
      parsed$sample_column_name <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--sample-id") {
      parsed$sample_id <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output") {
      parsed$output <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  return(parsed)
}

# Clean barcode suffix (remove "-1", "-2", etc. for 10X Genomics format)
clean_barcode_suffix <- function(barcode) {
  if (grepl("^-", barcode)) {
    return(barcode)
  }
  
  # Check if barcode has format with single dash followed by digits
  parts <- strsplit(barcode, "-")[[1]]
  if (length(parts) == 2) {
    prefix <- parts[1]
    suffix <- parts[2]
    
    # Check if suffix is purely numeric and prefix contains only ATCGN
    if (grepl("^[0-9]+$", suffix) && grepl("^[ATCGN]+$", toupper(prefix))) {
      return(prefix)
    }
  }
  
  return(barcode)
}

# Convert matrix to long format data.table (built directly as data.table for efficiency)
matrix_to_long_format <- function(counts_matrix, cleaned_cell_names, gene_names, count_column_name = "Count") {
  if (is(counts_matrix, "sparseMatrix")) {
    counts_summary <- summary(counts_matrix)
    # Build data.table directly from sparse matrix summary
    dt <- data.table(
      CellId = cleaned_cell_names[counts_summary$j],
      GeneId = gene_names[counts_summary$i]
    )
    dt[, (count_column_name) := counts_summary$x]
    # Filter out zeros (shouldn't be in sparse matrix summary, but just in case)
    dt <- dt[get(count_column_name) > 0]
  } else {
    # Dense matrix - convert to long format
    counts_dense <- as.matrix(counts_matrix)
    # Build data.table directly using CJ (cross join) which is faster than expand.grid
    dt <- CJ(CellId = cleaned_cell_names, GeneId = gene_names, sorted = FALSE)
    dt[, (count_column_name) := as.vector(counts_dense)]
    # Filter out zeros
    dt <- dt[get(count_column_name) > 0]
  }
  return(dt)
}

# Add sample_id column to data.table
add_sample_id_column <- function(dt, sample_id, count_column_name = "Count") {
  if (!is.null(sample_id)) {
    # Use data.table syntax for efficient column operations
    dt[, SampleId := sample_id]
    # Reorder columns
    if (count_column_name == "Count") {
      setcolorder(dt, c("SampleId", "CellId", "GeneId", "Count"))
    } else {
      setcolorder(dt, c("SampleId", "CellId", "GeneId", count_column_name))
    }
  }
  return(dt)
}

# Normalize counts matrix (scale to 10,000 per cell)
normalize_counts <- function(counts, seurat_obj) {
  start_time <- Sys.time()
  
  # Get normalized data (slot = "data" contains log-normalized data)
  # If not available, normalize manually
  normalized <- GetAssayData(seurat_obj, slot = "data")
  
  # If data slot is empty or same as counts, normalize manually
  if (length(normalized) == 0 || identical(normalized, counts)) {
    # Manual normalization: scale to 10,000 per cell
    # Optimized for sparse matrices - avoid creating dense intermediate matrices
    if (is(counts, "sparseMatrix")) {
      cell_totals <- colSums(counts)
      # Calculate scale factors
      scale_factors <- 10000 / cell_totals
      scale_factors[is.infinite(scale_factors) | is.nan(scale_factors)] <- 0
      # Efficient column scaling for sparse matrices using diagonal multiplication
      # Create diagonal matrix and multiply (preserves sparsity)
      diag_matrix <- Diagonal(x = scale_factors)
      normalized <- counts %*% diag_matrix
    } else {
      # Dense matrix - use vectorized operations
      cell_totals <- colSums(counts)
      normalized <- sweep(counts, 2, cell_totals, "/") * 10000
    }
  }
  
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(paste("Normalization completed in", round(elapsed_time, 2), "seconds\n"), file = stderr())
  
  return(normalized)
}

# Handle duplicate CellId-GeneId combinations (optimized using data.table)
handle_duplicate_combinations <- function(dt, count_column_name = "Count") {
  start_time <- Sys.time()
  
  # Check for duplicates efficiently
  n_before <- nrow(dt)
  
  # Determine grouping columns
  has_sample_id <- "SampleId" %in% colnames(dt)
  
  if (has_sample_id) {
    n_unique <- uniqueN(dt, by = c("SampleId", "CellId", "GeneId"))
    if (n_before > n_unique) {
      cat(paste("Warning: Found", n_before - n_unique, "duplicate CellId-GeneId combinations. Summing counts.\n"), file = stderr())
      # Use data.table aggregation (much faster than aggregate)
      dt <- dt[, .(Count = sum(get(count_column_name))), by = .(SampleId, CellId, GeneId)]
      setnames(dt, "Count", count_column_name)
      setcolorder(dt, c("SampleId", "CellId", "GeneId", count_column_name))
    }
  } else {
    n_unique <- uniqueN(dt, by = c("CellId", "GeneId"))
    if (n_before > n_unique) {
      cat(paste("Warning: Found", n_before - n_unique, "duplicate CellId-GeneId combinations. Summing counts.\n"), file = stderr())
      # Use data.table aggregation (much faster than aggregate)
      dt <- dt[, .(Count = sum(get(count_column_name))), by = .(CellId, GeneId)]
      setnames(dt, "Count", count_column_name)
      setcolorder(dt, c("CellId", "GeneId", count_column_name))
    }
  }
  
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(paste("Duplicate handling completed in", round(elapsed_time, 2), "seconds\n"), file = stderr())
  
  return(dt)
}

args_parsed <- parse_args(args)

# Validate required arguments
if (is.null(args_parsed$rds)) {
  stop("Error: --rds is required")
}

if (is.null(args_parsed$output)) {
  stop("Error: --output is required")
}

tryCatch({
  # Read the Seurat object
  cat(paste("Reading Seurat RDS file:", args_parsed$rds, "\n"), file = stderr())
  seurat_obj <- readRDS(args_parsed$rds)
  
  # Check if it's a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input file is not a Seurat object")
  }
  
  # Filter by sample if sample_name is provided
  if (!is.null(args_parsed$sample_name)) {
    # Require sample_column_name when filtering by sample
    if (is.null(args_parsed$sample_column_name)) {
      stop("Error: --sample-column-name is required when --sample-name is provided")
    }
    
    cat(paste("Filtering for sample:", args_parsed$sample_name, "\n"), file = stderr())
    cat(paste("Using sample column name:", args_parsed$sample_column_name, "\n"), file = stderr())
    
    metadata <- seurat_obj@meta.data
    
    # Validate that the column exists
    if (!(args_parsed$sample_column_name %in% colnames(metadata))) {
      stop(paste("Provided sample column '", args_parsed$sample_column_name, "' not found in metadata. Available columns: ", paste(colnames(metadata), collapse = ", "), sep = ""))
    }
    
    # Filter Seurat object for the target sample
    seurat_obj <- seurat_obj[, seurat_obj@meta.data[[args_parsed$sample_column_name]] == args_parsed$sample_name]
    
    if (ncol(seurat_obj) == 0) {
      stop(paste("No cells found for sample '", args_parsed$sample_name, "'.", sep = ""))
    }
    
    cat(paste("Filtered to", ncol(seurat_obj), "cells for sample '", args_parsed$sample_name, "'\n", sep = ""), file = stderr())
  }
  
  cat(paste("Seurat object shape:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "genes\n"), file = stderr())
  
  # Get counts matrix
  counts <- GetAssayData(seurat_obj, slot = "counts")
  
  # Get gene and cell identifiers
  gene_names <- rownames(counts)
  cell_names <- colnames(counts)
  
  # Clean cell barcodes if needed
  cleaned_cell_names <- sapply(cell_names, clean_barcode_suffix)
  if (!identical(cleaned_cell_names, cell_names)) {
    cat("Cleaning barcode suffixes...\n", file = stderr())
    cat(paste("Example:", cell_names[1], "->", cleaned_cell_names[1], "\n"), file = stderr())
  }
  
  # Convert to long format for raw counts (already returns data.table)
  cat(paste("Processing", ncol(counts), "cells ×", nrow(counts), "genes...\n"), file = stderr())
  dt_raw <- matrix_to_long_format(counts, cleaned_cell_names, gene_names, "Count")
  dt_raw <- add_sample_id_column(dt_raw, args_parsed$sample_id, "Count")
  
  cat(paste("Total non-zero entries:", nrow(dt_raw), "\n"), file = stderr())
  
  # Handle duplicate CellId-GeneId combinations
  dt_raw <- handle_duplicate_combinations(dt_raw, "Count")
  
  # Ensure output path ends with .parquet
  output_path <- args_parsed$output
  if (!grepl("\\.parquet$", output_path)) {
    # Remove .gz or .csv extension if present, then add .parquet
    output_path <- sub("\\.(csv|gz)$", "", output_path)
    output_path <- paste0(output_path, ".parquet")
  }
  
  # Write raw counts using arrow to Parquet format
  cat(paste("Writing raw count matrix to", output_path, "\n"), file = stderr())
  write_parquet(dt_raw, output_path, compression = "zstd")
  
  # Normalize counts
  cat("Normalizing counts...\n", file = stderr())
  normalized <- normalize_counts(counts, seurat_obj)
  
  # Convert normalized to long format (already returns data.table)
  dt_norm <- matrix_to_long_format(normalized, cleaned_cell_names, gene_names, "NormalizedCount")
  dt_norm <- add_sample_id_column(dt_norm, args_parsed$sample_id, "NormalizedCount")
  
  # Handle duplicate CellId-GeneId combinations in normalized data
  dt_norm <- handle_duplicate_combinations(dt_norm, "NormalizedCount")
  
  # Write normalized counts to Parquet format
  # Generate normalized output path with .parquet extension
  if (grepl("\\.parquet$", output_path)) {
    normalized_output <- sub("\\.parquet$", "_normalized.parquet", output_path)
  } else {
    normalized_output <- paste0(output_path, "_normalized.parquet")
  }
  cat(paste("Writing normalized count matrix to", normalized_output, "\n"), file = stderr())
  write_parquet(dt_norm, normalized_output, compression = "zstd")
  
  cat("Done!\n", file = stderr())
  
}, error = function(e) {
  cat(paste("Error:", conditionMessage(e), "\n"), file = stderr())
  quit(status = 1)
})

