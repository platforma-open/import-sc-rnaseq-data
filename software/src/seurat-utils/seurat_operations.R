#!/usr/bin/env Rscript
# Unified script for Seurat RDS operations: extract-counts or infer-species

library(SeuratObject)
library(Matrix)
library(dplyr)
library(data.table)
library(arrow)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description = "Unified script for Seurat RDS operations")
parser$add_argument("operation", choices = c("extract-counts", "infer-species", "check-input"), 
                    help = "Operation to perform: extract-counts, infer-species, or check-input")
parser$add_argument("--format", help = "Gene format (optional)", default = NULL)
parser$add_argument("--rds", help = "Path to input Seurat RDS file (required for extract-counts)")
parser$add_argument("--sample-name", help = "Sample name to filter (optional)", default = NULL)
parser$add_argument("--sample-column-name", help = "Sample column name in metadata (required if --sample-name is provided)", default = NULL)
parser$add_argument("--sample-id", help = "Sample ID to add to output (optional)", default = NULL)
parser$add_argument("--output", help = "Path to output Parquet file (required for extract-counts)")
parser$add_argument("--error-output", help = "Path to error output TSV file (required for check-input)", default = "check_error.tsv")
parser$add_argument("input_file", nargs = "?", help = "Path to input Seurat RDS file (required for infer-species and check-input)")

args <- parser$parse_args()

# Helper function to get assay data with compatibility for both old and new Seurat versions
# Seurat 5.0.0+ uses 'layer' argument, older versions use 'slot'
get_assay_data_compat <- function(seurat_obj, layer_name) {
  tryCatch({
    # Try new API first (Seurat 5.0.0+)
    GetAssayData(seurat_obj, layer = layer_name)
  }, error = function(e) {
    # Fall back to old API (Seurat < 5.0.0)
    GetAssayData(seurat_obj, slot = layer_name)
  })
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
  
  # Get normalized data (layer = "data" contains log-normalized data)
  # If not available, normalize manually
  normalized <- get_assay_data_compat(seurat_obj, "data")
  
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

# Operation: extract-counts
operation_extract_counts <- function() {
  if (is.null(args$rds)) {
    stop("Error: --rds is required for extract-counts operation")
  }
  if (is.null(args$output)) {
    stop("Error: --output is required for extract-counts operation")
  }
  
  tryCatch({
    # Read the Seurat object
    cat(paste("Reading Seurat RDS file:", args$rds, "\n"), file = stderr())
    seurat_obj <- readRDS(args$rds)
    
    # Check if it's a Seurat object
    if (!inherits(seurat_obj, "Seurat")) {
      stop("Input file is not a Seurat object")
    }
    
    # Filter by sample if sample_name is provided
    if (!is.null(args$sample_name)) {
      # Require sample_column_name when filtering by sample
      if (is.null(args$sample_column_name)) {
        stop("Error: --sample-column-name is required when --sample-name is provided")
      }
      
      cat(paste("Filtering for sample:", args$sample_name, "\n"), file = stderr())
      cat(paste("Using sample column name:", args$sample_column_name, "\n"), file = stderr())
      
      metadata <- seurat_obj@meta.data
      
      # Validate that the column exists
      if (!(args$sample_column_name %in% colnames(metadata))) {
        stop(paste("Provided sample column '", args$sample_column_name, "' not found in metadata. Available columns: ", paste(colnames(metadata), collapse = ", "), sep = ""))
      }
      
      # Filter Seurat object for the target sample
      seurat_obj <- seurat_obj[, seurat_obj@meta.data[[args$sample_column_name]] == args$sample_name]
      
      if (ncol(seurat_obj) == 0) {
        stop(paste("No cells found for sample '", args$sample_name, "'.", sep = ""))
      }
      
      cat(paste("Filtered to", ncol(seurat_obj), "cells for sample '", args$sample_name, "'\n", sep = ""), file = stderr())
    }
    
    cat(paste("Seurat object shape:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "genes\n"), file = stderr())
    
    # Get counts matrix
    counts <- get_assay_data_compat(seurat_obj, "counts")
    
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
    dt_raw <- add_sample_id_column(dt_raw, args$sample_id, "Count")
    
    cat(paste("Total non-zero entries:", nrow(dt_raw), "\n"), file = stderr())
    
    # Handle duplicate CellId-GeneId combinations
    dt_raw <- handle_duplicate_combinations(dt_raw, "Count")
    
    # Ensure output path ends with .parquet
    output_path <- args$output
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
    dt_norm <- add_sample_id_column(dt_norm, args$sample_id, "NormalizedCount")
    
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
}

# Operation: infer-species
operation_infer_species <- function() {
  input_file <- args$input_file
  if (is.null(input_file)) {
    stop("Error: input_file is required for infer-species operation")
  }
  
  # Output filenames
  output_file <- "species.txt"
  format_output_file <- "gene_format.txt"
  
  tryCatch({
    # Read the Seurat object
    cat(paste("Reading Seurat RDS file:", input_file, "\n"), file = stderr())
    seurat_obj <- readRDS(input_file)
    
    # Check if it's a Seurat object
    if (!inherits(seurat_obj, "Seurat")) {
      stop("Input file is not a Seurat object")
    }
    
    # Get gene names from counts matrix
    counts <- get_assay_data_compat(seurat_obj, "counts")
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
}

# Operation: check-input
operation_check_input <- function() {
  input_file <- args$input_file
  error_output_file <- args$error_output
  
  if (is.null(input_file)) {
    stop("Error: input_file is required for check-input operation")
  }
  
  # Function to write error to TSV file
  write_error_tsv <- function(index, message, file_path) {
    error_df <- data.frame(Index = as.character(index), Log = message, stringsAsFactors = FALSE)
    write.table(error_df, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Function to estimate RDS file memory requirement
  estimate_memory_needed <- function(file_path) {
    file_size <- file.info(file_path)$size
    estimated_memory <- file_size * 1.5
    return(estimated_memory)
  }
  
  # Function to get available memory using R's built-in functions
  get_available_memory <- function() {
    tryCatch({
      # Get R's maximum vector allocation size (returns MB)
      max_vsize_mb <- mem.maxVSize()
      
      cat(paste("Debug: mem.maxVSize() returned:", max_vsize_mb, "MB\n"), file = stderr())
      
      # Check if max_vsize is NA, Inf, or unreasonably large (unlimited)
      # On many Unix systems it returns Inf (unlimited)
      if (is.na(max_vsize_mb) || is.infinite(max_vsize_mb) || max_vsize_mb > 1e9) {
        cat("Debug: max_vsize is NA, Inf, or too large (unlimited)\n", file = stderr())
        return(NA)
      }
      
      # Check if max_vsize is unreasonably small (less than 100MB)
      if (max_vsize_mb < 100) {
        cat(paste("Debug: max_vsize is too small (", max_vsize_mb, "MB), likely invalid\n"), file = stderr())
        return(NA)
      }
      
      # Convert to bytes
      max_vsize_bytes <- max_vsize_mb * 1024^2
      
      # Get current memory usage by R
      gc_info <- gc(reset = FALSE, full = TRUE)
      # Sum of used memory (Ncells and Vcells in MB)
      used_mb <- sum(gc_info[, "(Mb)"])
      used_bytes <- used_mb * 1024^2
      
      cat(paste("Debug: Current R memory usage:", round(used_mb, 2), "MB\n"), file = stderr())
      
      # Available memory is max_vsize minus what's currently used
      available_bytes <- max_vsize_bytes - used_bytes
      
      cat(paste("Debug: Calculated available memory:", round(available_bytes / 1024^2, 2), "MB\n"), file = stderr())
      
      # Return available bytes (conservative estimate)
      # Return NA if the calculation doesn't make sense
      if (available_bytes < 0) {
        cat("Debug: Calculated available memory is negative\n", file = stderr())
        return(NA)
      }
      
      return(available_bytes)
    }, error = function(e) {
      cat(paste("Debug: Error in get_available_memory():", conditionMessage(e), "\n"), file = stderr())
      # If we can't determine available memory, return NA
      return(NA)
    })
  }
  
  tryCatch({
    # Check if file exists
    if (!file.exists(input_file)) {
      write_error_tsv(1, paste("File not found:", input_file), error_output_file)
      cat(paste("Error: File not found:", input_file, "\n"), file = stderr())
      quit(status = 0)  # Exit successfully, error was written
    }
    
    cat(paste("Checking Seurat RDS file:", input_file, "\n"), file = stderr())
    
    # Get file size
    file_size <- file.info(input_file)$size
    file_size_mb <- round(file_size / (1024^2), 2)
    cat(paste("File size:", file_size_mb, "MB\n"), file = stderr())
    
    # Estimate memory needed
    estimated_memory <- estimate_memory_needed(input_file)
    estimated_memory_mb <- round(estimated_memory / (1024^2), 2)
    cat(paste("Estimated memory needed:", estimated_memory_mb, "MB\n"), file = stderr())
    
    # Check available memory
    available_memory <- get_available_memory()
    if (!is.na(available_memory)) {
      available_memory_mb <- round(available_memory / (1024^2), 2)
      cat(paste("Available memory:", available_memory_mb, "MB\n"), file = stderr())
      
      # Check if we have enough memory (with 20% safety margin)
      if (available_memory < estimated_memory * 1.2) {
        error_msg <- sprintf("Possible memory issue: file requires ~%.2f MB, but only %.2f MB appears to be available. Loading may fail. Consider using a machine with more RAM if you encounter errors.", 
                           estimated_memory_mb, available_memory_mb)
        write_error_tsv(6, error_msg, error_output_file)
        cat(paste("Warning:", error_msg, "\n"), file = stderr())
        quit(status = 0)  # Exit successfully, warning was written
      }
    } else {
      cat("Note: Could not determine available memory, proceeding with load attempt\n", file = stderr())
    }
    
    # Try to read the RDS file
    cat("Attempting to load RDS file...\n", file = stderr())
    seurat_obj <- readRDS(input_file)
    
    # Check if it's a Seurat object
    if (!inherits(seurat_obj, "Seurat")) {
      write_error_tsv(7, "Input file is not a valid Seurat object. Please provide an RDS file containing a Seurat object.", error_output_file)
      cat("Error: Input file is not a valid Seurat object\n", file = stderr())
      quit(status = 0)  # Exit successfully, error was written
    }
    
    # Basic validation - check for counts matrix
    counts <- tryCatch({
      get_assay_data_compat(seurat_obj, "counts")
    }, error = function(e) {
      error_msg <- sprintf("Failed to extract counts matrix from Seurat object: %s", conditionMessage(e))
      write_error_tsv(8, error_msg, error_output_file)
      cat(paste("Error:", error_msg, "\n"), file = stderr())
      quit(status = 0)  # Exit successfully, error was written
    })
    
    if (is.null(counts) || length(counts) == 0) {
      write_error_tsv(8, "Seurat object does not contain a valid counts matrix. The 'counts' slot is empty or missing.", error_output_file)
      cat("Error: Seurat object does not contain a valid counts matrix\n", file = stderr())
      quit(status = 0)  # Exit successfully, error was written
    }
    
    # Check for gene identifiers
    gene_names <- rownames(counts)
    if (is.null(gene_names) || length(gene_names) == 0) {
      write_error_tsv(9, "No gene identifiers found in the counts matrix. Rows must have names (gene IDs).", error_output_file)
      cat("Error: No gene identifiers found in the counts matrix\n", file = stderr())
      quit(status = 0)  # Exit successfully, error was written
    }
    
    # Check for cell identifiers
    cell_names <- colnames(counts)
    if (is.null(cell_names) || length(cell_names) == 0) {
      write_error_tsv(10, "No cell identifiers found in the counts matrix. Columns must have names (cell barcodes).", error_output_file)
      cat("Error: No cell identifiers found in the counts matrix\n", file = stderr())
      quit(status = 0)  # Exit successfully, error was written
    }
    
    # If we got here, everything is OK
    cat(paste("Seurat object validated successfully:", nrow(counts), "genes ×", ncol(counts), "cells\n"), file = stderr())
    cat("File format check passed\n", file = stderr())
    
    # Write empty error file (success)
    write.table(data.frame(Index = character(0), Log = character(0)), 
                file = error_output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
  }, error = function(e) {
    # Catch any unexpected errors
    error_msg <- paste("Error reading or validating Seurat file:", conditionMessage(e))
    write_error_tsv(11, error_msg, error_output_file)
    cat(paste(error_msg, "\n"), file = stderr())
    quit(status = 0)  # Exit successfully, error was written
  })
}

# Main execution
tryCatch({
  if (args$operation == "extract-counts") {
    operation_extract_counts()
  } else if (args$operation == "infer-species") {
    operation_infer_species()
  } else if (args$operation == "check-input") {
    operation_check_input()
  } else {
    stop(paste("Unknown operation:", args$operation))
  }
}, error = function(e) {
  cat(paste("Error:", conditionMessage(e), "\n"), file = stderr())
  quit(status = 1)
})

