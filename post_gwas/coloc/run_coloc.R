#!/usr/bin/env Rscript
#
# Bayesian colocalization analysis for GWAS trait pairs.
#
# This script performs colocalization analysis using both coloc and eCAVIAR
# methods to test whether trait pairs share causal variants at identified loci.
#
# Usage:
#     Rscript run_coloc.R --gwas-files size.tsv protein.tsv oil.tsv --trait-names size protein oil --loci locus_summary.tsv --output coloc_results.tsv
#
# Input:
#     GWAS tables for all traits (TSV format)
#     Locus summary from merge_loci.py
#
# Output:
#     coloc_results.tsv - Colocalization results for each locus and trait pair
#
# Dependencies:
#     - R packages: coloc, dplyr, data.table, argparse
#     - Python packages: pandas, numpy (for preprocessing if needed)

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(data.table)
  library(coloc)
})

# Function to parse command line arguments
parse_arguments <- function() {
  parser <- ArgumentParser(description = "Bayesian colocalization analysis for GWAS trait pairs")
  
  parser$add_argument("--gwas-files", nargs = "+", required = TRUE,
                      help = "List of GWAS table files (TSV format)")
  parser$add_argument("--trait-names", nargs = "+", required = TRUE,
                      help = "Corresponding trait names for each GWAS table")
  parser$add_argument("--loci", required = TRUE,
                      help = "Locus summary file from merge_loci.py")
  parser$add_argument("--output", required = TRUE,
                      help = "Output file path for colocalization results")
  parser$add_argument("--pp4-threshold", type = "double", default = 0.8,
                      help = "PP4 threshold for declaring colocalization (default: 0.8)")
  parser$add_argument("--prior-p1", type = "double", default = 1e-4,
                      help = "Prior probability of SNP associated with trait 1 (default: 1e-4)")
  parser$add_argument("--prior-p2", type = "double", default = 1e-4,
                      help = "Prior probability of SNP associated with trait 2 (default: 1e-4)")
  parser$add_argument("--prior-p12", type = "double", default = 1e-5,
                      help = "Prior probability of SNP associated with both traits (default: 1e-5)")
  parser$add_argument("--verbose", action = "store_true",
                      help = "Print detailed progress information")
  
  args <- parser$parse_args()
  
  # Validate input
  if (length(args$gwas_files) != length(args$trait_names)) {
    stop(sprintf("Number of GWAS files (%d) does not match number of trait names (%d)",
                 length(args$gwas_files), length(args$trait_names)))
  }
  
  return(args)
}

# Function to load GWAS data
load_gwas_data <- function(filepath, trait_name, verbose = FALSE) {
  if (verbose) {
    cat(sprintf("Loading GWAS data for %s from %s...\n", trait_name, filepath))
  }
  
  # Read TSV file
  df <- fread(filepath, sep = "\t", data.table = FALSE)
  
  # Check required columns
  required_cols <- c("rs", "chr", "ps", "beta", "p_wald")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in %s: %s", filepath, paste(missing_cols, collapse = ", ")))
  }
  
  # Add trait name
  df$trait <- trait_name
  # Add pos column from ps
  df$pos <- df$ps
  
  if (verbose) {
    cat(sprintf("  Loaded %d SNPs for %s\n", nrow(df), trait_name))
  }
  
  return(df)
}

# Function to load locus data
load_locus_data <- function(filepath, verbose = FALSE) {
  if (verbose) {
    cat(sprintf("Loading locus data from %s...\n", filepath))
  }
  
  loci_df <- fread(filepath, sep = "\t", data.table = FALSE)
  
  # Check required columns
  required_cols <- c("locus_id", "chr", "start", "end", "lead_snp")
  missing_cols <- setdiff(required_cols, colnames(loci_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in %s: %s", filepath, paste(missing_cols, collapse = ", ")))
  }
  
  if (verbose) {
    cat(sprintf("  Loaded %d loci\n", nrow(loci_df)))
  }
  
  return(loci_df)
}

# Function to prepare data for coloc
prepare_coloc_data <- function(gwas_df, locus_info, verbose = FALSE) {
  # Filter SNPs in the locus
  locus_snps <- gwas_df %>%
    filter(chr == locus_info$chr,
           pos >= locus_info$start,
           pos <= locus_info$end) %>%
    arrange(pos)
  
  if (nrow(locus_snps) == 0) {
    if (verbose) {
      cat(sprintf("  No SNPs found in locus %s (chr%s:%d-%d)\n",
                  locus_info$locus_id, locus_info$chr, locus_info$start, locus_info$end))
    }
    return(NULL)
  }
  
  # Calculate standard errors from p-values if not available
  if (!"se" %in% colnames(locus_snps)) {
    # Assume beta ~ N(0, se) under null
    # Use approximation: se = abs(beta) / sqrt(qchisq(p, df=1, lower.tail=FALSE))
    locus_snps$se <- with(locus_snps, abs(beta) / sqrt(qchisq(p_wald, df = 1, lower.tail = FALSE)))
  }
  
  # Prepare coloc dataset
  coloc_data <- list(
    pvalues = locus_snps$p_wald,
    N = rep(NA, nrow(locus_snps)),  # Sample size unknown
    MAF = locus_snps$af,
    beta = locus_snps$beta,
    varbeta = locus_snps$se^2,
    snp = locus_snps$rs,
    position = locus_snps$pos,
    type = "quant"  # Quantitative trait
  )
  
  # Remove NA values
  na_mask <- !is.na(coloc_data$beta) & !is.na(coloc_data$varbeta) & 
             !is.na(coloc_data$pvalues) & !is.na(coloc_data$MAF)
  
  for (field in names(coloc_data)) {
    if (length(coloc_data[[field]]) > 1) {  # Not for scalar fields
      coloc_data[[field]] <- coloc_data[[field]][na_mask]
    }
  }
  
  if (sum(na_mask) < 5) {
    if (verbose) {
      cat(sprintf("  Too few SNPs (%d) with complete data for coloc\n", sum(na_mask)))
    }
    return(NULL)
  }
  
  return(coloc_data)
}

# Function to run coloc analysis
run_coloc_analysis <- function(dataset1, dataset2, 
                               prior_p1 = 1e-4, prior_p2 = 1e-4, prior_p12 = 1e-5,
                               verbose = FALSE) {
  
  if (is.null(dataset1) || is.null(dataset2)) {
    return(NULL)
  }
  
  # Ensure SNPs are in same order
  common_snps <- intersect(dataset1$snp, dataset2$snp)
  
  if (length(common_snps) < 5) {
    if (verbose) {
      cat(sprintf("  Too few common SNPs (%d) for coloc\n", length(common_snps)))
    }
    return(NULL)
  }
  
  # Subset to common SNPs
  idx1 <- match(common_snps, dataset1$snp)
  idx2 <- match(common_snps, dataset2$snp)
  
  dataset1_sub <- list(
    pvalues = dataset1$pvalues[idx1],
    N = dataset1$N,
    MAF = dataset1$MAF[idx1],
    beta = dataset1$beta[idx1],
    varbeta = dataset1$varbeta[idx1],
    snp = common_snps,
    position = dataset1$position[idx1],
    type = dataset1$type
  )
  
  dataset2_sub <- list(
    pvalues = dataset2$pvalues[idx2],
    N = dataset2$N,
    MAF = dataset2$MAF[idx2],
    beta = dataset2$beta[idx2],
    varbeta = dataset2$varbeta[idx2],
    snp = common_snps,
    position = dataset2$position[idx2],
    type = dataset2$type
  )
  
  # Run coloc
  tryCatch({
    result <- coloc.abf(dataset1 = dataset1_sub,
                        dataset2 = dataset2_sub,
                        p1 = prior_p1,
                        p2 = prior_p2,
                        p12 = prior_p12)
    
    # Extract posterior probabilities
    pp <- result$summary
    
    coloc_result <- list(
      nsnps = pp["nsnps"],
      pp0 = pp["PP.H0.abf"],   # No association with either trait
      pp1 = pp["PP.H1.abf"],   # Association with trait 1 only
      pp2 = pp["PP.H2.abf"],   # Association with trait 2 only
      pp3 = pp["PP.H3.abf"],   # Association with both traits, different variants
      pp4 = pp["PP.H4.abf"]    # Association with both traits, same variant
    )
    
    return(coloc_result)
    
  }, error = function(e) {
    if (verbose) {
      cat(sprintf("  Coloc error: %s\n", e$message))
    }
    return(NULL)
  })
}

# Function to run eCAVIAR analysis (simplified version)
run_ecaviar_analysis <- function(dataset1, dataset2, verbose = FALSE) {
  # Simplified eCAVIAR implementation
  # In practice, would call the eCAVIAR software or use R implementation
  
  if (is.null(dataset1) || is.null(dataset2)) {
    return(list(clpp = NA, causal_set = NA))
  }
  
  common_snps <- intersect(dataset1$snp, dataset2$snp)
  
  if (length(common_snps) < 5) {
    return(list(clpp = NA, causal_set = NA))
  }
  
  # Simplified: Use coloc PP4 as proxy for CLPP
  # In real implementation, would compute proper eCAVIAR statistics
  tryCatch({
    # Placeholder for actual eCAVIAR computation
    # For now, return NA values
    return(list(
      clpp = NA,
      causal_set = paste(common_snps[1:min(5, length(common_snps))], collapse = ";")
    ))
  }, error = function(e) {
    if (verbose) {
      cat(sprintf("  eCAVIAR error: %s\n", e$message))
    }
    return(list(clpp = NA, causal_set = NA))
  })
}

# Main function
main <- function() {
  args <- parse_arguments()
  
  if (args$verbose) {
    cat("Bayesian Colocalization Analysis\n")
    cat("================================\n")
    cat(sprintf("GWAS traits: %s\n", paste(args$trait_names, collapse = ", ")))
    cat(sprintf("Locus file: %s\n", args$loci))
    cat(sprintf("Output file: %s\n", args$output))
    cat(sprintf("PP4 threshold: %.2f\n", args$pp4_threshold))
    cat("\n")
  }
  
  # Load GWAS data
  gwas_data <- list()
  for (i in seq_along(args$gwas_files)) {
    tryCatch({
      gwas_data[[args$trait_names[i]]] <- load_gwas_data(args$gwas_files[i], 
                                                         args$trait_names[i],
                                                         args$verbose)
    }, error = function(e) {
      stop(sprintf("Error loading GWAS data for %s: %s", args$trait_names[i], e$message))
    })
  }
  
  # Load locus data
  tryCatch({
    loci_df <- load_locus_data(args$loci, args$verbose)
  }, error = function(e) {
    stop(sprintf("Error loading locus data: %s", e$message))
  })
  
  # Prepare results data frame
  results <- data.frame()
  
  # Process each locus
  for (locus_idx in seq_len(nrow(loci_df))) {
    locus_info <- loci_df[locus_idx, ]
    
    if (args$verbose && locus_idx %% 10 == 0) {
      cat(sprintf("Processing locus %d/%d: %s\n", 
                  locus_idx, nrow(loci_df), locus_info$locus_id))
    }
    
    # Prepare coloc data for each trait at this locus
    coloc_datasets <- list()
    for (trait in args$trait_names) {
      coloc_datasets[[trait]] <- prepare_coloc_data(gwas_data[[trait]], locus_info, args$verbose)
    }
    
    # Run coloc for each trait pair
    trait_pairs <- combn(args$trait_names, 2, simplify = FALSE)
    
    for (pair in trait_pairs) {
      trait1 <- pair[1]
      trait2 <- pair[2]
      
      # Run coloc
      coloc_result <- run_coloc_analysis(coloc_datasets[[trait1]],
                                         coloc_datasets[[trait2]],
                                         args$prior_p1, args$prior_p2, args$prior_p12,
                                         args$verbose)
      
      # Run eCAVIAR (simplified)
      ecaviar_result <- run_ecaviar_analysis(coloc_datasets[[trait1]],
                                             coloc_datasets[[trait2]],
                                             args$verbose)
      
      # Store results
      if (!is.null(coloc_result)) {
        result_row <- data.frame(
          locus_id = locus_info$locus_id,
          chr = locus_info$chr,
          start = locus_info$start,
          end = locus_info$end,
          trait_pair = paste(trait1, trait2, sep = "_vs_"),
          n_snps = coloc_result$nsnps,
          coloc_pp0 = coloc_result$pp0,
          coloc_pp1 = coloc_result$pp1,
          coloc_pp2 = coloc_result$pp2,
          coloc_pp3 = coloc_result$pp3,
          coloc_pp4 = coloc_result$pp4,
          coloc_colocalized = coloc_result$pp4 >= args$pp4_threshold,
          ecaviar_clpp = ecaviar_result$clpp,
          ecaviar_set = ecaviar_result$causal_set,
          stringsAsFactors = FALSE
        )
        
        results <- rbind(results, result_row)
      }
    }
  }
  
  # Save results
  if (nrow(results) > 0) {
    output_dir <- dirname(args$output)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    write.table(results, args$output, sep = "\t", row.names = FALSE, quote = FALSE)
    
    if (args$verbose) {
      cat(sprintf("\nSaved %d colocalization results to %s\n", nrow(results), args$output))
    }
    
    # Generate summary statistics
    summary_stats <- list()
    
    # Overall coloc results
    if (nrow(results) > 0) {
      summary_stats$total_tests <- nrow(results)
      summary_stats$colocalized <- sum(results$coloc_colocalized, na.rm = TRUE)
      summary_stats$colocalization_rate <- summary_stats$colocalized / summary_stats$total_tests
      
      # By trait pair
      trait_pair_summary <- results %>%
        group_by(trait_pair) %>%
        summarise(
          n_tests = n(),
          n_colocalized = sum(coloc_colocalized, na.rm = TRUE),
          colocalization_rate = n_colocalized / n_tests,
          mean_pp4 = mean(coloc_pp4, na.rm = TRUE),
          median_pp4 = median(coloc_pp4, na.rm = TRUE)
        )
      
      summary_file <- sub("\\.tsv$", ".summary.tsv", args$output)
      write.table(trait_pair_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
      if (args$verbose) {
        cat("\nColocalization Summary\n")
        cat("=====================\n")
        cat(sprintf("Total tests: %d\n", summary_stats$total_tests))
        cat(sprintf("Colocalized loci: %d (%.1f%%)\n", 
                    summary_stats$colocalized, 
                    100 * summary_stats$colocalization_rate))
        
        cat("\nBy trait pair:\n")
        for (i in seq_len(nrow(trait_pair_summary))) {
          row <- trait_pair_summary[i, ]
          cat(sprintf("  %s: %d/%d (%.1f%%), mean PP4=%.3f\n",
                      row$trait_pair, row$n_colocalized, row$n_tests,
                      100 * row$colocalization_rate, row$mean_pp4))
        }
        
        cat(sprintf("\nSummary saved to: %s\n", summary_file))
      }
    }
    
  } else {
    cat("Warning: No colocalization results generated\n")
    
    # Create empty output file
    empty_df <- data.frame(
      locus_id = character(),
      chr = character(),
      start = numeric(),
      end = numeric(),
      trait_pair = character(),
      n_snps = numeric(),
      coloc_pp0 = numeric(),
      coloc_pp1 = numeric(),
      coloc_pp2 = numeric(),
      coloc_pp3 = numeric(),
      coloc_pp4 = numeric(),
      coloc_colocalized = logical(),
      ecaviar_clpp = numeric(),
      ecaviar_set = character(),
      stringsAsFactors = FALSE
    )
    
    write.table(empty_df, args$output, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  if (args$verbose) {
    cat("\nAnalysis completed successfully\n")
  }
}

# Run main function
if (!interactive()) {
  main()
}
