#' ZymoR Pipeline: Batch Process Genomes
#'
#' Scans a directory for FASTA files, extracts amplicons, and identifies haplotypes
#' for all samples, merging everything into one final report.
#'
#' @param input_dir String. Path to folder containing .fasta, .fna, or .txt files.
#' @param forward_primers Character vector. Primer sequences.
#' @param reverse_primers Character vector. Primer sequences.
#' @param reference_seq DNAString. Reference sequence for alignment.
#' @param target_positions Named vector. Positions to check (e.g., c(L50=50)).
#' @param cds_ranges IRanges. Ranges of exons.
#' @param haplotype_db List. Database of known haplotypes.
#' @param output_dir String. Where to save amplicons and the final table.
#' @param ... Additional arguments passed to analyze_genome (e.g., max_mismatch).
#'
#' @export

run_zymor_pipeline <- function(input_dir, 
                               forward_primers, 
                               reverse_primers, 
                               reference_seq, 
                               target_positions, 
                               cds_ranges, 
                               haplotype_db,
                               output_dir = "zymor_results",
                               ...) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  fasta_files <- list.files(input_dir, pattern = "\\.(fasta|fna|txt|fa)$", full.names = TRUE)
  if (length(fasta_files) == 0) stop("No sequence files found in input_dir.")
  
  all_mutation_results <- list()
  not_found_list <- c()
  
  global_amplicons_with <- list()
  global_amplicons_no   <- list()

  cat("=== Starting ZymoR Pipeline ===\n")
  cat("Total files to process:", length(fasta_files), "\n\n")

  for (f_path in fasta_files) {
    file_label <- tools::file_path_sans_ext(basename(f_path))
    
    sample_folder <- file.path(output_dir, file_label)
    if (!dir.exists(sample_folder)) dir.create(sample_folder, recursive = TRUE)
    
    cat("Processing:", file_label, "... ")

    amp_results <- analyze_genome(
      genome_path = f_path,
      forward_primers = forward_primers,
      reverse_primers = reverse_primers,
      output_dir = sample_folder, 
      ...
    )
    
    if (is.null(amp_results)) {
      cat("[NOT FOUND]\n")
      not_found_list <- c(not_found_list, file_label)
      unlink(sample_folder, recursive = TRUE)
      next
    }
    
    cat("[FOUND]\n")
    
    for(res in amp_results) {
      global_amplicons_with[[res$amplicon_id]] <- res$with_p
      global_amplicons_no[[res$amplicon_id]]   <- res$no_p
    }
    
    mut_matrix <- get_mutation_matrix(
      analyze_genome_results = amp_results,
      reference_seq = reference_seq,
      target_positions = target_positions,
      cds_ranges = cds_ranges,
      haplotype_db = haplotype_db
    )
    
    write.csv(mut_matrix, file.path(sample_folder, paste0(file_label, "_mutations.csv")), row.names = FALSE)
    all_mutation_results[[file_label]] <- mut_matrix
  }
  
  
  cat("\n--- Finalizing Global Results ---\n")

  if(length(global_amplicons_with) > 0) {
    all_with <- Biostrings::DNAStringSet(global_amplicons_with)
    Biostrings::writeXStringSet(all_with, file.path(output_dir, "ALL_SAMPLES_with_primers.fasta"))
    cat("- Saved: ALL_SAMPLES_with_primers.fasta\n")
  }

  if(length(global_amplicons_no) > 0) {
    all_no <- Biostrings::DNAStringSet(global_amplicons_no)
    Biostrings::writeXStringSet(all_no, file.path(output_dir, "ALL_SAMPLES_without_primers.fasta"))
    cat("- Saved: ALL_SAMPLES_without_primers.fasta\n")
  }

  if (length(all_mutation_results) > 0) {
    final_table <- do.call(rbind, all_mutation_results)
    rownames(final_table) <- NULL
    write.csv(final_table, file.path(output_dir, "FINAL_MUTATION_REPORT.csv"), row.names = FALSE)
    cat("- Saved: FINAL_MUTATION_REPORT.csv\n")
  }

  if (length(not_found_list) > 0) {
    writeLines(not_found_list, file.path(output_dir, "MISSING_AMPLICONS_LOG.txt"))
    cat("- Warning: Gene not found in", length(not_found_list), "samples (see MISSING_AMPLICONS_LOG.txt)\n")
  }
  
  cat("\n=== Pipeline Complete ===\n")
  return(final_table)
}