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

run_zymo_pipeline <- function(input_dir, 
                              forward_primers, 
                              reverse_primers, 
                              reference_seq, 
                              target_positions, 
                              cds_ranges, 
                              haplotype_db,
                              output_dir = "zymor_results",
                              ...) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  fasta_files <- list.files(input_dir, pattern = "\\.(fasta|fna|txt)$", full.names = TRUE)
  
  if (length(fasta_files) == 0) stop("No sequence files found in input_dir.")
  
  all_mutation_results <- list()
  not_found_list <- c()

  cat("Starting ZymoR Pipeline on", length(fasta_files), "files...\n")

  for (f_path in fasta_files) {
    file_name <- basename(f_path)
    cat("\n--- Processing:", file_name, "---\n")

    amp_results <- analyze_genome(
      genome_path = f_path,
      forward_primers = forward_primers,
      reverse_primers = reverse_primers,
      output_dir = file.path(output_dir, "amplicons"),
      ...
    )
    
    if (is.null(amp_results)) {
      cat("! Warning: No amplicon found for", file_name, "\n")
      not_found_list <- c(not_found_list, file_name)
      next
    }
    
    mut_matrix <- get_mutation_matrix(
      analyze_genome_results = amp_results,
      reference_seq = reference_seq,
      target_positions = target_positions,
      cds_ranges = cds_ranges,
      haplotype_db = haplotype_db
    )
    
    all_mutation_results[[file_name]] <- mut_matrix
  }
  
  if (length(all_mutation_results) > 0) {
    final_table <- do.call(rbind, all_mutation_results)
    rownames(final_table) <- NULL
    
    write.csv(final_table, file.path(output_dir, "final_haplotype_report.csv"), row.names = FALSE)
    cat("\nSuccess! Final table saved to:", file.path(output_dir, "final_haplotype_report.csv"), "\n")
  } else {
    final_table <- NULL
    cat("\nNo haplotypes identified in any file.\n")
  }

  if (length(not_found_list) > 0) {
    writeLines(not_found_list, file.path(output_dir, "not_found_log.txt"))
    cat("Files with no amplicons logged in 'not_found_log.txt'\n")
  }
  
  return(final_table)
}