#' Genome Analysis for In Silico Amplicons
#'
#' The primary function for processing a FASTA file. It searches through contigs 
#' for amplicons products using specified primers, accounting for strand orientation 
#' and mismatches.
#'
#' @param genome_path String. Path to the FASTA file containing genome sequences or contigs.
#' @param forward_primers Character vector. Forward primer sequences (5'->3').
#' @param reverse_primers Character vector. Reverse primer sequences (5'->3').
#' @param max_mismatch Integer. Maximum allowed mismatches for primer binding (default is 1).
#' @param max_amplicon_length Integer. Maximum allowed length of the amplicons (default is 5000 bp).
#' @param min_amplicon_length Integer. Minimum allowed length of the amplicons (default is 100 bp).
#' @param output_dir String. Directory where the resulting FASTA files will be saved.
#' @param all Logical. If FALSE (default), the function stops after finding an amplicon in the first matching contig.
#' @param min_contig_length Integer. Minimum contig length to be considered for analysis.
#'
#' @return A list of objects containing extracted sequences (with and without primers) and their identifiers.
#' @export
#'

analyze_genome <- function(genome_path,
                           forward_primers,
                           reverse_primers,
                           max_mismatch = 1,
                           max_amplicon_length = 5000,
                           min_amplicon_length = 100,
                           output_dir = "contigs_results",
                           all = FALSE,
                           min_contig_length = 1000) {
  
  cat("=== GENOME ANALYSIS ===\n")
  genome_seqs <- readDNAStringSet(genome_path)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  raw_results <- list()
  file_label <- tools::file_path_sans_ext(basename(genome_path))
  
  for (i in seq_along(genome_seqs)) {
    contig_name <- names(genome_seqs)[i]
    if (is.null(contig_name) || contig_name == "") contig_name <- paste0("contig_", i)
    if (width(genome_seqs[i]) < min_contig_length) next
    
    current_seq <- DNAString(toupper(as.character(genome_seqs[[i]])))
    
    contig_results <- analyze_single_contig(
      sequence = current_seq,
      contig_name = contig_name,
      forward_primers = forward_primers,
      reverse_primers = reverse_primers,
      max_mismatch = max_mismatch,
      max_amplicon_length = max_amplicon_length,
      min_amplicon_length = min_amplicon_length
    )
    
    if (length(contig_results) > 0) {
      raw_results[[contig_name]] <- contig_results
      cat("*** FOUND", length(contig_results), "amplicons in:", contig_name, "***\n")
      if (!all) break
    }
  }
  
  if (length(raw_results) == 0) {
    cat("No amplicons found.\n")
    return(NULL)
  }
  
  final_list <- unlist(raw_results, recursive = FALSE)
  
  n_total <- length(final_list)
  new_names <- if(n_total == 1) file_label else paste0(file_label, "_", seq_len(n_total))
  names(final_list) <- new_names
  
  for(i in seq_along(final_list)) {
    final_list[[i]]$amplicon_id <- new_names[i]
  }
  
  save_final_results(final_list, output_dir)
  return(final_list)
}

#' Analyze Single Contig for Amplicons
#'
#' An internal helper function that searches for amplicons within a single DNAString sequence.
#' It handles four primer pairing combinations (forward-reverse on both strands).
#'
#' @param sequence A DNAString object representing the contig sequence.
#' @param contig_name String. The name of the contig.
#' @param forward_primers Character vector. Forward primers.
#' @param reverse_primers Character vector. Reverse primers.
#' @param max_mismatch Integer. Maximum allowed mismatches.
#' @param max_amplicon_length Integer. Maximum amplicon length.
#' @param min_amplicon_length Integer. Minimum amplicon length.
#'
#' @return A list of lists, each containing \code{with_p} (sequence with primers) and \code{no_p} (sequence without primers).
#' @keywords internal

analyze_single_contig <- function(sequence, contig_name,
                                  forward_primers, reverse_primers,
                                  max_mismatch = 1,
                                  max_amplicon_length = 5000,
                                  min_amplicon_length = 100) {
  
  amplicons <- list()
  
  f_p <- DNAString(forward_primers[1]); f_rc <- reverseComplement(f_p)
  r_p <- DNAString(reverse_primers[1]); r_rc <- reverseComplement(r_p)
  
  hits <- list(
    f_dir = matchPattern(f_p,  sequence, max.mismatch = max_mismatch),
    f_rc  = matchPattern(f_rc, sequence, max.mismatch = max_mismatch),
    r_dir = matchPattern(r_p,  sequence, max.mismatch = max_mismatch),
    r_rc  = matchPattern(r_rc, sequence, max.mismatch = max_mismatch)
  )
  
  find_pairs <- function(m1, m2, is_reverse) {
    if (length(m1) == 0 || length(m2) == 0) return()
    
    for (i in seq_along(m1)) {
      for (j in seq_along(m2)) {
        start_pos <- min(start(m1[i]), start(m2[j]))
        end_pos <- max(end(m1[i]), end(m2[j]))
        len <- end_pos - start_pos + 1
        
        if (len >= min_amplicon_length && len <= max_amplicon_length) {
          frag <- subseq(sequence, start_pos, end_pos)
          
          if (is_reverse) {
            frag_with <- reverseComplement(frag)
            no_p <- reverseComplement(subseq(sequence, min(end(m1[i]), end(m2[j])) + 1, max(start(m1[i]), start(m2[j])) - 1))
          } else {
            frag_with <- frag
            no_p <- subseq(sequence, min(end(m1[i]), end(m2[j])) + 1, max(start(m1[i]), start(m2[j])) - 1)
          }
          
          amplicons[[length(amplicons) + 1]] <<- list(
            with_p = frag_with,
            no_p = no_p
          )
        }
      }
    }
  }
  
  find_pairs(hits$f_dir, hits$r_rc, FALSE)
  find_pairs(hits$r_dir, hits$f_rc, TRUE)
  find_pairs(hits$r_rc, hits$f_rc, TRUE)
  find_pairs(hits$f_dir, hits$r_dir, FALSE)
  
  return(amplicons)
}

#' Export Results to FASTA Files
#'
#' Saves identified amplicon sequences into two separate FASTA files:
#' one with primer sequences included and one with only the internal sequence.
#'
#' @param final_list The resulting list from the \code{analyze_genome} function.
#' @param output_dir String. The output directory.
#'
#' @return No return value; writes files to disk.
#' @export

save_final_results <- function(final_list, output_dir) {
  headers <- names(final_list)
  seqs_with <- DNAStringSet(lapply(final_list, function(x) x$with_p))
  seqs_without <- DNAStringSet(lapply(final_list, function(x) x$no_p))
  
  names(seqs_with) <- headers
  names(seqs_without) <- headers
  
  writeXStringSet(seqs_with, file.path(output_dir, "amplicons_with_primers.fasta"))
  writeXStringSet(seqs_without, file.path(output_dir, "amplicons_without_primers.fasta"))
  
  cat("Analysis complete. Results saved in:", output_dir, "\n")
}