#' Generate MFS1 Promoter Analysis Table
#'
#' Processes a single FASTA file or a directory of FASTA files to analyze the MFS1
#' promoter region for the presence of retrotransposon insertions (Types I, II, III).
#'
#' @param input_path String. Path to a single file (.fasta, .fna, .txt) or a folder 
#'   containing sequence files.
#' @param output_dir String. Where to save results and tables. Defaults to "zymor_results".
#' @param ... Additional arguments passed to analyze_genome (e.g., max_mismatch).
#'
#' @return A data frame where rows are samples and columns include 
#'   Sample_ID, Insert_Type oraz Description.
#' @export
#'

get_MFS1 <- function(input_path, output_dir = "zymor_results", ...) {
  
  MFS_F <- "GGAGCGGAGCGGGACGAATCGAGA"
  MFS_R <- "GAGGCCAGTCGAGGCTTTTACATA"
  
  MFS_reference <- Biostrings::DNAString("GGAGCGGAGCGGGACGAATCGAGAACATGATCCCTGATCCGTTCCCGTTCGATCTCGCTTATTACGAGACGAGAAAAAGACAGCCACAGCCGCAAGGATTCGGACTTGACGACTTTTCGGTGGACAGGACAGATTCCGGACGGGGTAGGAATGCGAGGAGAGGGTAAGGTAGGTGAACACCTTATACTCCGTTTCTTTCCCATTCCTCTCCTCCATCATCACAGTCATCAGATCCGCACAACAATCATTCGGGTTCGGTTCCTTCCCGTCATACTCATTGCCATTGGAGCTGAAAGGTCGAATTGCTGCTGCCTTGGGTCTCGCTTCGCTTCTGTGATTCGACTTGTTCCGCCGTTGGATCTGGTGGAACGGAGATATTCGACCCCGGAATCCCCGAGTCCCGTCTGCTCGACAAACAGCCCGCCTTTGCATACAAGTAAGGTATCACCTCACTTCTCAGGTTTGTAGCATTCACCACTTCGACCAGGTCAGGTCAGGCCCACCCCCGTCAACCCTCCAACAACCATACCTGGTAGTTCTACCCATCCACCGACTTCGTCATCGACGATACCGGCAGATCGCTCATACGAGCACGACTTCACCTCTCATCACCACCCCACACCGCCACCACAACCATGCGCTTCTTCAGCAAGAAGTCGCCGGCCGCGGTCGTCGACTCCACCGCACCCTCATCATCCCACCACTTGGACAGCGACGCCGAGCACTACAACGAAAAGCCCGAACCCGTCGAAGACATCACCACGCCCGTGCGAAACTCTTCCGAAGTCAACAATCCCGCCGAAGCCAAAGTCGAAGCCATCGAAAATGAAGCCTTGTCGAGACAATTCTCCAAGGACTACCCGACTGGTCTGAAACTCACCCTCATCACTCTGGCATTATGTTTGTCCGTGTTGTGTATGGCGTTGGACAATACCATCATCAGCACGGCTATTCCGAGGATCACGGATCAGTTTAATTCGTATGTTCGCATTGAAGCAGCCTGCTTGCTATGACTTACTAACATCCACCTCCTCCAGCATCGACGACGTAGGCTGGTACGGCTCCTCCTACCTCCTCACAACCTGCGCACTACAGCTCTTCTTCGGTCAGTCCTCTCAAGCTGCATGAGGCCTACGGCAAAGGTCAACCATCGCTCACACACTCTTTCACGTCAGGCAAACTCTACACCTTCTACTCCATCAAGACCATCTACCTCACCGCCCTCGTCATCTTCGAGATCGGCTCGGCGGTATGCGGTGCAGCGCCCAACTCGATCGCGCTCATCATCGGCCGTGCGATTGCCGGTGTCGGGTCCGCGGGTATTTTCTCGGGCGCCATTCTCATCGTCGCCAATACCGTTCCGCTCGAGAAACGTCCCATTTACACCGGGTTTATTGGTGGCATGTATGGTATTGCTAGTGTCGTTGGTCCTTTGATGGGTGGTGGTGAGTACAGTCCCAGCCTACTCTTGCAGTAGCAATCCTACTGACACGAGTTTAGCATTCACCGACCACGTCTCATGGAGATGGTGAGTATTGTTCCAGTATTCTCTTGCAGCCTCGCCCTTACTAACATGACATGCAGGTGCTTCTACATCAACCTCCCTCTAGGCGCCATCACCATCGCCTTCATAGTCTTCTTCTACCACCCCACCAAATCCGCCGCCCAAACCATCGCCACAACCTGGCGCTCCAAACTCCAACAATTCGACATCTTCGGCACCCTCGTCTTCCTCCCCATGATCGTCTGCCTCCTCCTCGCCCTCCAATGGGGCGGCAGCAAATACCCCTGGTCGGACGGAAAAGTCATCGCCCTCTTCGTCGTCTTCGGCGTCCTCCTCCTCATCTTCGCCACCATCCAAGTCCTCAAAGCCGACCACGCCACCGTCCCACCCCGCCTCCTCCGCAACCGCACCGTCGCAGGAGCCTGCTGGTTCGTCCTCTGCCTCGGCGGCGCCTTTTTCATCTTCATCTACTACCTCCCCCTCTGGTTCCAGGCCATCAAAGGCGCCGCTCCCACCAAATCCGGAATCATGAACATCCCCCTCGTCCTCGGCTTGGTCCTCATGTCCATCCTCTCCGGCGGGCTGATTACCGTAATCGGATACTACACGCCCTTTGTCATCGCCTCGTCAATCTTCATGGCCATTGGAGCGGGTATGCTGTCGACATTCACACCCGCCACCGGCTCCCCGAAATGGATCGGCTTCCAAGTCCTCTTCGGCATCGGCGTGGGACTGGGAATGCAACAACCTCTCATCGCCGTGCAGGCCGTTTTACCCGCTAAAGACGTCCCGATCGGCACGGCGATGATCATGTTCAGCCAAACCTTGGGAGGGGCGCTGTTCATCGCTGTCGCGCAGAATGTCTTCTCCACGCAACTACTCAAAAACTTGCGATCGATCGTCCCCGATTTGGATCCGAGGATTGTGTTGGAAGCGGGCGCCACGACGTTGGTGAAGAGTGTACCCGGAGAGTCATTGATGGGGGTGAGGTTGGCGTACAATGATTCGATCATGTCTGCGTTTTTGGTGGCGACGGCGATGGCGACGATTAGTGTTCTGGGAGCGGTGGCGTTTGAGTGGAAGAGCGTCAAGGGGAAGAAGATTGAGATGGGGGGTGCGTGAGGGAGGAATGATCCAGGAGGAATGACCCATGAGGAAGAGGTGGAGGGGTTGGGTGAGCGTGGAGGTGGAGACTGTATGATGATACCCTTGTTTTGCTTTTTCTTTCTGCGACGATGGAGTATGAGCGGGAGGTTTGGAGTTTTGTCGGGACCTCGCATTGGCGCCGCGATACAGCGCGGGTCTAGAGGCCAGTCGAGGCTTTTACATA")
  mfs1_db <- list(
    "Type_I"   = "519 bp Insert (LTR retrotransposon), strong inducer",
    "Type_II"  = "150–369 bp Insert, moderate inducer",
    "Type_III" = "149 bp Insert, weak inducer",
    "Rare_Or_Indel" = "Rare variant or specific Indel in 5' UTR",
    "None"     = "No insertion"
  )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  is_dir <- dir.exists(input_path)
  fasta_files <- c()
  
  if (is_dir) {
    fasta_files <- list.files(input_path, pattern = "\\.(fasta|fna|txt|fa)$", full.names = TRUE)
  } else {
    if (file.exists(input_path)) {
      fasta_files <- c(input_path)
    } else {
      stop("Provided file or directory does not exist.")
    }
  }
  
  if (length(fasta_files) == 0) {
    stop("No sequence files found in the specified input path.")
  }
  
  all_results <- list()
  not_found_list <- c()
  
  global_amplicons_with <- list()
  global_amplicons_no   <- list()
  
  cat("=== Starting MFS1 Pipeline ===\n")
  cat("Total files to process:", length(fasta_files), "\n\n")
  
  for (f_path in fasta_files) {
    file_label <- tools::file_path_sans_ext(basename(f_path))
    cat("Processing:", file_label, "... ")
    
    amp_results <- analyze_genome(
      genome_path = f_path,
      forward_primers = c(MFS_F),
      reverse_primers = c(MFS_R),
      ...
    )
    
    if (is.null(amp_results)) {
      cat("[NOT FOUND]\n")
      not_found_list <- c(not_found_list, file_label)
      next
    }
    
    cat("[FOUND]\n")
    
    if (length(amp_results) > 1) {
      for (i in seq_along(amp_results)) {
        res <- amp_results[[i]]
        sample_amp_id <- paste0(file_label, "_", i)
        
        if (!is.null(res$with_p)) {
          global_amplicons_with[[sample_amp_id]] <- res$with_p
        }
        if (!is.null(res$no_p)) {
          global_amplicons_no[[sample_amp_id]] <- res$no_p
        }
      }
    } else {
      res <- amp_results[[1]]
      sample_amp_id <- file_label
      
      if (!is.null(res$with_p)) {
        global_amplicons_with[[sample_amp_id]] <- res$with_p
      }
      if (!is.null(res$no_p)) {
        global_amplicons_no[[sample_amp_id]] <- res$no_p
      }
    }
    
    curr_len <- nchar(as.character(amp_results[[1]]$with_p))
    ref_len <- nchar(as.character(MFS_reference))
    diff_len <- curr_len - ref_len
    
    if (diff_len > 480 && diff_len < 550) {
      insert_type <- "Type_I"
    } else if (diff_len >= 148 && diff_len <= 370) {
      if (abs(diff_len - 149) < 5) {
        insert_type <- "Type_III"
      } else {
        insert_type <- "Type_II"
      }
    } else if (abs(diff_len) < 20) {
      insert_type <- "None"
    } else {
      insert_type <- "Rare_Or_Indel"
    }
    
    res_df <- data.frame(
      Sample_ID = file_label,
      Insert_Type = insert_type,
      Description = mfs1_db[[insert_type]],
      Diff_Length = diff_len,
      stringsAsFactors = FALSE
    )
    
    all_results[[file_label]] <- res_df
  }
  
  cat("\n--- Finalizing Global Results ---\n")
  
  if (length(global_amplicons_with) > 0) {
    seq_list_with <- Biostrings::DNAStringSet(global_amplicons_with)
    Biostrings::writeXStringSet(seq_list_with, file.path(output_dir, "MFS_with_primers.fasta"))
    cat("- Saved: MFS_with_primers.fasta\n")
  }
  
  if (length(global_amplicons_no) > 0) {
    seq_list_no <- Biostrings::DNAStringSet(global_amplicons_no)
    Biostrings::writeXStringSet(seq_list_no, file.path(output_dir, "MFS_without_primers.fasta"))
    cat("- Saved: MFS_without_primers.fasta\n")
  }
  
  if (length(all_results) > 0) {
    final_table <- do.call(rbind, all_results)
    rownames(final_table) <- NULL
    write.csv(final_table, file.path(output_dir, "MFS1_analysis_table.csv"), row.names = FALSE)
    cat("- Saved: MFS1_analysis_table.csv\n")
  } else {
    final_table <- NULL
  }
  
  if (length(not_found_list) > 0) {
    writeLines(not_found_list, file.path(output_dir, "MFS1_missing_amplicons.txt"))
    cat("- Warning: Amplicon not found in", length(not_found_list), "samples\n")
  }
  
  cat("\n=== Complete ===\n")
  return(final_table)
}