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
  
  MFS_F <- "AGCGGAGCGGGACGAATCGA"
  MFS_R <- "GGTTCCTTCCCGTCATACTCATTGCCATTGGAG"
  
  MFS_reference <- Biostrings::DNAString("AGCGGAGCGGGACGAATCGAGAACATGATCCCTGATCCGTTCCCGTTCGATCTCGCTTATTACGAGACGAGAAAAAGACAGCCACAGCCGCAAGGATTCGGACTTGACGACTTTTCGGTGGACAGGACAGATTCCGGACGGGGTAGGAATGCGAGGAGAGGGTAAGGTAGGTGAACACCTTATACTCCGTTTCTTTCCCATTCCTCTCCTCCATCATCACAGTCATCAGATCCGCACAACAATCATTCGGGTTCGGTTCCTTCCCGTCATACTCATTGCCATTGGAG")
  
  mfs1_db <- list(
    "Type_I"        = "519 bp Insert",
    "Type_II"       = "150–369 bp Insert",
    "Type_III"      = "149 bp Insert",
    "Rare_Or_Indel" = "Rare variant or specific Indel",
    "None"          = "No insertion"
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
  
  cat("=== Starting Comprehensive MFS1 Pipeline ===\n")
  cat("Total files to process:", length(fasta_files), "\n\n")
  
  for (f_path in fasta_files) {
    file_label <- tools::file_path_sans_ext(basename(f_path))
    cat("Processing:", file_label, "... ")
    
    amp_results <- analyze_genome(
      genome_path = f_path,
      forward_primers = c(MFS_F),
      reverse_primers = c(MFS_R),
      min_contig_length = 100,
      max_mismatch = 2,
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
    
    curr_seq <- amp_results[[1]]$with_p
    res_indels <- analyze_indels(curr_seq, MFS_reference)
    
    res_df <- data.frame(
      Sample_ID = file_label,
      Insert_Type = paste(res_indels$Insert_Types, collapse = ", "),
      Description = paste(sapply(res_indels$Insert_Types, function(x) mfs1_db[[x]] %||% x), collapse = "; "),
      Ins_Count = res_indels$Ins_Count,
      Ins_Lengths = paste(res_indels$Ins_Lengths, collapse = ", "),
      Del_Count = res_indels$Del_Count,
      Del_Lengths = paste(res_indels$Del_Lengths, collapse = ", "),
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

`%||%` <- function(x, y) if (length(x) == 0) y else x

analyze_indels <- function(amp_seq, ref_seq) {
  
  aln <- pwalign::pairwiseAlignment(
    pattern = amp_seq,
    subject = ref_seq,
    type = "global"
  )
  
  indel_info <- pwalign::indel(aln)
  
  insertions <- indel_info@insertion
  deletions  <- indel_info@deletion
  
  ins_widths <- Biostrings::width(insertions)
  del_widths <- Biostrings::width(deletions)
  
  num_insertions <- length(insertions)
  num_deletions  <- length(deletions)
  
  classify_single_indel <- function(w) {
    if (length(w) == 1) {
      if (w >= 480 & w <= 550) {
        return("Type_I")
      } else if (w >= 130 & w <= 370) {
        return("Type_II")
      } else if (abs(w - 149) < 15) {
        return("Type_III")
      } else if (w > 20) {
        return("Rare_Or_Indel")
      } else {
        return("None")
      }
    } else {
      sapply(w, function(x) {
        if (x >= 480 & x <= 550) {
          "Type_I"
        } else if (x >= 130 & x <= 370) {
          "Type_II"
        } else if (abs(x - 149) < 15) {
          "Type_III"
        } else if (x > 20) {
          "Rare_Or_Indel"
        } else {
          "None"
        }
      })
    }
  }
  
  if (num_insertions > 0) {
    insert_types <- classify_single_indel(ins_widths)
  } else {
    insert_types <- "None"
  }
  
  return(list(
    Insert_Types = insert_types,
    Ins_Count = num_insertions,
    Ins_Lengths = ins_widths,
    Del_Count = num_deletions,
    Del_Lengths = del_widths
  ))
}