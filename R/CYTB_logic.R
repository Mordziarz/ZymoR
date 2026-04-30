#' Generate Mutation Matrix for CYTB
#'
#' Processes a single FASTA file or a directory of FASTA files to detect mutations 
#' at specific amino acid positions within the CYTB gene and maps them to a reference genome.
#'
#' @param input_path String. Path to a single file (.fasta, .fna, .txt) or a folder 
#'   containing sequence files.
#' @param output_dir String. Where to save results and tables. Defaults to "zymor_results".
#' @param ... Additional arguments passed to analyze_genome (e.g., max_mismatch).
#'
#' @return A data frame where rows are samples and columns include 
#'   Sample_ID and individual columns for each analyzed CYTB mutation position.
#' @export
#'

get_CYTB <- function(input_path, output_dir = "zymor_results", ...) {

CYTB_F <- "ATGTCAATAGGCATGCAACCGT"
CYTB_R <- "ATACATAGATCATGTAACTATA"
  
CYTB_reference <- Biostrings::DNAString("TTATTTTACGTTATTTGTGTTGTTTAAGTGCATATCCACTAAAGTATTCTCGATTAAGCTTACTAGAGGCACTATGATTAAGAAATGTGAAAAGTACAGTACTGTACTAATTTGCCCAAATTCTATAAATGGTGACTCAACGTGTTTAGCACCTAATTGCATTAATACTAAGAAGTTAGCAACAAATATGTAAAATGTGATCTTGCTTAAAGGTCTAAATTGTAATCCTCTACTTCTACCTAAATCGGTAATTGGCATTATCATTATAATCAATATAGCAGAAAACATTGCGATCACACCTAATAGTTTGTTAGGTATTGATCTTAGTATAGCGTAGAAAGGTAATAAATATCATTCAGGCACTATAGCAGGTGGAGTTTGCATAGGGTTAGCCATAACATAATTCTCGCTGTCACCTAAAACGTTAGGCATAAAGAAAACAAATATTGATAACACTATAATGAATAAAAATATTGTTATTAAATCTTTGAATATAAAGTAAGGGGCGAATGGTAATCTATCGTAGTTACCTGATACACCTAAAGGATTTCCTGAACCCGCTGTATCGTGTAAAGCTATTAGATGCATTAAAACTAATGCAGCTAACACAAACGGTAAAACGAAATGTAGAGCAAAGAATCTGTTCAATGTTGCATTGTTAACAGAAAATCCACCTCATACGAATTCAACTATGTCTTGTCCAACTCAAGGTATTGCACTCAATAAGTTAGTTATAACTGTTGCTCCTCATAAAGACATTTGACCATAAGGTAATACATACCCTAAGAATGCGGTTGCCATCATCAGAACTAGTATTATAGTACCGATTGTTCATGTTAATGTTCTAGGGGCTTTGTATGACCCGTAGTATAACCCTCTACCCACGTGCAGGTAAACTAAAAAGAAGAAAGCTGAAGCAGTATTAGAGTGAAGATAACGTATAAGTCATCCATTGTTTACATCTCTCATTATGTGTTCAACTGAGTTAAATGCTTCAGATACGCTAGGGTTATAATGCATAGCAAGTGTAACACCAGTAACGATTTGTATTACAAGACAAAATCCCAGTAATGAACCGAAATTTCATAGATAGCTTAGATTACTTGGTTGTGGTGAATCGATAAGATAACCATTAACCAGACTAAATAAAGGGTGACTCTTTCAAATCCTCAT")

CYTB_CDS <- IRanges::IRanges(
  start = c(1),
  end   = c(1173),
  names = c("CDS1")
)
  
CYTB_target_positions <- c(
F129L=129,Y132C=132,G143A=143
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
  
  all_mutation_results <- list()
  not_found_list <- c()
  
  global_amplicons_with <- list()
  global_amplicons_no   <- list()
  
  cat("=== Starting get_CYTB Pipeline ===\n")
  cat("Total files to process:", length(fasta_files), "\n\n")
  
  for (f_path in fasta_files) {
    file_label <- tools::file_path_sans_ext(basename(f_path))
    
    cat("Processing:", file_label, "... ")
    
    amp_results <- analyze_genome(
      genome_path = f_path,
      forward_primers = c(CYTB_F),
      reverse_primers = c(CYTB_R),
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
    
    mut_matrix <- .internal_cytb_processing(
      analyze_genome_results = amp_results,
      CYTB_reference = CYTB_reference,
      cds_ranges = CYTB_CDS,
      CYTB_target_positions = CYTB_target_positions
    )
    
    all_mutation_results[[file_label]] <- mut_matrix
  }
  
  cat("\n--- Finalizing Global Results ---\n")
  
  if (length(global_amplicons_with) > 0) {
    seq_list_with <- Biostrings::DNAStringSet(global_amplicons_with)
    Biostrings::writeXStringSet(seq_list_with, file.path(output_dir, "CYTB_with_primers.fasta"))
    cat("- Saved: CYTB_with_primers.fasta\n")
  }
  
  if (length(global_amplicons_no) > 0) {
    seq_list_no <- Biostrings::DNAStringSet(global_amplicons_no)
    Biostrings::writeXStringSet(seq_list_no, file.path(output_dir, "CYTB_without_primers.fasta"))
    cat("- Saved: CYTB_without_primers.fasta\n")
  }
  
  if (length(all_mutation_results) > 0) {
    final_table <- do.call(rbind, all_mutation_results)
    rownames(final_table) <- NULL
    write.csv(final_table, file.path(output_dir, "CYTB_mutation_table.csv"), row.names = FALSE)
    cat("- Saved: CYTB_mutation_table.csv\n")
  } else {
    final_table <- NULL
  }
  
  if (length(not_found_list) > 0) {
    writeLines(not_found_list, file.path(output_dir, "CYTB_missing_amplicons.txt"))
    cat("- Warning: Gene not found in", length(not_found_list), "samples (see CYTB_missing_amplicons.txt)\n")
  }
  
  cat("\n=== Complete ===\n")
  return(final_table)
}

.internal_cytb_processing <- function(analyze_genome_results, CYTB_reference, cds_ranges, CYTB_target_positions) {
  
  get_neg_pos <- function(nt_idx, ranges) {
    all_pts <- c()
    for(i in seq_along(ranges)) {
      all_pts <- c(all_pts, seq(end(ranges[i]), start(ranges[i])))
    }
    return(all_pts[nt_idx])
  }
  
  check_mut <- function(amp, ref, t_pos) {
    aln <- pwalign::pairwiseAlignment(
      pattern = DNAString(as.character(amp)), 
      subject = DNAString(as.character(ref)), 
      type = "global-local"
    )
    
    sub_aln <- as.character(pwalign::subject(aln))
    pat_aln <- as.character(pwalign::pattern(aln))
    ref_idx <- which(strsplit(sub_aln, "")[[1]] != "-")
    
    if(max(t_pos) > length(ref_idx)) return(list(aa = "EMPTY"))
    
    t_aln <- ref_idx[t_pos]
    pat_vec <- strsplit(pat_aln, "")[[1]]
    raw_ex <- paste0(pat_vec[t_aln], collapse = "")
    
    if(nchar(raw_ex) < 3 || grepl("-", raw_ex)) return(list(aa = "DEL"))
    
    codon_dna <- Biostrings::complement(DNAString(raw_ex))
    
    amino <- as.character(Biostrings::translate(
      codon_dna, 
      genetic.code = Biostrings::getGeneticCode("1"), 
      if.fuzzy.codon = "solve",
      no.init.codon = TRUE
    ))
    
    return(list(aa = amino))
  }
  
  mut_names <- names(CYTB_target_positions)
  col_names <- c("Sample_ID", mut_names)
  
  res_table <- data.frame(matrix(NA, nrow = length(analyze_genome_results), ncol = length(col_names)))
  colnames(res_table) <- col_names
  
  for (i in seq_along(analyze_genome_results)) {
    amp_data <- analyze_genome_results[[i]]
    curr_amp <- amp_data[["with_p"]]
    curr_id  <- as.character(amp_data[["amplicon_id"]])
    
    res_table[i, "Sample_ID"] <- curr_id
    
    for (m_name in mut_names) {
      aa_num <- CYTB_target_positions[[m_name]]
      dna_idx <- get_neg_pos(c((aa_num*3)-2, (aa_num*3)-1, (aa_num*3)), cds_ranges)
      
      mut_res <- check_mut(curr_amp, CYTB_reference, dna_idx)
      wild_aa <- substr(m_name, 1, 1)
      
      if (mut_res$aa != wild_aa && mut_res$aa != "EMPTY") {
        res_table[i, m_name] <- mut_res$aa
      } else {
        res_table[i, m_name] <- "wt"
      }
    }
  }
  
  return(res_table)
}