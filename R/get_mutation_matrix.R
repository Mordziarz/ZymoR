#' Generate Mutation Matrix and Haplotype Identification
#'
#' Processes amplicon sequences to detect mutations at specific amino acid positions, 
#' maps them to a reference genome, and identifies the resulting haplotype based 
#' on a provided database. This function handles negative strand orientations 
#' by reversing coordinates and calculating reverse complements for translation.
#'
#' @param analyze_genome_results List. The output from \code{analyze_genome}, containing 
#'   amplicon sequences (\code{with_p}) and identifiers (\code{amplicon_id}).
#' @param reference_seq DNAString or DNAStringSet. The reference sequence used 
#'   for coordinate mapping and pairwise alignment.
#' @param target_positions Named numeric vector. Key-value pairs where names are 
#'   mutation labels (e.g., "L50") and values are the amino acid positions (e.g., 50).
#' @param cds_ranges IRanges. The coding sequence ranges. Coordinates should 
#'   reflect the genomic positions of exons.
#' @param haplotype_db List. A database of known haplotypes where each entry 
#'   is a named character vector of expected mutations.
#'
#' @details 
#' The function performs several steps for each sample:
#' \enumerate{
#'   \item Maps amino acid positions to genomic DNA coordinates using \code{cds_ranges}.
#'   \item Performs a global-local pairwise alignment between the amplicon and reference.
#'   \item Extracts the corresponding codon, handles reverse complementation (for negative strands), 
#'         and translates it to an amino acid.
#'   \item Compares the result against the wild-type (inferred from the first letter of the position name).
#'   \item Matches the final mutation profile against the \code{haplotype_db}.
#' }
#'
#' @return A data frame (matrix) where rows are samples and columns include 
#'   Sample_ID, Haplotype name, Status (KNOWN/NEW), and individual columns 
#'   for each analyzed mutation position containing the detected amino acid or "wt".
#' @export
#'

get_mutation_matrix <- function(analyze_genome_results, 
                                reference_seq, 
                                target_positions, 
                                cds_ranges, 
                                haplotype_db) {
  
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
    if(max(t_pos) > length(ref_idx)) return(list(aa="EMPTY"))
    t_aln <- ref_idx[t_pos]
    raw_ex <- substr(pat_aln, min(t_aln), max(t_aln))
    if(nchar(raw_ex) < 3 || grepl("-", raw_ex)) return(list(aa="DEL"))
    rc_codon <- as.character(Biostrings::reverseComplement(DNAString(raw_ex)))
    amino <- as.character(Biostrings::translate(DNAString(rc_codon)))
    return(list(aa = amino))
  }
  
  mut_names <- names(target_positions)
  indel_col <- "INDEL_Promoter"
  col_names <- c("Sample_ID", "Haplotype", "Status", indel_col, mut_names)
  
  res_table <- data.frame(matrix(NA, nrow = length(analyze_genome_results), ncol = length(col_names)))
  colnames(res_table) <- col_names
  
  is_cyp51 <- identical(target_positions, CYP51_target_positions)

  for (i in seq_along(analyze_genome_results)) {
    amp_data <- analyze_genome_results[[i]]
    curr_amp <- amp_data[["with_p"]]
    curr_id  <- as.character(amp_data[["amplicon_id"]])
    
    res_table[i, "Sample_ID"] <- curr_id
    detected_muts <- list()
    
    if (is_cyp51) {
      aln_for_indel <- pwalign::pairwiseAlignment(
        pattern = DNAString(as.character(curr_amp)), 
        subject = reference_seq, 
        type = "global-local"
      )
      end_in_amp <- Biostrings::end(pwalign::pattern(aln_for_indel))
      amp_total_len <- nchar(as.character(curr_amp))
      trailing_seq_len <- amp_total_len - end_in_amp
      
      if (trailing_seq_len > 150) {
        res_table[i, indel_col] <- "YES"
      } else {
        res_table[i, indel_col] <- "NO"
      }
    } else {
      res_table[i, indel_col] <- "N/A"
    }

    for (m_name in mut_names) {
      aa_num <- target_positions[[m_name]]
      dna_idx <- get_neg_pos(c((aa_num*3)-2, (aa_num*3)-1, (aa_num*3)), cds_ranges)
      
      mut_res <- check_mut(curr_amp, reference_seq, dna_idx)
      wild_aa <- substr(m_name, 1, 1)
      
      if (mut_res$aa != wild_aa && mut_res$aa != "EMPTY") {
        res_table[i, m_name] <- mut_res$aa
        detected_muts[[m_name]] <- mut_res$aa
      } else {
        res_table[i, m_name] <- "wt"
      }
    }
    
    final_detected <- unlist(detected_muts)
    
    h_name <- identify_haplotype(final_detected, haplotype_db)
    res_table[i, "Haplotype"] <- h_name
    
    is_known <- h_name %in% names(haplotype_db) || h_name == "WildType (A0)"
    res_table[i, "Status"] <- if(is_known) "KNOWN" else "NEW"
  }
  
  return(res_table)
}

#' Identify Haplotype Name from Detected Mutations
#'
#' Internal helper to match a vector of detected mutations against the database.
#'
#' @param detected_muts Named character vector of mutations found in the sample.
#' @param haplotype_db List of known haplotypes.
#' @return String. The name of the haplotype or "WildType (A0)" / "Unknown".
#' @keywords internal

identify_haplotype <- function(detected_muts, haplotype_db) {
  
  if (length(detected_muts) == 0) {
    return("WildType (A0)")
  }
  
  detected_muts <- detected_muts[order(names(detected_muts))]
  
  for (h_name in names(haplotype_db)) {
    db_muts <- haplotype_db[[h_name]]
    db_muts <- db_muts[order(names(db_muts))]
    
    if (length(detected_muts) == length(db_muts)) {
      if (all(names(detected_muts) == names(db_muts)) && 
          all(detected_muts == db_muts)) {
        return(h_name)
      }
    }
  }
  
  return("Unknown")
}