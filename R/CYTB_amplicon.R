#' Analyze CYTB gene mutations from BAM files
#'
#' This function processes BAM files to identify specific mutations in the CYTB gene.
#' It compares observed amino acids against a defined reference and a known mutation 
#' profile for each position. It provides a summary count of reads per classification status.
#'
#' @param input_bam A character string representing the file path to a single 
#'   BAM file or a directory containing multiple BAM files.
#'
#' @return A \code{data.table} containing the following columns:
#'   \item{Sample_ID}{The name of the analyzed sample.}
#'   \item{Status}{Classification status: "KNOWN" (matches reference or known mutation), 
#'   "NEW" (any other amino acid or deletion), or "Cant classify" (if reads fall outside coverage).}
#'   \item{N}{The number of reads assigned to the specific mutation combination and status.}
#'   \item{F129L, Y132C, G143A}{Columns representing the detected amino acid (or "wt", "Del", "OUT") 
#'   at each target position.}
#'
#' @details
#' The function maps CYTB gene positions and performs translation of codons using the standard genetic code.
#' - "wt" is assigned if the amino acid matches the reference.
#' - The specific mutation (e.g., "L" for F129L) is labeled as defined in the target profile.
#' - Any other amino acid or deletion ("Del") is classified as "NEW".
#' - If a read does not cover all target positions, it is marked as "OUT" and the status 
#'   is set to "Cant classify" for that specific read entry.
#'
#' @import data.table
#' @import GenomicAlignments
#' @import Biostrings
#' @import Rsamtools
#' @export

get_CYTB_amplicon <- function(input_bam) {
  
  .datatable.aware <- TRUE
  
  CYTB_CDS <- IRanges::IRanges(start = 1, end = 185)
  CYTB_target_positions <- c(F129L=24, Y132C=27, G143A=38)
  
  full_map <- unlist(lapply(seq_along(CYTB_CDS), function(i) seq(end(CYTB_CDS[i]), start(CYTB_CDS[i]))))
  gen_code <- Biostrings::getGeneticCode("1")
  target_names <- names(CYTB_target_positions)
  
  bam_files <- if(dir.exists(input_bam)) list.files(input_bam, pattern = "\\.bam$", full.names = TRUE) else input_bam
  final_results <- data.table()
  
  for(bam_file in bam_files) {
    sample_name <- tools::file_path_sans_ext(basename(bam_file))
    message("Analyzing sample: ", sample_name)
    
    aln <- readGAlignments(bam_file, param = ScanBamParam(
      what = c("seq", "cigar", "pos"), 
      which = GRanges("CYTB", IRanges(1, 200))
    ))
    
    if ("cigarillo" %in% loadedNamespaces() || requireNamespace("cigarillo", quietly = TRUE)) {
      seqs_aligned <- cigarillo::project_sequences(mcols(aln)$seq, cigar(aln), from = "query", to = "reference", D.letter = "-", N.letter = "N")
    } else {
      seqs_aligned <- sequenceLayer(mcols(aln)$seq, cigar(aln), from = "query", to = "reference", D.letter = "-", N.letter = "N")
    }
    
    res_dt <- data.table(read_idx = seq_along(aln))
    seq_widths <- width(seqs_aligned)
    read_starts <- start(aln)
    
    for (m_name in target_names) {
      pos <- CYTB_target_positions[[m_name]]
      wt_aa <- substr(m_name, 1, 1)          # F
      mut_aa <- substr(m_name, nchar(m_name), nchar(m_name)) # L
      
      idx_glob <- full_map[c((pos*3)-2, (pos*3)-1, (pos*3))]
      l1 <- idx_glob[1] - read_starts + 1; l2 <- idx_glob[2] - read_starts + 1; l3 <- idx_glob[3] - read_starts + 1
      in_range <- (l1 >= 1 & l3 <= seq_widths)
      
      res_vec <- rep("OUT", length(aln))
      valid_idx <- which(in_range)
      
      if (length(valid_idx) > 0) {
        c1 <- subseq(seqs_aligned[valid_idx], l1[valid_idx], l1[valid_idx])
        c2 <- subseq(seqs_aligned[valid_idx], l2[valid_idx], l2[valid_idx])
        c3 <- subseq(seqs_aligned[valid_idx], l3[valid_idx], l3[valid_idx])
        codons <- paste0(c1, c2, c3)
        
        is_del <- grepl("-", codons)
        res_vec[valid_idx[is_del]] <- "Del"
        
        idx_non_del <- valid_idx[!is_del]
        if (length(idx_non_del) > 0) {
          codon_dnas <- complement(DNAStringSet(codons[!is_del]))
          aas <- as.character(translate(codon_dnas, genetic.code = gen_code, if.fuzzy.codon = "solve", no.init.codon = TRUE))
          
          res_vec[idx_non_del] <- ifelse(aas == wt_aa, "wt", aas)
        }
      }
      set(res_dt, j = m_name, value = res_vec)
    }
    
    is_incomplete <- apply(res_dt[, ..target_names], 1, function(x) any(x == "OUT"))
    res_dt[is_incomplete == TRUE, (target_names) := "OUT"]
    
    is_new <- apply(res_dt[, ..target_names], 1, function(x) {
      any(sapply(seq_along(x), function(i) {
        val <- x[i]
        m_name <- target_names[i]
        mut_char <- substr(m_name, nchar(m_name), nchar(m_name))
        !(val == "wt" | val == mut_char | val == "OUT")
      }))
    })
    
    res_dt[, Status := "KNOWN"]
    res_dt[is_new == TRUE, Status := "NEW"]
    res_dt[is_incomplete == TRUE, Status := "Cant classify"]
    
    sample_summary <- res_dt[, .(N = .N), by = c("Status", target_names)]
    sample_summary[, Sample_ID := sample_name]
    
    setcolorder(sample_summary, c("Sample_ID", "Status", "N", target_names))
    final_results <- rbind(final_results, sample_summary)
  }
  
  return(final_results)
}
