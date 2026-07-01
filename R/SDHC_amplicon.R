#' Analyze SDHC gene mutations from BAM files
#'
#' This function processes BAM files to identify specific mutations in the SDHC gene.
#' It maps genomic coordinates to coding sequences (CDS) by accounting for the 
#' exon structure of the SDHC gene, performs codon translation using the 
#' mitochondrial genetic code (code 4), and categorizes read-level mutations 
#' based on a reference profile.
#'
#' @param input_bam A character string representing the file path to a single 
#'   BAM file or a directory containing multiple BAM files.
#'
#' @return A \code{data.table} containing the following columns:
#'   \item{Sample_ID}{The name of the analyzed sample.}
#'   \item{Status}{Classification status: "KNOWN" (matches reference or known mutation), 
#'   "NEW" (any other amino acid or deletion), or "Cant classify" (if reads fall outside coverage).}
#'   \item{N}{The number of reads assigned to the specific mutation combination and status.}
#'   \item{T79N, W80S, N86S, H152R, V166M}{Columns representing the detected amino acid 
#'   (or "wt", "Del", "OUT") at each target position.}
#'
#' @details
#' The function maps SDHC gene positions (on the minus strand) and performs translation 
#' of codons using the mitochondrial genetic code.
#' - "wt" is assigned if the amino acid matches the reference.
#' - The specific mutation (e.g., "N" for T79N) is labeled as defined in the target profile.
#' - Any other amino acid or deletion ("Del") is classified as "NEW".
#' - If a read does not cover all target positions, it is marked as "OUT" and the status 
#'   is set to "Cant classify" for that specific read entry.
#'
#' @import data.table
#' @import GenomicAlignments
#' @import Biostrings
#' @import Rsamtools
#' @export

get_SDHC_amplicon <- function(input_bam) {
  
  .datatable.aware <- TRUE
  
  SDHC_CDS <- IRanges::IRanges(start = c(656, 502, 1), end = c(695, 583, 439))
  full_map <- unlist(lapply(seq_along(SDHC_CDS), function(i) seq(end(SDHC_CDS[i]), start(SDHC_CDS[i]))))
  
  gen_code <- Biostrings::getGeneticCode("4")
  SDHC_target_positions <- c(T79N=79, W80S=80, N86S=86, H152R=152, V166M=166)
  target_names <- names(SDHC_target_positions)
  
  bam_files <- if(dir.exists(input_bam)) list.files(input_bam, pattern = "\\.bam$", full.names = TRUE) else input_bam
  final_results <- data.table()
  
  for(bam_file in bam_files) {
    sample_name <- tools::file_path_sans_ext(basename(bam_file))
    message("Analyzing sample: ", sample_name)
    
    aln <- readGAlignments(bam_file, param = ScanBamParam(
      what = c("seq", "cigar", "pos"), 
      which = GRanges("SDHC", IRanges(1, 700))
    ))
    
    if (length(aln) == 0) next
    
    if ("cigarillo" %in% loadedNamespaces() || requireNamespace("cigarillo", quietly = TRUE)) {
      seqs_aligned <- cigarillo::project_sequences(mcols(aln)$seq, cigar(aln), from = "query", to = "reference", D.letter = "-", N.letter = "N")
    } else {
      seqs_aligned <- sequenceLayer(mcols(aln)$seq, cigar(aln), from = "query", to = "reference", D.letter = "-", N.letter = "N")
    }
    
    res_dt <- data.table(read_id = names(aln))
    seq_widths <- width(seqs_aligned)
    read_starts <- start(aln)
    
    for (m_name in target_names) {
      pos <- SDHC_target_positions[[m_name]]
      idx_glob <- full_map[c((pos*3)-2, (pos*3)-1, (pos*3))]
      
      l1 <- idx_glob[1] - read_starts + 1; l2 <- idx_glob[2] - read_starts + 1; l3 <- idx_glob[3] - read_starts + 1
      in_range <- (l1 >= 1 & l2 >= 1 & l3 >= 1 & l1 <= seq_widths & l2 <= seq_widths & l3 <= seq_widths)
      
      res_vec <- rep("OUT", length(aln))
      valid_idx <- which(in_range)
      
      if (length(valid_idx) > 0) {
        c1 <- subseq(seqs_aligned[valid_idx], l1[valid_idx], l1[valid_idx])
        c2 <- subseq(seqs_aligned[valid_idx], l2[valid_idx], l2[valid_idx])
        c3 <- subseq(seqs_aligned[valid_idx], l3[valid_idx], l3[valid_idx])
        codons <- paste0(c1, c2, c3)
        
        is_del <- grepl("-", codons)
        res_vec[valid_idx[is_del]] <- "Del"
        
        if (any(!is_del)) {
          idx_non_del <- valid_idx[!is_del]

          codon_dnas <- complement(DNAStringSet(codons[!is_del]))
          aas <- as.character(translate(codon_dnas, genetic.code = gen_code, if.fuzzy.codon = "solve", no.init.codon = TRUE))
          
          wt_aa <- substr(m_name, 1, 1)
          res_vec[idx_non_del] <- ifelse(aas == wt_aa, "wt", aas)
        }
      }
      set(res_dt, j = m_name, value = res_vec)
    }
    
    is_incomplete <- apply(res_dt[, ..target_names], 1, function(x) any(x == "OUT"))
    res_dt[, Status := "KNOWN"]
    res_dt[is_incomplete == TRUE, (target_names) := "OUT"]
    res_dt[is_incomplete == TRUE, Status := "Cant classify"]
    
    is_new <- apply(res_dt[, ..target_names], 1, function(x) {
      any(sapply(seq_along(x), function(i) {
        val <- x[i]; mut_char <- substr(target_names[i], nchar(target_names[i]), nchar(target_names[i]))
        !(val == "wt" | val == mut_char | val == "OUT")
      }))
    })
    res_dt[is_new == TRUE, Status := "NEW"]
    
    sample_summary <- res_dt[, .(N = .N), by = c("Status", target_names)]
    sample_summary[, Sample_ID := sample_name]
    setcolorder(sample_summary, c("Sample_ID", "Status", "N", target_names))
    final_results <- rbind(final_results, sample_summary)
  }
  
  return(final_results)
}
