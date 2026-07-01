#' Analyze SDHB gene mutations from BAM files
#'
#' This function processes BAM files to identify specific mutations in the SDHB gene.
#' It maps genomic coordinates to coding sequences (CDS) by accounting for exon 
#' structure, performs codon translation using the mitochondrial genetic code (code 4),
#' and categorizes read-level mutations based on a reference profile.
#'
#' @param input_bam A character string representing the file path to a single 
#'   BAM file or a directory containing multiple BAM files.
#'
#' @return A \code{data.table} containing:
#'   \item{Sample_ID}{The name of the analyzed sample.}
#'   \item{Status}{Classification status: "KNOWN" (matches reference or defined variants), 
#'   "NEW" (unexpected amino acid or deletion), or "Cant classify" (insufficient read coverage).}
#'   \item{N}{The number of reads assigned to the specific mutation combination.}
#'   \item{N225I/T, H267Y/R/L, T268I, I269V/T}{Columns containing the detected amino acid 
#'   code ("wt" for wild-type, specific variant letters, "Del" for deletion, or "OUT" for missing coverage).}
#'
#' @details
#' The function uses a mapping strategy to bridge intron gaps based on the SDHB exon 
#' structure (Exon 1: 1-415, Exon 2: 471-622, Exon 3: 682-1008). 
#' 
#' \itemize{
#'   \item \strong{Translation:} Uses the mitochondrial genetic code (NCBI table 4).
#'   \item \strong{Mapping:} Automatically handles genomic BAM alignments using 
#'   \code{cigarillo} (preferred) or \code{sequenceLayer}.
#'   \item \strong{Classification:} 
#'     \itemize{
#'       \item \emph{wt}: Matches the first amino acid defined in the column name.
#'       \item \emph{KNOWN}: Matches any other variant specified in the column name.
#'       \item \emph{NEW}: Any other amino acid or deletion not defined in the reference.
#'     }
#' }
#' 
#' @note 
#' This function assumes that the BAM files are aligned to the SDHB gene region. 
#' The intron structure is handled via a genomic-to-CDS coordinate translation.
#'
#' @import data.table
#' @import GenomicAlignments
#' @import Biostrings
#' @import Rsamtools
#' @export

get_SDHB_amplicon <- function(input_bam) {
  
  .datatable.aware <- TRUE
  
  exons <- list(c(1, 415), c(471, 622), c(682, 1008))
  
  get_genomic_positions <- function(aa_pos) {
    target_nt <- (aa_pos - 1) * 3 + c(1, 2, 3)
    res <- c()
    current_offset <- 0
    
    for (exon in exons) {
      len <- exon[2] - exon[1] + 1
      for (nt in target_nt) {
        if (nt > current_offset && nt <= current_offset + len) {
          res <- c(res, exon[1] + (nt - current_offset - 1))
        }
      }
      current_offset <- current_offset + len
    }
    return(res)
  }
  
  SDHB_target_positions <- c(`N225I/T`=225, `H267Y/R/L`=267, `T268I`=268, `I269V/T`=269)
  target_names <- names(SDHB_target_positions)
  gen_code <- Biostrings::getGeneticCode("4")
  
  bam_files <- if(dir.exists(input_bam)) list.files(input_bam, pattern = "\\.bam$", full.names = TRUE) else input_bam
  final_results <- data.table()
  
  for(bam_file in bam_files) {
    sample_name <- tools::file_path_sans_ext(basename(bam_file))
    message("Analyzing sample: ", sample_name)
    
    aln <- readGAlignments(bam_file, param = ScanBamParam(
      what = c("seq", "cigar", "pos"), 
      which = GRanges("SDHB", IRanges(1, 1008))
    ))
    
    if ("cigarillo" %in% loadedNamespaces() || requireNamespace("cigarillo", quietly = TRUE)) {
      seqs_aligned <- cigarillo::project_sequences(mcols(aln)$seq, cigar(aln), from = "query", to = "reference", D.letter = "-", N.letter = "N")
    } else {
      seqs_aligned <- sequenceLayer(mcols(aln)$seq, cigar(aln), from = "query", to = "reference", D.letter = "-", N.letter = "N")
    }
    
    res_dt <- data.table(read_idx = seq_along(aln))
    read_starts <- start(aln)
    
    for (m_name in target_names) {
      idx_glob <- get_genomic_positions(SDHB_target_positions[[m_name]])
      
      l1 <- idx_glob[1] - read_starts + 1
      l2 <- idx_glob[2] - read_starts + 1
      l3 <- idx_glob[3] - read_starts + 1
      in_range <- (l1 >= 1 & l3 <= width(seqs_aligned))
      
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
          aas <- as.character(translate(DNAStringSet(codons[!is_del]), genetic.code = gen_code, if.fuzzy.codon = "solve", no.init.codon = TRUE))
          
          raw_chars <- gsub("[^A-Z]", "", m_name)
          wt_aa <- substr(raw_chars, 1, 1)
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
        raw_chars <- gsub("[^A-Z]", "", target_names[i])
        known_aas <- unlist(strsplit(raw_chars, ""))
        !(val == "wt" | val %in% known_aas | val == "OUT")
      }))
    })
    
    res_dt[, Status := "KNOWN"]
    res_dt[is_new == TRUE, Status := "NEW"]
    res_dt[is_incomplete == TRUE, Status := "Cant classify"]
    
    sample_summary <- res_dt[, .(N = .N), by = c("Status", target_names)]
    sample_summary[, Sample_ID := sample_name]
    final_results <- rbind(final_results, sample_summary)
  }
  return(final_results)
}
