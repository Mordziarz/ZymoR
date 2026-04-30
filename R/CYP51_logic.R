#' Generate Mutation Matrix and Haplotype Identification for CYP51
#'
#' Processes a single FASTA file or a directory of FASTA files to detect mutations 
#' at specific amino acid positions within the CYP51 gene, maps them to a reference genome, 
#' and identifies the resulting haplotype based on a provided database.
#'
#' @param input_path String. Path to a single file (.fasta, .fna, .txt) or a folder 
#'   containing sequence files.
#' @param output_dir String. Where to save results and tables. Defaults to "zymor_results".
#' @param ... Additional arguments passed to analyze_genome (e.g., max_mismatch).
#'
#' @return A data frame (matrix) where rows are samples and columns include 
#'   Sample_ID, Haplotype name, Status (KNOWN/NEW), and individual columns 
#'   for each analyzed CYP51 mutation position.
#' @export
#'

get_CYP51 <- function(input_path, output_dir = "zymor_results", ...) {

  CYP51_F <- "ACAGGATGTCGTCTGGATAGTTAAGGTA"
  CYP51_R <- "TCTTCCCCAGACCAAATTATCGTTCGA"

  reference_seq <- Biostrings::DNAString("TCAGTTCTTCTCCTCCTTCTCCTCCCTCCTCTCCCACTTTACTACTGCCGGCGACAGCGGCCGGCTGAACAAACTGCTGTAATCCGTACCCACCACGTTGTCGCTGCCATCCACATTGTAAAACTTGAAATCGCGAACCATCGTCGCTGTAATGGTCTGCAATTGCACGTACGCGAATTGCTCGCCGATACATCTGTGTCGTCCCGCACCAAAGGGTAAGTATGGCGACGCCGCGCCCTTGCTTACCAGGCCGTAGCCATAGTCTTCTTTCTCCTCGGCGATGCTTCCTAGGGCAGTCGTCGGGGACAGGTGCTTGTATTTCTCGGACGGGCTCTCGTCCCATCGATGCGGCTCCCAATGGAGGCAGTCGGGAAAGTGCTCGTCCATGCGGCTCGTTGTGCCCGGAGCGGCCAGAAGAGTGTGGGTGGTTGGAATGACGTATGCCGTACCTTCGATGGGCATGGGAGACTTGACCTTGCGCAGAATGGAGTGGATTGGAGCGTGAATACGAAGGGTTTCTTTGACGACTTGATTGAGGAGGGTGAGTTTCGAGAGGTTGGCGTATGTGAGCTCCTTGATACTGCCGTCGGCGTTCACACCGAGCATATCCTTTTGTTCTTGGAGGAGTTCGTCTTGGATGTCGGGGCGGGATGCGAGGCGGAGAGTGATCCAGGACTCGGTCGCAGATGAAGAGTGCTGGCCGGCCATGAGCAGCGCAATCATCATATGAGCAATCTCCTTGTCGGGAATGGCATTGCCGTCCTTGTATTTGCACCTAGAAGTCGAGGGTCAGCATCGTCTCCAGGGAGGAACTTGGGCCGCATGGCCGTGGGACACGCGTGGCAACACTGTATCTGCACGCCTTCTGCTTGTGACTTTACTAACTGCATCAAGTTGTGGATCACTGATGGTACAAGACGTTAGCAAATGCTCGGCGATGTCATGGCCACTAACAGGCCAACTGCAAGACTCCGAGAGCTGGACGCCTCCTGTGCCTGACTTCACTTACTGTCTTCCTCATGTTCGCCCGTTTTGGACTCTCGTCTCTTCTGAATGATCGACATGTATGTCTCGGACATCTTCTTCTGCGCATAATCGCGGCGTCGGTTCTGGGGAAGGGGCGCCCACGGAAGCATAAAGTTGATCGGTGTGAATCCCATATCGAGGTAGTGGTAGAGGTCCGCGAAAGACGAGTCGAAGCCCTCGCGGACTTCCTTTCCTTGCAATGATCGGCTGGCAGTATAGATCGTAAGTTCGGCGAGGGCTGGTGGGAGATCGATCGTGCCGCTGGTCGATGCGAACTTCTTATGAGGGTTGTTGCGGTCGAAGAACTGGCGGGTCTCGGCGGCGATCAAGGTCACGTAGGACTGGAGGGCAGAGGTTGTGAGGCCGTACTTGACGAACTAAGTGACAAGTCAATATAGCTGCACTTTTGCCCGTATGTTCTCAATTCTCCGTACCTTCTTCTGCTCCATGAGCTTCGAATTGGGACAATCATAAACCACATCCTTGCCAAAGACAGGAGTGGTCAGCGGGCTGTATATCTCCTCCGCGTTGACGTCCTTCAGTTTTCCATTCAAAATAAAATCATTGCCCTTGGTGCCCAAGCACACCGTCGTCTTCTTTCCCAGCAGGATGAATGTAAAGACATCTCCATACTTTTCCCGACAGGAGAAGAAGAACTTGTATGGGTCGATGCCGTAGGTGATGGTGCTTCCGATGAAGGGCACCCAGTGGAATACGAGTGGCGGATCGGACAACTTGCCACGGAAGAGAAGTTGTGAGAGGACATTGAGGAGGATGGCGAGCGTGCTGAAGGCGAGGAATCCAAGTCCGACAAGTTTCCAGAGGCTGGTCTGGCCGAATTGCGCGTCGAACTGCGCGAGGACTTCCTGGAGGAGACCCAT")
  
  cds_ranges <- IRanges::IRanges(
    start = c(1461, 1010, 886, 1),
    end   = c(1907, 1403, 904, 775),
    names = c("CDS1", "CDS2", "CDS3", "CDS4")
  )

  CYP51_target_positions <- c(
    L50=50, D107=107, D134=134, V136=136, Y137=137, N178=178, S188=188, S208=208,
    S259=259, N284=284, H303=303, A311=311, G312=312, A379=379, I381=381, A410=410,
    G412=412, Y459=459, G460=460, Y461=461, G476=476, V490=490, G510=510, N513=513, S524=524
  )

  CYP51_db <- list(
    "A1"  = c(G460 = "D"), "A2"  = c(L50 = "S"), "A3"  = c(Y137 = "F"), "A4"  = c(V136 = "C"),
    "A5"  = c(Y461 = "Del"), "A6"  = c(S259 = "F"), "A7"  = c(G476 = "S"), "A8"  = c(Y459 = "C"),
    "A9"  = c(Y459 = "D"), "A10" = c(D107 = "V"), "A11" = c(Y461 = "H"), "B1"  = c(L50 = "S", Y137 = "F"),
    "B2"  = c(V136 = "A", Y461 = "S"), "B3"  = c(I381 = "V", Y459 = "D"), "B4"  = c(I381 = "V", Y459 = "S"),
    "B5"  = c(L50 = "S", Y459 = "D"), "B6"  = c(L50 = "S", G460 = "D"), "B7"  = c(L50 = "S", S188 = "N"),
    "B8"  = c(L50 = "S", Y461 = "S"), "B9"  = c(L50 = "S", Y461 = "H"), "B10" = c(V136 = "C", Y461 = "H"),
    "B11" = c(Y137 = "F", S524 = "T"), "B12" = c(I381 = "V", Y461 = "H"), "B13" = c(L50 = "S", Y459 = "C"),
    "C1"  = c(I381 = "V", Y461 = "H", G510 = "C"), "C2"  = c(L50 = "S", I381 = "V", Y459 = "S"),
    "C3"  = c(L50 = "S", I381 = "V", Y459 = "C"), "C4"  = c(L50 = "S", I381 = "V", Y459 = "D"),
    "C5"  = c(L50 = "S", V136 = "C", Y461 = "S"), "C6"  = c(L50 = "S", V136 = "A", Y461 = "S"),
    "C7"  = c(L50 = "S", V136 = "A", Y461 = "H"), "C8"  = c(L50 = "S", I381 = "V", Y461 = "H"),
    "C9"  = c(L50 = "S", V136 = "A", Y461 = "N"), "C10" = c(V136 = "A", Y459 = "Del", G460 = "Del"),
    "C11" = c(V136 = "A", Y461 = "S", S524 = "T"), "C12" = c(V136 = "A", I381 = "V", Y461 = "H"),
    "C13" = c(V136 = "A", S188 = "N", Y459 = "C"), "C14" = c(L50 = "S", G312 = "A", Y459 = "D"),
    "C15" = c(L50 = "S", N178 = "S", S188 = "N"), "C16" = c(D107 = "V", I381 = "V", Y461 = "H"),
    "C17" = c(V136 = "C", I381 = "V", Y461 = "H"), "C18" = c(L50 = "S", S188 = "N", N513 = "K"),
    "C19" = c(L50 = "S", S188 = "N", G460 = "D"), "C20" = c(L50 = "S", Y461 = "S", S524 = "T"),
    "C21" = c(L50 = "S", V136 = "G", Y461 = "S"), "D1"  = c(S188 = "N", Y459 = "Del", G460 = "Del", N513 = "K"),
    "D2"  = c(L50 = "S", S188 = "N", G460 = "D", N513 = "K"), "D3"  = c(L50 = "S", S188 = "N", Y461 = "H", N513 = "K"),
    "D4"  = c(L50 = "S", S188 = "N", Y459 = "Del", G460 = "Del"), "D5"  = c(L50 = "S", Y459 = "Del", G460 = "Del", N513 = "K"),
    "D6"  = c(L50 = "S", D134 = "G", V136 = "G", Y461 = "S"), "D7"  = c(L50 = "S", V136 = "A", Y461 = "S", S524 = "T"),
    "D8"  = c(V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del"), "D9"  = c(L50 = "S", D134 = "G", V136 = "A", Y461 = "H"),
    "D10" = c(L50 = "S", V136 = "A", I381 = "V", Y461 = "H"), "D11" = c(L50 = "S", V136 = "A", Y461 = "H", S524 = "T"),
    "D12" = c(V136 = "A", I381 = "V", Y459 = "Del", G460 = "Del"), "D13" = c(V136 = "C", I381 = "V", Y461 = "H", S524 = "T"),
    "D14" = c(L50 = "S", V136 = "A", G460 = "Del", S524 = "T"), "D15" = c(D134 = "G", V136 = "A", S188 = "N", Y461 = "L"),
    "D16" = c(L50 = "S", I381 = "V", G460 = "Del", Y461 = "H"), "D17" = c(D134 = "G", S188 = "N", Y459 = "Del", G460 = "Del"),
    "D18" = c(D134 = "G", V136 = "A", I381 = "V", Y461 = "H"), "D19" = c(D107 = "V", I381 = "V", N513 = "K", S524 = "T"),
    "D20" = c(L50 = "S", G312 = "A", I381 = "V", Y461 = "H"), "D21" = c(L50 = "S", S188 = "N", I381 = "V", Y461 = "H"),
    "D22" = c(L50 = "S", G312 = "A", I381 = "V", Y459 = "D"), "D23" = c(L50 = "S", S188 = "N", Y459 = "C", N513 = "K"),
    "D24" = c(L50 = "S", A311 = "G", Y461 = "S", V490 = "L"), "D25" = c(L50 = "S", I381 = "V", Y461 = "H", S524 = "T"),
    "E1"  = c(L50 = "S", S188 = "N", Y459 = "Del", G460 = "Del", N513 = "K"), "E2"  = c(L50 = "S", D134 = "G", V136 = "A", Y461 = "S", S524 = "T"),
    "E3"  = c(L50 = "S", V136 = "A", I381 = "V", Y461 = "S", S524 = "T"), "E4"  = c(L50 = "S", D134 = "G", V136 = "A", I381 = "V", Y461 = "H"),
    "E5"  = c(L50 = "S", V136 = "A", I381 = "V", Y461 = "H", S524 = "T"), "E6"  = c(L50 = "S", V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del"),
    "E7"  = c(L50 = "S", V136 = "C", I381 = "V", Y461 = "H", S524 = "T"), "E8"  = c(L50 = "S", V136 = "C", S188 = "N", I381 = "V", Y461 = "H"),
    "E9"  = c(L50 = "S", D134 = "G", V136 = "A", I381 = "V", Y459 = "S"), "E10" = c(L50 = "S", D134 = "G", V136 = "A", I381 = "V", Y461 = "S"),
    "E11" = c(V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del", N513 = "K"), "E12" = c(L50 = "S", S188 = "N", A379 = "G", Y459 = "Del", G460 = "Del"),
    "E13" = c(L50 = "S", D134 = "G", V136 = "G", Y461 = "S", S524 = "T"), "E14" = c(S188 = "N", Y459 = "Del", G460 = "Del", G510 = "C", N513 = "K"),
    "E15" = c(L50 = "S", V136 = "A", Y459 = "Del", G460 = "Del", N513 = "K"), "E16" = c(L50 = "S", V136 = "C", Y459 = "Del", G460 = "Del", N513 = "K"),
    "E17" = c(L50 = "S", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"), "E18" = c(L50 = "S", S188 = "N", I381 = "V", G460 = "Del", Y461 = "H"),
    "E19" = c(L50 = "S", D107 = "V", I381 = "V", Y461 = "H", S524 = "T"), "E20" = c(L50 = "S", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del"),
    "E21" = c(D134 = "G", V136 = "A", S188 = "N", G460 = "Del", Y461 = "Del"), "E22" = c(S188 = "N", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"),
    "E23" = c(V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del", S524 = "T"), "F1"  = c(L50 = "S", S188 = "N", A379 = "G", I381 = "V", Y459 = "D", S524 = "T"),
    "F2"  = c(L50 = "S", S188 = "N", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"), "F3"  = c(L50 = "S", V136 = "C", S188 = "N", Y459 = "Del", G460 = "Del", N513 = "K"),
    "F4"  = c(L50 = "S", V136 = "C", S188 = "N", I381 = "V", Y461 = "H", S524 = "T"), "F5"  = c(L50 = "S", V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del", N513 = "K"),
    "F6"  = c(L50 = "S", V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del", S524 = "T"), "F7"  = c(L50 = "S", V136 = "A", A379 = "G", I381 = "V", Y461 = "S", S524 = "T"),
    "F8"  = c(L50 = "S", D134 = "G", V136 = "A", I381 = "V", Y461 = "H", S524 = "T"), "F9"  = c(V136 = "C", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", S524 = "T"),
    "F10" = c(L50 = "S", D134 = "G", V136 = "A", I381 = "V", Y461 = "S", S524 = "T"), "F11" = c(L50 = "S", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"),
    "F12" = c(L50 = "S", D134 = "G", V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del"), "F13" = c(V136 = "C", S188 = "N", Y459 = "Del", G460 = "Del", G510 = "C", N513 = "K"),
    "F14" = c(S188 = "N", I381 = "V", Y459 = "Del", G460 = "Del", Y461 = "G", N513 = "K"), "F15" = c(L50 = "S", V136 = "A", S188 = "N", A379 = "G", Y459 = "Del", G460 = "Del"),
    "F16" = c(S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"), "G1"  = c(L50 = "S", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"),
    "G2"  = c(L50 = "S", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", S524 = "T"), "G3"  = c(L50 = "S", V136 = "A", S188 = "N", Y459 = "Del", G460 = "Del", N513 = "K", S524 = "T"),
    "G4"  = c(L50 = "S", S188 = "N", S208 = "T", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del"), "G5"  = c(L50 = "S", D134 = "G", V136 = "A", S208 = "T", I381 = "V", Y461 = "H", S524 = "T"),
    "G6"  = c(L50 = "S", D134 = "G", V136 = "A", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"), "G7"  = c(L50 = "S", V136 = "A", S188 = "N", A379 = "G", I381 = "V", Y461 = "S", S524 = "T"),
    "G8"  = c(S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", Y461 = "G", N513 = "K"), "G9"  = c(L50 = "S", V136 = "A", S188 = "N", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"),
    "G10" = c(L50 = "S", D134 = "G", V136 = "A", A379 = "G", I381 = "V", Y461 = "S", S524 = "T"), "H1"  = c(L50 = "S", S188 = "N", A379 = "G", I381 = "V", G412 = "A", Y459 = "Del", G460 = "Del", N513 = "K"),
    "H2"  = c(L50 = "S", S188 = "N", S208 = "T", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"), "H3"  = c(L50 = "S", S188 = "N", N284 = "H", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"),
    "H4"  = c(L50 = "S", V136 = "A", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", S524 = "T"), "H5"  = c(L50 = "S", S188 = "N", A379 = "G", I381 = "V", A410 = "T", Y459 = "Del", G460 = "Del", N513 = "K"),
    "H6"  = c(L50 = "S", V136 = "C", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", S524 = "T"), "H7"  = c(L50 = "S", D134 = "G", V136 = "A", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", S524 = "T"),
    "H8"  = c(L50 = "S", S188 = "N", H303 = "Y", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K"), "H9"  = c(L50 = "S", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K", S524 = "T"),
    "I1"  = c(L50 = "S", V136 = "A", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K", S524 = "T"), "I2"  = c(L50 = "S", D134 = "G", V136 = "A", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K", S524 = "T"),
    "I3"  = c(L50 = "S", V136 = "C", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K", S524 = "T"), "J01" = c(L50 = "S", D134 = "G", V136 = "A", S188 = "N", A379 = "G", I381 = "V", Y459 = "Del", G460 = "Del", N513 = "K", S524 = "T")
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

  cat("=== Starting get_CYP51 Pipeline ===\n")
  cat("Total files to process:", length(fasta_files), "\n\n")

  for (f_path in fasta_files) {
    file_label <- tools::file_path_sans_ext(basename(f_path))

    sample_folder <- file.path(output_dir, file_label)
    if (!dir.exists(sample_folder)) dir.create(sample_folder, recursive = TRUE)

    cat("Processing:", file_label, "... ")

    amp_results <- analyze_genome(
      genome_path = f_path,
      forward_primers = c(CYP51_F),
      reverse_primers = c(CYP51_R),
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

    mut_matrix <- .internal_cyp_processing(
      analyze_genome_results = amp_results,
      reference_seq = reference_seq,
      cds_ranges = cds_ranges,
      CYP51_target_positions = CYP51_target_positions,
      CYP51_db = CYP51_db
    )

    write.csv(mut_matrix, file.path(sample_folder, paste0(file_label, "_mutations.csv")), row.names = FALSE)
    all_mutation_results[[file_label]] <- mut_matrix
  }

  cat("\n--- Finalizing Global Results ---\n")

  if(length(global_amplicons_with) > 0) {
    all_with <- Biostrings::DNAStringSet(global_amplicons_with)
    Biostrings::writeXStringSet(all_with, file.path(output_dir, "CYP51_with_primers.fasta"))
    cat("- Saved: CYP51_with_primers.fasta\n")
  }

  if(length(global_amplicons_no) > 0) {
    all_no <- Biostrings::DNAStringSet(global_amplicons_no)
    Biostrings::writeXStringSet(all_no, file.path(output_dir, "CYP51_without_primers.fasta"))
    cat("- Saved: CYP51_without_primers.fasta\n")
  }

  if (length(all_mutation_results) > 0) {
    final_table <- do.call(rbind, all_mutation_results)
    rownames(final_table) <- NULL
    write.csv(final_table, file.path(output_dir, "CYP51_haplotype_table.csv"), row.names = FALSE)
    cat("- Saved: CYP51_haplotype_table.csv\n")
  }

  if (length(not_found_list) > 0) {
    writeLines(not_found_list, file.path(output_dir, "CYP51_missing_amplicons.txt"))
    cat("- Warning: Gene not found in", length(not_found_list), "samples (see CYP51_missing_amplicons.txt)\n")
  }

  cat("\n=== Complete ===\n")
  return(final_table)
}

.internal_cyp_processing <- function(analyze_genome_results, reference_seq, cds_ranges, CYP51_target_positions, CYP51_db) {
  
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

  mut_names <- names(CYP51_target_positions)
  col_names <- c("Sample_ID", "Haplotype", "Status", "INDEL_Promoter", mut_names)

  res_table <- data.frame(matrix(NA, nrow = length(analyze_genome_results), ncol = length(col_names)))
  colnames(res_table) <- col_names

  for (i in seq_along(analyze_genome_results)) {
    amp_data <- analyze_genome_results[[i]]
    curr_amp <- amp_data[["with_p"]]
    curr_id  <- as.character(amp_data[["amplicon_id"]])

    res_table[i, "Sample_ID"] <- curr_id
    detected_muts <- list()

    aln_for_indel <- pwalign::pairwiseAlignment(
      pattern = DNAString(as.character(curr_amp)), 
      subject = reference_seq, 
      type = "global-local"
    )
    end_in_amp <- Biostrings::end(pwalign::pattern(aln_for_indel))
    amp_total_len <- nchar(as.character(curr_amp))
    trailing_seq_len <- amp_total_len - end_in_amp

    if (trailing_seq_len > 150) {
      res_table[i, "INDEL_Promoter"] <- "YES"
    } else {
      res_table[i, "INDEL_Promoter"] <- "NO"
    }

    for (m_name in mut_names) {
      aa_num <- CYP51_target_positions[[m_name]]
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

    h_name <- .internal_identify_haplotype(final_detected, CYP51_db)
    res_table[i, "Haplotype"] <- h_name

    is_known <- h_name %in% names(CYP51_db) || h_name == "WildType (A0)"
    res_table[i, "Status"] <- if (is_known) "KNOWN" else "NEW"
  }

  return(res_table)
}

.internal_identify_haplotype <- function(detected_muts, haplotype_db) {
  detected_muts <- detected_muts[detected_muts != "wt" & !is.na(detected_muts)]
  
  if (length(detected_muts) == 0) return("WildType (A0)")
  
  det_names <- names(detected_muts)
  det_vals  <- toupper(as.character(detected_muts))
  
  for (h_name in names(haplotype_db)) {
    db_muts <- haplotype_db[[h_name]]
    if (length(det_vals) != length(db_muts)) next
    
    db_names <- names(db_muts)
    db_vals  <- toupper(as.character(db_muts))
    
    if (all(det_names[order(det_names)] == db_names[order(db_names)]) &&
        all(det_vals[order(det_names)] == db_vals[order(det_names)])) {
      return(h_name)
    }
  }
  return("Unknown")
}