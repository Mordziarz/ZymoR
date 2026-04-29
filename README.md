# ZymoR
ZymoR is a package that enables rapid identification of haplotypes in Zymoseptoria.


# Installation 
```r
install.packages("devtools")
library(devtools)
devtools::install_github('Mordziarz/ZymoR')
library(ZymoR)
```

## Install BiocManager if not installed

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```
## Install missing Bioconductor packages

```r
BiocManager::install(setdiff(c("Biostrings","IRanges","pwalign"), 
                             installed.packages()[,"Package"]))
```

## Load Bioconductor packages

```r
lapply(c("Biostrings","IRanges","pwalign"), 
       library, character.only = TRUE)
```

## Packages
```r
library(Biostrings)
library(IRanges)
library(pwalign)
```

# Find gene

The first step is to locate and extract the CYP51 gene from your genomic data (FASTA files) using pre-defined primers.

| Gene | Forward primer name | Reverse primer name |
| :--- | :--- | :--- |
| **CYP51** | CYP51_F | CYP51_R |

```r
analyze_genome_results <- analyze_genome("path_to_your_fasta",
                                          forward_primers = CYP51_F,
                                          reverse_primers = CYP51_R,
                                          all = T,
                                          max_mismatch = 1,
                                          output_dir = "contigs_results")

analyze_genome_results$`your_fasta_name` # If you want to view the sequences in R 
```
You can also:
```r
# Define your primers
Primer_F <- "ATGGGTCTCCTCACTGAAATC" # replace with your sequence
Primer_R <- "TTAGTTTCTCCCTTGGTACGC" # replace with your sequence

AND:

analyze_genome_results <- analyze_genome("path_to_your_fasta",
                                          forward_primers = Primer_F,
                                          reverse_primers = Primer_R,
                                          all = T,
                                          max_mismatch = 1,
                                          output_dir = "contigs_results")
```

# Compare to database

```r
get_mutation_matrix_results <- get_mutation_matrix(analyze_genome_results = analyze_genome_results,
                                                 reference_seq = CYP51_reference,
                                                 target_positions = target_positions_CYP51,
                                                 cds_ranges = CYP51_CDS,
                                                 haplotype_db = CYP51_db
)

print(get_mutation_matrix_results)
```


| Sample_ID | Haplotype | Status | L50 | D107 | D134 | V136 | Y137 | N178 | S188 | S208 | S259 | N284 | H303 | A311 | G312 | A379 | I381 | A410 | G412 | Y459 | G460 | Y461 | G476 | V490 | G510 | N513 | S524 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1WIDPa** | E5 | KNOWN | S | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |

# Batch analysis of multiple FASTA files

ZymoR enables the analysis of multiple FASTA files at once:

```r
run_zymor_pipeline_results <- run_zymor_pipeline("/dane/Septorie/",
                                                 forward_primers = CYP51_F,
                                                 reverse_primers = CYP51_R,
                                                 all = T,
                                                 max_mismatch = 1,
                                                 output_dir = "zymo_test",
                                                 reference_seq = CYP51_reference,
                                                 target_positions = CYP51_target_positions,
                                                 cds_ranges = CYP51_CDS,
                                                 haplotype_db = CYP51_db)


print(run_zymor_pipeline_results)
```

| Sample_ID | Haplotype | Status | L50 | D107 | D134 | V136 | Y137 | N178 | S188 | S208 | S259 | N284 | H303 | A311 | G312 | A379 | I381 | A410 | G412 | Y459 | G460 | Y461 | G476 | V490 | G510 | N513 | S524 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1M3a** | Unknown | **NEW** | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | K | wt |
| **1WIDPa** | E5 | **KNOWN** | S | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| **1WIDPe** | Unknown | **NEW** | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | wt | T |
| **2m** | Unknown | **NEW** | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | wt | T |
| **3mb** | Unknown | **NEW** | S | wt | G | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | wt | T |
| **DP152** | Unknown | **NEW** | M | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | wt | wt | S | wt | wt | wt | wt | T |
| **IPO323** | Unknown | **NEW** | M | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt |
| **IRE30** | C6 | **KNOWN** | S | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | S | wt | wt | wt | wt | wt |
| **NT321** | Unknown | **NEW** | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | K | T |
| **pgp** | Unknown | **NEW** | M | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| **scaffolds**| Unknown | **NEW** | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | K | wt |
| **Z1WIDPd** | Unknown | **NEW** | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | K | wt |
| **Z25B36** | Unknown | **NEW** | M | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| **Z7Ma** | Unknown | **NEW** | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | DEL | wt | wt | wt | K | wt |
