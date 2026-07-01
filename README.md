# ZymoR
ZymoR is a package for rapid identification of Cyp51 haplotypes in Zymoseptoria tritici. The package also detects mutations in CYTB and SDH genes associated with fungicide resistance, as well as insertions in the promoter regions of CYP51 and MFS1. 

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
BiocManager::install(setdiff(c("Biostrings","IRanges","pwalign","Rsamtools","GenomicAlignments","Biostrings"), 
                             installed.packages()[,"Package"]))
```

## Load Bioconductor packages

```r
lapply(c("Biostrings","IRanges","pwalign","Rsamtools","GenomicAlignments","Biostrings"), 
       library, character.only = TRUE)
```

## Install data.table

```r
install.packages("data.table")
```

## Packages
```r
library(Biostrings)
library(IRanges)
library(pwalign)
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(data.table)
```
# Identification from genome assembly
## CYP51

CYP51 haplotypes in Zymoseptoria tritici (based on Huf et. al., 2018 and Gaab et al., 2024).

```r
get_CYP51_res <- ZymoR::get_CYP51(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_CYP51_res)
```

| Sample_ID | Haplotype | Status | INDEL_Promoter | L50 | D107 | D134 | V136 | Y137 | N178 | S188 | S208 | S259 | N284 | H303 | A311 | G312 | A379 | I381 | A410 | G412 | Y459 | G460 | Y461 | G476 | V490 | G510 | N513 | S524 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| strain_13 | WildType (A0) | KNOWN | NO | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt |
| strain_14 | C6 | KNOWN | NO | S | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | S | wt | wt | wt | wt | wt |
| strain_6 | D13 | KNOWN | NO | wt | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| strain_16 | D13 | KNOWN | NO | wt | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| strain_19 | D13 | KNOWN | NO | wt | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| strain_2 | E5 | KNOWN | YES | S | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| strain_1 | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| strain_17 | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| strain_18 | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| strain_20 | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| strain_9 | H3 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | H | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| strain_4 | H5 | KNOWN | NO | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | T | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| strain_3 | H6 | KNOWN | YES | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| strain_10 | H6 | KNOWN | NO | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| strain_8 | H7 | KNOWN | NO | S | wt | G | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| strain_11 | H7 | KNOWN | NO | S | wt | G | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| strain_7 | I2 | KNOWN | NO | S | wt | G | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | T |
| strain_15 | I3 | KNOWN | NO | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | T |
| strain_5 | J01 | KNOWN | NO | S | wt | G | A | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | T |
| strain_12 | Unknown | NEW | NO | wt | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | wt | wt | S | wt | wt | wt | wt | T |


## CYTB

Mutations in the mitochondrial cytochrome b (cytb) gene associated with resistance to QoI (strobilurin) fungicides (FRAC Group 11).

```r
get_CYTB_res <- ZymoR::get_CYTB(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_CYTB_res)
```

| Mutation | Amino acid substitution | Fungicide group affected | Resistance level | Status in *Z. tritici* |
| :--- | :--- | :--- | :--- | :--- |
| **F129L** | Phenylalanine &rarr; Leucine (codon 129) | QoIs (strobilurins) | Moderate / partial | Rare or sporadic in field populations (PMC) |
| **Y132C** | Tyrosine &rarr; Cysteine (codon 132) | metyltetraprole (QoL subgroup 11A) | Reduced sensitivity / resistance | Laboratory-selected; associated with fitness penalty (ScienceDirect) |
| **G143A** | Glycine &rarr; Alanine (codon 143) | QoIs (strobilurins) | High / near-complete resistance | Major global field resistance mutation (PMC) |
| **G137R** | Glycine &rarr; Arginine (codon 137) | QoIs (strobilurins) | Moderate / partial | Very rare; mostly reported in other fungi, only sporadic relevance to *Z. tritici* |

## SDHD

SDHD mutations associated with reduced sensitivity to fungicides.


```r
get_SDHD_res <- ZymoR::get_SDHD(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_SDHD_res)
```

| Mutation | Amino acid substitution | SDHI resistance phenotype | Notes |
| :--- | :--- | :--- | :--- |
| **R47W** | Arginine &rarr; Tryptophan at codon 47 | Reduced sensitivity / resistance to multiple SDHIs | Detected in field isolates together with SdhC-H152R; associated with strong reduction in sensitivity to penthiopyrad and several other SDHIs. |
| **I50L** | Isoleucine &rarr; Leucine at codon 50 | Reduced SDHI sensitivity | Reported as a less common SDHD substitution associated with altered SDHI sensitivity; generally considered lower impact than major SDHB/C mutations. |


## SDHC

SDHC mutations associated with reduced sensitivity to fungicides.

```r
get_SDHC_res <- ZymoR::get_SDHC(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_SDHC_res)
```

| Mutation | Amino acid substitution | SDHI resistance phenotype | Notes |
| :--- | :--- | :--- | :--- |
| **T79N/I** | Threonine &rarr; Asparagine / Isoleucine at codon 79 | Low to moderate resistance to multiple SDHIs | Among the most common SDHC mutations in European field populations; associated with reduced sensitivity to fluxapyroxad, bixafen, penthiopyrad, and related SDHIs. |
| **W80S** | Tryptophan &rarr; Serine at codon 80 | Moderate SDHI resistance | First detected in UK populations around 2012; confers reduced sensitivity to several commercial SDHIs. |
| **S83G** | Serine &rarr; Glycine at codon 83 | Reduced SDHI sensitivity (less commonly reported) | Rare mutation; discussed less frequently than T79N or H152R and considered of lower epidemiological importance. |
| **N86S** | Asparagine &rarr; Serine at codon 86 | Low to moderate resistance | One of the most prevalent SDHC mutations in Europe; widespread in western and northern European populations. |
| **H152R** | Histidine &rarr; Arginine at codon 152 | High resistance / major sensitivity shift | The most impactful SDHC mutation identified so far in *Z. tritici*; associated with strong cross-resistance to most commercial SDHIs and reduced field efficacy. |
| **V166M** | Valine &rarr; Methionine at codon 166 | Reduced SDHI sensitivity | Less common field mutation; reported from European monitoring studies, usually with lower resistance impact than H152R. |

## SDHB

SDHB mutations associated with reduced sensitivity to fungicides.

```r
get_SDHB_res <- ZymoR::get_SDHB(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_SDHB_res)
```

| Mutation | Amino acid substitution | SDHI resistance phenotype | Notes |
| :--- | :--- | :--- | :--- |
| **N225I/T** | Asparagine &rarr; Isoleucine / Threonine at codon 225 | Low to moderate resistance to SDHIs | Common field mutation; associated with reduced sensitivity but generally weaker than major SDHC mutations like H152R. Often detected in European monitoring populations. |
| **H267Y/R/L** | Histidine &rarr; Tyrosine / Arginine / Leucine at codon 267 | Moderate resistance (variable cross-resistance) | Lab- and field-documented mutation cluster; can affect sensitivity depending on SDHI compound. Sometimes shows differential cross-resistance patterns. |
| **T268I** | Threonine &rarr; Isoleucine at codon 268 | Low to moderate resistance | Frequently observed in field populations; often co-occurs with other SDH mutations; contributes to additive resistance. |
| **I269V/T** | Isoleucine &rarr; Valine / Threonine at codon 269 | Low resistance / minor sensitivity shift | Less frequent and typically lower impact; usually part of multi-mutation genotypes rather than a strong standalone resistance driver. |

## MFS1

Insertions found in promoter region of MFS1.

```r
get_MFS1_res <- ZymoR::get_MFS1(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_MFS1_res)
```

* **Type I (519 bp Insert)**
  * **Description:** This is a 519-base pair sequence derived from an ancient Ty1/Copia retrotransposon (an LTR retro-element). It contains four conserved MCB (MluI cell cycle box) motifs that act as strong transcriptional activators.
  * **Impact:** It is the strongest known driver of *MFS1* overexpression, leading to high resistance levels across multiple fungicide classes, including DMIs (azoles), SDHIs, and QoIs.

* **Type II (150–369 bp Inserts)**
  * **Description:** These are shorter insertions (variants include IIA, IIB, and others up to IIF) that share conserved domains with retrotransposons of the RLX_Lard_Gridr family.
  * **Impact:** Like Type I, these lead to *MFS1* overexpression and an MDR (multidrug resistance) phenotype, although the level of resistance and expression can be slightly lower than Type I.

* **Type III (149 bp Insert)**
  * **Description:** A significantly shorter insertion of 149 base pairs.
  * **Impact:** This mutation confers weaker resistance compared to Types I and II, as it results in a lower level of *MFS1* overexpression.

* **Other Rare Variants**
  * **Description & Impact:** Includes larger inserts (>519 bp) and specific insertion-deletion (indel) events in the 5&prime; UTR. Some of these, such as Indel VIII, are linked to terbinafine and metyltetraprole resistance, while others (like Type VI) appear not to confer resistance at all, suggesting the promoter is a highly plastic "hot-spot" for mutations.

# Nanopore long-read haplotype identification

The analysis requires amplicon sequencing; so far, we have tested the method using Nanopore sequencing. Below are the primers for the individual genes:

| Gene  | Forward              | Reverse               |
|-------|----------------------|-----------------------|
| Cyp51 | GACCTGCAGGCAGAACTAAG | ACAGGATGTCGTCTGGATAGT |
| CytB  | GCACGTGGGTAGAGGGTTAT | CAGGTGGAGTTTGCATAGGG  |
| SDHB  | ATACCACACAATGGCTCTTCG | GTCTTCCGTCGATTTCGAGAC |
| SDHC  | AGAAGAGGAGAGGGAGAGAAGC | GCACTCCCTTGGGTCCTGAT  |
| SDHD  | CACTCCTCCAAACCGTATCCT | GGCATCATCGTCAAGCAAG   |

## Reference sequence

ZymoR includes a built-in database of reference sequences for the `CYP51`, `CYTB`, `SDHB`, `SDHC`, and `SDHD` genes. You can export these to your working directory for downstream analysis.

To save the reference sequences (FASTA format) to your current working directory, run in R:

```r
get_reference_sequences()

```

## Mapping

Once the reference file (reference_zymoseptoria.fasta) is saved, you can map your sequencing data against it using minimap2, followed by sorting and indexing the resulting BAM files.

You can perform the mapping, sorting, and indexing in a single command line operation in a Linux terminal:

```bash
minimap2 -x map-ont -a -t 8 reference_zymoseptoria.fasta barcode26.fastq.gz | samtools sort -o maps.bam - && samtools index maps.bam

```

## Variant Analysis

ZymoR enables you to analyze the `CYP51`, `CYTB`, `SDHB`, `SDHC`, and `SDHD` genes directly from BAM files. By performing the mapping, sorting, and indexing steps described above, you can effectively assess genetic variation and diversity within a single sample using the aligned reads.

### CYP51

```R
results <- get_CYP51_amplicon(input_bam = "path_to_your_bam_or_directory")

print(results)

write.csv2(results,"get_CYP51_amplicon.csv")
```
The result table generated by the get_CYP51_amplicon() function provides a summary of haplotype counts for each analyzed sample.

| Sample_ID | Haplotype | N   | Status | L50 | D107 | D134 | V136 | Y137 | N178 | S188 | S208 | S259 | N284 | H303 | A311 | G312 | A379 | I381 | A410 | G412 | Y459 | G460 | Y461 | G476 | V490 |
|-----------|-----------|-----|--------|-----|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| barcode46 | F8        | 135 | KNOWN  | S   | wt   | G    | A    | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | V    | wt   | wt   | wt   | wt   | H    | wt   | wt   |
| barcode38 | G10       | 96  | KNOWN  | S   | wt   | G    | A    | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | G    | V    | wt   | wt   | wt   | wt   | S    | wt   | wt   |
| barcode46 | F4        | 59  | KNOWN  | S   | wt   | wt   | C    | wt   | wt   | N    | wt   | wt   | wt   | wt   | wt   | wt   | wt   | V    | wt   | wt   | wt   | wt   | H    | wt   | wt   |
| barcode30 | F8        | 48  | KNOWN  | S   | wt   | G    | A    | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | V    | wt   | wt   | wt   | wt   | H    | wt   | wt   |
| barcode94 | F8        | 46  | KNOWN  | S   | wt   | G    | A    | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | V    | wt   | wt   | wt   | wt   | H    | wt   | wt   |
| barcode86 | E1        | 1   | KNOWN  | S   | wt   | wt   | wt   | wt   | wt   | N    | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | Del  | Del  | wt   | wt   | wt   |
| barcode94 | H7        | 1   | KNOWN  | S   | wt   | G    | A    | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | G    | V    | wt   | wt   | Del  | Del  | wt   | wt   | wt   |
| barcode94 | G2        | 1   | KNOWN  | S   | wt   | wt   | wt   | wt   | wt   | N    | wt   | wt   | wt   | wt   | wt   | wt   | G    | V    | wt   | wt   | Del  | Del  | wt   | wt   | wt   |
| barcode94 | D25       | 1   | KNOWN  | S   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | V    | wt   | wt   | wt   | wt   | H    | wt   | wt   |
| barcode94 | C8        | 1   | KNOWN  | S   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | wt   | V    | wt   | wt   | wt   | wt   | H    | wt   | wt   |

### CYTB

```R
results <- get_CYTB_amplicon(input_bam = "path_to_your_bam_or_directory")

print(results)

write.csv2(results,"get_CYTB_amplicon.csv")
```

| Sample_ID | Status | N | F129L | Y132C | G143A |
| :--- | :--- | :--- | :--- | :--- | :--- |
| barcode26 | KNOWN | 41381 | wt | wt | A |
| barcode50 | KNOWN | 25794 | wt | wt | A |
| barcode74 | KNOWN | 24207 | wt | wt | A |
| barcode82 | KNOWN | 23503 | wt | wt | A |
| barcode90 | KNOWN | 22816 | wt | wt | A |

### SDHB

```R
results <- get_SDHB_amplicon(input_bam = "path_to_your_bam_or_directory")

print(results)

write.csv2(results,"get_SDHB_amplicon.csv")
```

| Sample_ID | Status | N | N225I/T | H267Y/R/L | T268I | I269V/T |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| barcode27 | KNOWN | 10867 | wt | wt | wt | wt |
| barcode67 | KNOWN | 7178 | wt | wt | wt | wt |
| barcode35 | KNOWN | 6687 | wt | wt | wt | wt |
| barcode43 | KNOWN | 6022 | wt | wt | wt | wt |
| barcode91 | KNOWN | 5518 | wt | wt | wt | wt |

### SDHC

```R
results <- get_SDHC_amplicon(input_bam = "path_to_your_bam_or_directory")

print(results)

write.csv2(results,"get_SDHC_amplicon.csv")
```

| Sample_ID | Status | N | T79N | W80S | N86S | H152R | V166M |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| barcode28 | KNOWN | 8065 | wt | wt | wt | wt | wt |
| barcode55 | KNOWN | 7202 | wt | wt | wt | wt | wt |
| barcode28 | KNOWN | 6334 | wt | wt | S | wt | wt |
| barcode44 | KNOWN | 4186 | wt | wt | wt | wt | wt |
| barcode60 | KNOWN | 3824 | wt | wt | wt | wt | wt |


### SDHD


```R
results <- get_SDHD_amplicon(input_bam = "path_to_your_bam_or_directory")

print(results)

write.csv2(results,"get_SDHD_amplicon.csv")
```

| Sample_ID | Status | N | R47W | I50L |
| :--- | :--- | :--- | :--- | :--- |
| barcode37 | KNOWN | 16044 | wt | wt |
| barcode45 | KNOWN | 11047 | wt | wt |
| barcode29 | KNOWN | 10305 | wt | wt |
| barcode53 | KNOWN | 8651 | wt | wt |
| barcode61 | KNOWN | 7625 | wt | wt |