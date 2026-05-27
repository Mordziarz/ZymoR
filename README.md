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

# CYP51

CYP51 haplotypes in Zymoseptoria tritici (based on Huf et. al., 2018 and Gaab et al., 2024)

```r
get_CYP51_res <- ZymoR::get_CYP51(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_CYP51_res)
```

| Sample_ID | Haplotype | Status | INDEL_Promoter | L50 | D107 | D134 | V136 | Y137 | N178 | S188 | S208 | S259 | N284 | H303 | A311 | G312 | A379 | I381 | A410 | G412 | Y459 | G460 | Y461 | G476 | V490 | G510 | N513 | S524 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1M3a | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| 1WIDPa | E5 | KNOWN | YES | S | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| 1WIDPe | H6 | KNOWN | YES | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| 25B16 | H5 | KNOWN | NO | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | T | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| 25B17 | J01 | KNOWN | NO | S | wt | G | A | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | T |
| 25B19 | D13 | KNOWN | NO | wt | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| 25B25 | I2 | KNOWN | NO | S | wt | G | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | T |
| 25B30 | H7 | KNOWN | NO | S | wt | G | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| 25B33 | H3 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | H | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| 2m | H6 | KNOWN | NO | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| 3mb | H7 | KNOWN | NO | S | wt | G | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | wt | T |
| DP152 | Unknown | NEW | NO | wt | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | G | V | wt | wt | wt | wt | S | wt | wt | wt | wt | T |
| IPO323 | WildType (A0) | KNOWN | NO | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt |
| IRE30 | C6 | KNOWN | NO | S | wt | wt | A | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | S | wt | wt | wt | wt | wt |
| NT321 | I3 | KNOWN | NO | S | wt | wt | C | wt | wt | N | wt | wt | wt | wt | wt | wt | G | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | T |
| pgp | D13 | KNOWN | NO | wt | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| scaffolds | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| Z1WIDPd | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |
| Z25B36 | D13 | KNOWN | NO | wt | wt | wt | C | wt | wt | wt | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | wt | wt | H | wt | wt | wt | wt | T |
| Z7Ma | F2 | KNOWN | YES | S | wt | wt | wt | wt | wt | N | wt | wt | wt | wt | wt | wt | wt | V | wt | wt | DEL | DEL | wt | wt | wt | wt | K | wt |


# CYTB

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

# SDHD

SDHD mutations associated with reduced sensitivity to fungicides


```r
get_SDHD_res <- ZymoR::get_SDHD(input_path="path_to_your_folder_or_fasta",output_dir = "your_output_folder_name")


print(get_SDHD_res)
```

| Mutation | Amino acid substitution | SDHI resistance phenotype | Notes |
| :--- | :--- | :--- | :--- |
| **R47W** | Arginine &rarr; Tryptophan at codon 47 | Reduced sensitivity / resistance to multiple SDHIs | Detected in field isolates together with SdhC-H152R; associated with strong reduction in sensitivity to penthiopyrad and several other SDHIs. |
| **I50L** | Isoleucine &rarr; Leucine at codon 50 | Reduced SDHI sensitivity | Reported as a less common SDHD substitution associated with altered SDHI sensitivity; generally considered lower impact than major SDHB/C mutations. |


# SDHC

SDHC mutations associated with reduced sensitivity to fungicides

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

# SDHB

SDHB mutations associated with reduced sensitivity to fungicides

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

# MFS1

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