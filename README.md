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

# CYP51

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
