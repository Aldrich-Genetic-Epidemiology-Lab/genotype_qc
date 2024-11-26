# Genotype QC and imputation pipeline
A standardized pipeline for quality control and imputation of array-based genotype data.

## Getting started
### Software
You should have the following command line tools installed:
* ```R``` (v4.1.0)
* ```R data.table``` (v1.16.2)
* ```R ggplot2``` (v3.5.1)
* ```R ggpmisc``` (v0.6.0)
* ```PLINK``` (v1.9 and v2.0)
* ```flashpca``` (v2.0)
* ```picard``` (v3.3.0)
* ```bcftools``` (v1.8)
* ```bgzip``` (v1.11)
* ```perl``` (v5.32.1)
* ```gzip``` (v1.5)
* ```Python``` (v3.7.12)
* ```Python dask``` (v2022.02.0)
* ```Python pandas``` (v1.3.5)
* ```KING``` (v2.2.7)

An Anaconda environment containing all of these tools can be installed using the included ```genotype_qc_env.yml``` file.

### Data
Due to access restrictions for most individual-level genotype datasets, we will be using a multi-ancestry dataset consisting of 265 individuals of African and European ancestry that was compiled from participants in the Human Genome Diversity Project (HGDP). These data, as well as the code used to compile them, can be found in the ```sample_dataset``` folder. To simulate an array-based dataset, these sequencing data have been filtered to retain only SNPs directly genotyped on the Illumina MEGA array.

## QC workflow
An outline of the entire QC workflow is depicted in the following flowchart:

Starting from the raw genotype data, the entire QC workflow can be broken down into three distinct sub-steps:
* Pre-imputation QC
* Imputation
* Post-imputation QC

### Pre-imputation QC
1. Convert the VCF file to plink BED/BIM/FAM files:
```
plink2 --vcf hgdp_afr_eur_mega_snps.hg19.vcf.gz --make-bed --out hgdp_afr_eur_mega_snps.hg19
```

2. Print out the individual missingness of each individual in the dataset:
```
plink --bfile hgdp_afr_eur_mega_snps.hg19 --missing
```

Next, based on the degree of missingness, the a list was compiled of the duplicate individuals to be removed. Of each of the duplicate individuals, the one with a lower SNP missing rate was retained:
```
library("data.table")

path <- "plink.imiss"
file <- fread(path, header = TRUE, sep = " ", quote = "")
df <- as.data.frame(file)

df_dup <- df[duplicated(df$IID),]

new_df <- data.frame(matrix(NA, nrow = 0, ncol = 2))
names(new_df) <- c("FID", "IID")
unique_names <- unique(df_dup$IID)

for (name in unique_names) {
    temp_df <- df_dup[(df_dup$IID == name),]
    temp_df <- temp_df[temp_df$F_MISS == min(temp_df$F_MISS),]
    temp_df <- data.frame(temp_df$FID, temp_df$IID)
    names(temp_df) <- c("FID", "IID")
    new_df <- rbind(new_df, temp_df)
}

#Append the remaining individuals not duplicated
df_nodup <- df[!(df$IID %in% unique_names),]
df_nodup <- data.frame(df_nodup$FID, df_nodup$IID)
names(df_nodup) <- c("FID", "IID")

new_df <- rbind(new_df, df_nodup)

write.table(new_df, file = "hgdp_afr_eur_mega_snps_noduplicates.keep", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
```

Next, using this newly constructed list, keep all individuals in the file and filter out the remaining duplicates:
```
plink \
--bfile hgdp_afr_eur_mega_snps.hg19 \
--keep hgdp_afr_eur_mega_snps_noduplicates.keep \
--make-bed \
--out hgdp_afr_eur_mega_snps.hg19.nodup
```

In this case, there are no duplicate individuals. So we should see a printout similar to the following:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to hgdp_afr_eur_mega_snps.hg19.nodup.log.
Options in effect:
  --bfile hgdp_afr_eur_mega_snps.hg19
  --keep hgdp_afr_eur_mega_snps_noduplicates.keep
  --make-bed
  --out hgdp_afr_eur_mega_snps.hg19.nodup

64224 MB RAM detected; reserving 32112 MB for main workspace.
1270557 variants loaded from .bim file.
266 people (0 males, 0 females, 266 ambiguous) loaded from .fam.
Ambiguous sex IDs written to hgdp_afr_eur_mega_snps.hg19.nodup.nosex .
--keep: 266 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 266 founders and 0 nonfounders present.
Calculating allele freque