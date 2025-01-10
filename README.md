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

### QC overview
Starting from the raw genotype data, the entire QC workflow can be broken down into three distinct sub-steps:
* Pre-imputation QC
* Imputation
* Post-imputation QC

An outline of the entire QC workflow is depicted in the following flowchart:
![alt text](https://github.com/mjbetti/genotype_qc/blob/main/FigureS1.png?raw=true)

## Pre-imputation QC
**1.** Print out the individual missingness of each individual in the dataset:
```
plink --bfile gnomad_afr_eur_mega_snps.hg19 --missing
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

write.table(new_df, file = "gnomad_afr_eur_mega_snps_noduplicates.keep", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
```

Next, using this newly constructed list, keep all individuals in the file and filter out the remaining duplicates:
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19 \
--keep gnomad_afr_eur_mega_snps_noduplicates.keep \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup
```

In this case, there are no duplicate individuals. So we should see a printout similar to the following:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19
  --keep gnomad_afr_eur_mega_snps_noduplicates.keep
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup

64224 MB RAM detected; reserving 32112 MB for main workspace.
1305883 variants loaded from .bim file.
1435 people (0 males, 0 females, 1435 ambiguous) loaded from .fam.
Ambiguous sex IDs written to gnomad_afr_eur_mega_snps.hg19.nodup.nosex .
--keep: 1435 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1435 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.984844.
1305883 variants and 1435 people pass filters and QC.
Note: No phenotypes present.
```
\
**2.** Sometimes, samples are genotyped with reference individuals, such as those from HapMap or the CHM13 cell line, for assessing the quality of the assay. These reference individuals should be filtered out prior to continuing QC. In our case, a single CHM13 sample is included in our dataset.

Make file of CHM13 sample using known sample ID:
```
library("data.table")

fam_path <- "gnomad_afr_eur_mega_snps.hg19.nodup.fam"

#Open each of these files as a data frame
fam_file <- fread(fam_path, header = FALSE, sep = " ", quote = "")
fam_df <- as.data.frame(fam_file)
chm13 <- fam_df[startsWith(fam_df[,2], "CHM"),]
chm13 <- data.frame(chm13[,1], chm13[,2])

write.table(chm13, file = "gnomad_afr_eur_mega_snps.hg19.nodup.chm13.remove", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
```

The following plink command was next run to filter out the CHM13 individual from the dataset, retaining only HGDP individuals:
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup \
--remove gnomad_afr_eur_mega_snps.hg19.nodup.chm13.remove \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13
```
\
**3.** If available, sex (coded as 1 = male, 2 = female, 0 = unknown) and case-control status (coded as 1 = control, 2 = case, -9/0/non-numeric = missing) can be added to the fifth and sixth columns of the fam file, respectively. In this example, we have a covariate file called ```gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.sex.tsv``` with sex information compiled for all gnomAD+HGDP individuals. We will use this file to recode sex in the fam file (using R), while case-control status will stay coded as missing.
```
fam_path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.fam"
cov_path <- "gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.sex.tsv"

fam_file <- read.table(fam_path, header = FALSE, sep = " ", stringsAsFactors = FALSE)
fam_df <- as.data.frame(fam_file)

cov_file <- read.table(cov_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cov_df <- as.data.frame(cov_file)

merged_df <- merge(fam_df, cov_df, by.x = "V2", by.y = "sample")

merged_df[(merged_df$sex == "M"),"V5"] <- 1
merged_df[(merged_df$sex == "F"),"V5"] <- 2
merged_df[(merged_df$sex == "UNKNOWN"),"V5"] <- 0

merged_df <- merged_df[,c("V1", "V2", "V3", "V4", "V5", "V6")]

write.table(merged_df, file = fam_path, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
\
**4.** SNP pre-cleaning. Now that we are set with the files required for QC, we are able to begin the filtering phase. Using plink 1.9, the SNPs were pre-cleaned to retain only biallelic SNPs and filter by individual-level and SNP-level call rate. It is recommended to test multiple combinations of different call thresholds (i.e. 0.95 and 0.98) to find the optimal one that retains the highest number of individuals after filtering.

**SNP-level call rate of 0.98 (```--geno 0.02```) and individual-level call rate of 0.98 (```--mind 0.02```)**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13 \
--biallelic-only \
--geno 0.02 \
--mind 0.02 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02
```

Output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13
  --biallelic-only
  --geno 0.02
  --make-bed
  --mind 0.02
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02

64224 MB RAM detected; reserving 32112 MB for main workspace.
1305883 variants loaded from .bim file.
1434 people (753 males, 680 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02.nosex
.
85 people removed due to missing genotype data (--mind).
IDs written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02.irem
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1349 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3315 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.985585.
71830 variants removed due to missing genotype data (--geno).
1234053 variants and 1349 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02.bed
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02.bim
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.02.fam
... done.
```

**SNP-level call rate of 0.95 (```--geno 0.05```) and individual-level call rate of 0.95 (```--mind 0.05```)**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13 \
--biallelic-only \
--geno 0.05 \
--mind 0.05 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05
```

Output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13
  --biallelic-only
  --geno 0.05
  --make-bed
  --mind 0.05
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05

64224 MB RAM detected; reserving 32112 MB for main workspace.
1305883 variants loaded from .bim file.
1434 people (753 males, 680 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.nosex
.
0 people removed due to missing genotype data (--mind).
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1434 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3315 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.984818.
59605 variants removed due to missing genotype data (--geno).
1246278 variants and 1434 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.bed
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.bim
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.fam
... done.
```

**SNP-level call rate of 0.98 (```--geno 0.02```) and individual-level call rate of 0.95 (```--mind 0.05```)**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13 \
--biallelic-only \
--geno 0.02 \
--mind 0.05 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.05
```

Output:
```
64224 MB RAM detected; reserving 32112 MB for main workspace.
1305883 variants loaded from .bim file.
1434 people (753 males, 680 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.05.nosex
.
0 people removed due to missing genotype data (--mind).
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1434 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3315 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.05.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.984818.
80243 variants removed due to missing genotype data (--geno).
1225640 variants and 1434 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.05.bed
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.05.bim
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.02_mind.0.05.fam
... done.
```

**SNP-level call rate of 0.95 (```--geno 0.05```) and individual-level call rate of 0.98 (```--mind 0.02```)**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13 \
--biallelic-only \
--geno 0.05 \
--mind 0.02 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02
```

Output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13
  --biallelic-only
  --geno 0.05
  --make-bed
  --mind 0.02
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02

64224 MB RAM detected; reserving 32112 MB for main workspace.
1305883 variants loaded from .bim file.
1434 people (753 males, 680 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02.nosex
.
85 people removed due to missing genotype data (--mind).
IDs written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02.irem
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1349 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3315 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.985585.
55594 variants removed due to missing genotype data (--geno).
1250289 variants and 1349 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02.bed
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02.bim
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02.fam
... done.
```

We should find that a both SNP-level and individual-level call rate of 95% works best, as we retain both the largest number of individuals and SNPs. We should see that 59,605 variants and 0 individuals were removed.

**5.** Next, we will assign genetic ancestry-based population identifiers (African or European ancestry) to these individuals using principal components (PCs).

***a.*** A pre-compiled binary of FlashPCA2 can be downloaded from GitHub:
```
wget https://github.com/gabraham/flashpca/releases/download/v2.0/flashpca_x86-64.gz
gunzip flashpca_x86-64.gz
chmod 777 flashpca_x86-64
```
A copy is included in the ```scripts``` sub-directory of this repisotory.

***b.*** Before computing PCs, we will need to merge our filtered PLINK bed/bim/fam files with those containing reference individuals from the CEU, YRI, and CHB/JPT 1000 Genomes reference populations (representing European, African, and East Asian ancestries, respectively). A compiled set of PLINK files can be downloaded from Dropbox. A detailed explanation of how they were generated can be found at X.

***c.*** SNP pre-cleaning should be performed on the target set, which consists of retaining only biallelic, autosomal SNPs that meet specified SNP-level missingness, HWE, MAF, and LD thresholds:
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02 \
--biallelic-only \
--geno 0.02 \
--hwe 0.05 \
--maf 0.05 \
--indep-pairwise 80 8 0.15 \
--autosome \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15
 ```

We should observe the following output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.log.
Options in effect:
--autosome
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02
--biallelic-only
--geno 0.02
--hwe 0.05
--indep-pairwise 80 8 0.15
--maf 0.05
--make-bed
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 91816 MB successfully, after larger attempt(s) failed.
1214980 out of 1250289 variants loaded from .bim file.
1349 people (684 males, 664 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.nosex
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1349 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.999236.
15375 variants removed due to missing genotype data (--geno).
--hwe: 282665 variants removed due to Hardy-Weinberg exact test.
685110 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
231830 variants and 1349 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.bed
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.bim
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.fam
... done.
Pruned 8127 variants from chromosome 1, leaving 9150.
Pruned 9965 variants from chromosome 2, leaving 9478.
Pruned 9070 variants from chromosome 3, leaving 8134.
Pruned 7504 variants from chromosome 4, leaving 7585.
Pruned 6287 variants from chromosome 5, leaving 7248.
Pruned 10561 variants from chromosome 6, leaving 7139.
Pruned 6664 variants from chromosome 7, leaving 6635.
Pruned 5931 variants from chromosome 8, leaving 6178.
Pruned 4992 variants from chromosome 9, leaving 5502.
Pruned 5604 variants from chromosome 10, leaving 6154.
Pruned 5455 variants from chromosome 11, leaving 5701.
Pruned 5028 variants from chromosome 12, leaving 5669.
Pruned 4023 variants from chromosome 13, leaving 4415.
Pruned 3364 variants from chromosome 14, leaving 4030.
Pruned 3142 variants from chromosome 15, leaving 3838.
Pruned 3606 variants from chromosome 16, leaving 4289.
Pruned 2718 variants from chromosome 17, leaving 3831.
Pruned 3150 variants from chromosome 18, leaving 3842.
Pruned 2187 variants from chromosome 19, leaving 2854.
Pruned 2743 variants from chromosome 20, leaving 3339.
Pruned 1403 variants from chromosome 21, leaving 1730.
Pruned 1545 variants from chromosome 22, leaving 2020.
Pruning complete.  113069 of 231830 variants removed.
Marker lists written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.prune.in
and
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.prune.out
.
```

***d.*** We should reformat the variant IDs in our PLINK bim file to match the formatting in the reference dataset:
```
library("data.table")

#Declare the path of the bim file
bim_path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.bim"

#Open the bim file as a data frame
bim_file <- fread(bim_path, header = FALSE, sep = "\t", quote = "")
bim_df <- as.data.frame(bim_file)

#Generate new variant IDs for each SNP using the coordinate and allele information in the other columns
bim_df[,2] <- paste(paste0("chr", bim_df[,1]), bim_df[,4], bim_df[,6], bim_df[,5], "b37", sep = "_")

#Write the modified bim out
write.table(bim_df, file = bim_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

***e.*** We also need to prune the reference 1000 Genomes individuals using the same thresholds as the target individuals:
```
plink \
--bfile ceu_yri_phase3_ceu_yri_chb_jpt \
--biallelic-only \
--geno 0.02 \
--hwe 0.05 \
--maf 0.05 \
--indep-pairwise 80 8 0.15 \
--autosome \
--make-bed \
--out ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15 \
--allow-extra-chr
```

We should observe the following output:
```
PLINK v1.90b5.2 64-bit (9 Jan 2018)            www.cog-genomics.org/plink/1.9/
(C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.log.
Options in effect:
--allow-extra-chr
--autosome
--bfile ceu_yri_phase3_ceu_yri_chb_jpt
--biallelic-only
--geno 0.02
--hwe 0.05
--indep-pairwise 80 8 0.15
--maf 0.05
--make-bed
--out ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 91816 MB successfully, after larger attempt(s) failed.
80855722 out of 84358431 variants loaded from .bim file.
414 people (203 males, 211 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 414 founders and 0 nonfounders present.
Calculating allele frequencies... done.
0 variants removed due to missing genotype data (--geno).
--hwe: 4488195 variants removed due to Hardy-Weinberg exact test.
72038340 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
4329187 variants and 414 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.bed
+
ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.bim
+
ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.fam
... done.
Pruned 309315 variants from chromosome 1, leaving 21012.
Pruned 325243 variants from chromosome 2, leaving 20145.
Pruned 290870 variants from chromosome 3, leaving 17691.
Pruned 303989 variants from chromosome 4, leaving 17080.
Pruned 263311 variants from chromosome 5, leaving 15576.
Pruned 281146 variants from chromosome 6, leaving 16053.
Pruned 246538 variants from chromosome 7, leaving 14854.
Pruned 219387 variants from chromosome 8, leaving 12998.
Pruned 176877 variants from chromosome 9, leaving 12008.
Pruned 210289 variants from chromosome 10, leaving 13161.
Pruned 208694 variants from chromosome 11, leaving 12134.
Pruned 194353 variants from chromosome 12, leaving 13029.
Pruned 150637 variants from chromosome 13, leaving 9320.
Pruned 130650 variants from chromosome 14, leaving 8921.
Pruned 113705 variants from chromosome 15, leaving 8674.
Pruned 123389 variants from chromosome 16, leaving 9317.
Pruned 108079 variants from chromosome 17, leaving 9073.
Pruned 116516 variants from chromosome 18, leaving 8265.
Pruned 94067 variants from chromosome 19, leaving 7659.
Pruned 85942 variants from chromosome 20, leaving 6535.
Pruned 59627 variants from chromosome 21, leaving 4112.
Pruned 54142 variants from chromosome 22, leaving 4804.
Pruning complete.  4066766 of 4329187 variants removed.
Marker lists written to
ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.prune.in
and
ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.prune.out
.
```

***f.*** The reference dataset (1000 Genomes CEU, YRI, and CHB/JPT populations) can now be pruned and merged with the target data. The LD cutoff for pruning was an r2 threshold of 0.1. The input files for the HIGHLDvariable are included in this repository in the ```pcs``` sub-folder and are borrowed from the plinkQC R package:
```
#Declare the paths of the processed AALC cohort and 1000 Genomes data within the Aldrich Lab folder on ACCRE
BED_IN_DIR=./
BED_IN_NAME=gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15
BED_1KG_DIR=./1kg_ceu_yri_chb_jpt_refs
BED_1KG_NAME=ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15

#Declare the path of file containing genomic ranges of high-LD structure. This file is included with the plinkQC R package.
HIGHLD=high-LD-regions-hg19-GRCh37.txt

#Create the necessary output directories for outputs
ROOT_DIR=./

QCDIR=$ROOT_DIR\/qcdir

mkdir $QCDIR
mkdir $QCDIR\/plink_log

#Merge the cohort genotype data with the 1000 Genomes reference data using plink --merge
    #Filter reference and study data for non A-T ot G-C SNPs. This is because these SNPs are more difficult to align, and only a subset of SNPs is required for this type of analysis.
        #Input cohort
awk 'BEGIN {OFS="\t"}($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA"){print $2}' \
$BED_IN_DIR\/$BED_IN_NAME\.bim > \
$QCDIR\/$BED_IN_NAME\.ac_gt_snps

plink --bfile $BED_IN_DIR\/$BED_IN_NAME \
    --exclude $QCDIR/$BED_IN_NAME\.ac_gt_snps \
    --make-bed \
    --out $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
    --allow-extra-chr

mv $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps.log $QCDIR\/plink_log/$BED_IN_NAME\.no_ac_gt_snps.log

        #1000 Genomes
awk 'BEGIN {OFS="\t"}($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA"){print $2}' \
$BED_1KG_DIR\/$BED_1KG_NAME\.bim > \
$QCDIR\/$BED_1KG_NAME\.ac_gt_snps

plink --bfile $BED_1KG_DIR\/$BED_1KG_NAME \
    --exclude $QCDIR\/$BED_1KG_NAME\.ac_gt_snps \
    --make-bed \
    --allow-extra-chr \
    --out $QCDIR\/$BED_1KG_NAME\.no_ac_gt_snps

mv $QCDIR\/$BED_1KG_NAME\.no_ac_gt_snps.log $QCDIR\/plink_log/$BED_1KG_NAME\.no_ac_gt_snps.log

    #Prune study data. Conduct PCA on genetic variants that are pruned for variants in LD with r2 > 0.1 in a 50 kb window.
        #Input cohort
plink --bfile $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
    --exclude range $HIGHLD \
    --indep-pairwise 80 8 0.15 \
    --out $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
    --allow-extra-chr
mv $QCDIR\/$BED_IN_NAME\.prune.log $QCDIR\/plink_log/$BED_IN_NAME\.prune.log

plink --bfile $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps \
    --extract $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps.prune.in \
    --make-bed \
    --out $QCDIR\/$BED_IN_NAME\.pruned \
    --allow-extra-chr
mv $QCDIR\/$BED_IN_NAME\.pruned.log $QCDIR\/plink_log/$BED_IN_NAME\.pruned.log

    #Filter reference data for the same SNP set as in study. Use the list of pruned variants from the study sample to reduce the reference dataset to the size of the study samples
plink --bfile $BED_1KG_DIR\/$BED_1KG_NAME \
    --extract $QCDIR\/$BED_IN_NAME\.no_ac_gt_snps.prune.in \
    --make-bed \
    --allow-extra-chr \
    --out $QCDIR\/$BED_1KG_NAME\.pruned

mv $QCDIR\/$BED_1KG_NAME\.pruned.log $QCDIR\/plink_log/$BED_1KG_NAME\.pruned.log

    #Check and correct chromosome mismatch. Check that the variant IDs of the reference data have the same chromosome ID as the study data. Merging the files via plink will only work for variants with perfectly matching attributes. Because sex chromosomes and are often encoded differently and might make the matching more difficult, we will ignore sex chromosomes.
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)    {print a[$2],$2}' \
    $QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.pruned.bim | \
    sed -n '/^[XY]/!p' > $QCDIR\/$BED_1KG_NAME\.toUpdateChr

plink --bfile $QCDIR\/$BED_1KG_NAME\.pruned \
    --update-chr $QCDIR\/$BED_1KG_NAME\.toUpdateChr 1 2 \
    --make-bed \
    --out $QCDIR\/$BED_1KG_NAME\.updateChr \
    --allow-extra-chr
mv $QCDIR\/$BED_1KG_NAME\.updateChr.log $QCDIR\/plink_log/$BED_1KG_NAME\.updateChr.log

    #Position mismatch - find variants with mis-matching chromosome positions
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)    {print a[$2],$2}' \
    $QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.pruned.bim > \
    $QCDIR\/$BED_1KG_NAME\.toUpdatePos

    #Possible allele flips - check if non-matching allele codes are a simple case of allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.pruned.bim > \
    $QCDIR\/$BED_1KG_NAME\.toFlip

    #Update positions and flip alleles
plink --bfile $QCDIR\/$BED_1KG_NAME\.updateChr \
    --update-map $QCDIR\/$BED_1KG_NAME\.toUpdatePos 1 2 \
    --flip $QCDIR\/$BED_1KG_NAME\.toFlip \
    --make-bed \
    --out $QCDIR\/$BED_1KG_NAME\.flipped \
    --allow-extra-chr

mv $QCDIR\/$BED_1KG_NAME\.flipped.log $QCDIR\/plink_log/$BED_1KG_NAME\.flipped.log

    #Remove mismatches - any alleles that do not match after allele flipping are identified and removed from the reference dataset.
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $QCDIR\/$BED_IN_NAME\.pruned.bim $QCDIR\/$BED_1KG_NAME\.flipped.bim > \
    $QCDIR\/$BED_1KG_NAME\.mismatch

plink --bfile $QCDIR\/$BED_1KG_NAME\.flipped \
    --exclude $QCDIR\/$BED_1KG_NAME\.mismatch \
    --make-bed \
    --out $QCDIR\/$BED_1KG_NAME\.clean \
    --allow-extra-chr

mv $QCDIR\/$BED_1KG_NAME\.clean.log $QCDIR\/plink_log/$BED_1KG_NAME\.clean.log

    #Merge study genotypes and reference data - merge the AALC study data and 1000 Genomes reference dataset into a combined dataset
plink --bfile $QCDIR\/$BED_IN_NAME\.pruned  \
    --bmerge $QCDIR\/$BED_1KG_NAME\.clean.bed $QCDIR\/$BED_1KG_NAME\.clean.bim \
        $QCDIR\/$BED_1KG_NAME\.clean.fam  \
    --make-bed \
    --out $ROOT_DIR\/$BED_IN_NAME\.merge.$BED_1KG_NAME \
    --allow-extra-chr

mv $ROOT_DIR\/$BED_IN_NAME\.merge.$BED_1KG_NAME\.log $QCDIR/plink_log
```

***g.*** We also needed to LD prune the SNPS again (–indep-pairwise 80 8 0.15 (80-SNP sliding window, 8 SNPs each time, R2 > 0.15)):
```
mv gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.merge.ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.bed for_ld_pruning.bed
mv gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.merge.ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.fam for_ld_pruning.fam
mv gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.02_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.merge.ceu_yri_phase3_ceu_yri_chb_jpt_biallelic_geno0.02_mind0.05_hwe0.05_maf0.05_indep80_8_0.15.bim for_ld_pruning.bim

#Import plink 1.9, which is installed as a module on ACCRE
module load PLINK/1.9b_5.2

plink \
--bfile for_ld_pruning \
--indep-pairwise 80 8 0.15 \
--make-bed \
--out ld_pruned
```

We should observe the following output:
```
PLINK v1.90b5.2 64-bit (9 Jan 2018)            www.cog-genomics.org/plink/1.9/
(C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ld_pruned.log.
Options in effect:
--bfile for_ld_pruning
--indep-pairwise 80 8 0.15
--make-bed
--out ld_pruned

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 6893 MB successfully, after larger attempt(s) failed.
108191 variants loaded from .bim file.
1763 people (887 males, 875 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to ld_pruned.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 1763 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.900592.
108191 variants and 1763 people pass filters and QC.
Note: No phenotypes present.
--make-bed to ld_pruned.bed + ld_pruned.bim + ld_pruned.fam ... done.
Pruned 138 variants from chromosome 1, leaving 8284.
Pruned 125 variants from chromosome 2, leaving 8400.
Pruned 133 variants from chromosome 3, leaving 7303.
Pruned 117 variants from chromosome 4, leaving 6926.
Pruned 103 variants from chromosome 5, leaving 6474.
Pruned 105 variants from chromosome 6, leaving 6173.
Pruned 97 variants from chromosome 7, leaving 5943.
Pruned 82 variants from chromosome 8, leaving 5303.
Pruned 81 variants from chromosome 9, leaving 4971.
Pruned 92 variants from chromosome 10, leaving 5563.
Pruned 78 variants from chromosome 11, leaving 5119.
Pruned 90 variants from chromosome 12, leaving 5082.
Pruned 99 variants from chromosome 13, leaving 3984.
Pruned 65 variants from chromosome 14, leaving 3651.
Pruned 52 variants from chromosome 15, leaving 3459.
Pruned 47 variants from chromosome 16, leaving 3833.
Pruned 59 variants from chromosome 17, leaving 3473.
Pruned 63 variants from chromosome 18, leaving 3486.
Pruned 36 variants from chromosome 19, leaving 2610.
Pruned 47 variants from chromosome 20, leaving 2988.
Pruned 37 variants from chromosome 21, leaving 1561.
Pruned 21 variants from chromosome 22, leaving 1838.
Pruning complete.  1767 of 108191 variants removed.
Marker lists written to ld_pruned.prune.in and ld_pruned.prune.out .
```

***h.*** Next, we need to generate a .pop file indicating the target samples and population identifiers of the reference samples:
```
#Declare the file path for the merged fam file with select founder populations and cohort individuals labeled
fam_path <- "ld_pruned.fam"

#Open the fam file as a data frame and slice out the first column containing founder populations
fam_file <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
fam_df <- as.data.frame(fam_file, stringsAsFactors = FALSE)

pops <- fam_df[,1]

#For any individual that does not have a population listed as ASW or CEU, replace the value in pops with "-"
is_not_founder <- !(pops == "CEU" | pops == "YRI" | pops == "CHB" | pops == "JPT")

pops[is_not_founder] <- "-"

#Write out the pops data frame to a .pop file
out_name <- "ld_pruned.pop"
write.table(pops, file = out_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

***i.*** After generating the pop file, the fam file should finally be modified one last time to replace the family IDs (FIDs) with the individual IDs (IIDs), so that the first and second column have identical values. This can be performed using the following R script:
```
#Declare the path of the fam file (for plink set of SCCS cohort merged with CEU and YRI 1000 Genomes individuals)
fam_path <- "ld_pruned.fam"

#Open the fam file as a data frame
fam_file <- read.table(fam_path, header = FALSE)
fam_df <- as.data.frame(fam_file)

fam_df[,1] <- fam_df[,2]

write.table(fam_df, file = fam_path, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

h. We can now run FlashPCA2 using this LD-pruned dataset:
```
flashpca_x86-64 --bfile ld_pruned
```

The PCAs should be output to a file called ```pcs.txt```.

***j.*** Once PCs are computed, we can use the ```rescale_and_filter_pca.R``` script in the ```pca``` directory to assign population identifiers to the target individuals based on specified PC cutoff thresholds.
    
First, we need to make a file for the ```--persongroupfilepath``` argument using the following R script:
```
pca_path <- "pcs.txt"
pops_path <- "ld_pruned.pop"

pca_file <- read.table(pca_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pca_df <- as.data.frame(pca_file)

pops_file <- read.table(pops_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
pops_df <- as.data.frame(pops_file)

new_df <- data.frame(pca_df$IID, pops_df[,1])
new_df[new_df[,2] == "-",2] <- "gnomad"

new_df <- new_df[new_df[,2] == "gnomad",]

write.table(new_df, file = "gnomad_person_group_nopop.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

Next, we need to make a file for the ```--genomes1ksamplepoppath``` argument:
```
pca_path <- "pcs.txt"
pops_path <- "ld_pruned.pop"

pca_file <- read.table(pca_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pca_df <- as.data.frame(pca_file)

pops_file <- read.table(pops_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
pops_df <- as.data.frame(pops_file)

new_df <- data.frame(pca_df$IID, pops_df[,1])
new_df[new_df[,2] == "-",2] <- "gnomad"
new_df[new_df[,2] == "-",2] <- "gnomad"

new_df <- new_df[!(new_df[,2] == "gnomad"),]
new_df[new_df[,2] == "CHB",2] <- "CHBJPT"
new_df[new_df[,2] == "JPT",2] <- "CHBJPT"

new_df <- data.frame(new_df[,1], "GROUP", new_df[,2])
new_df[,2][(new_df[,3] == "CEU")] <- "EUR"
new_df[,2][(new_df[,3] == "YRI")] <- "AFR"
new_df[,2][(new_df[,3] == "CHBJPT")] <- "EAS"

names(new_df) <- c("ID", "GROUP", "POP")

write.table(new_df, file = "genomes1k_sample_pop.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
```

Finally, we need to make a file for the ```--genomes1krefgroups``` file. Use ```nano``` to open a new file:
```
nano genomes1krefgroups.txt
```
    
...and add the following test:
```
EUR-CEU,AFR-YRI,EAS-CHB
```

Finally, using these files we just made, we can re-scale the PCA results relative to the reference individuals:
```
Rscript ~/rescale_and_filter_pca.R \
--workingdir=./ \
--pcaresultspath=pcs.txt \
--persongroupfilepath=gnomad_person_group_nopop.txt \
--genomes1ksamplepoppath=genomes1k_sample_pop.txt \
--genomes1krefgroups='EUR-CEU,AFR-YRI,EAS-CHBJPT' \
--plottinggroups='gnomad'
```
    
Running this script should output the following PNG image showing the target gnomAD individuals plotted relative to the CEU, YRI, and CHB/JPT reference populations:
![alt text](https://github.com/mjbetti/genotype_qc/blob/main/pcs.txt.nonselected_scaled_PC1_vs_PC2.gnomad.png?raw=true)

I could now use this script again to assign individuals to European and African ancestry clusters. Specifically, we will classify individuals within 50% of the European ancestry reference cluster on both the CEU-YRI and CEU-CHB/JPT axes as having European ancestry. Individuals within 50% of the African ancestry reference cluster on both the YRI-CEU and YRI-CHB/JPT axes will be classified as having African ancestry.

**Identifying individuals of European ancestry**

First, we should make a ```pcs.txt.nonselected_scaled_PC1_vs_PC2_assigned_group_boundaries.EUR.50.50.txt``` file to classify European ancestry individuals and set the datapoints to the following values:
```
nano pcs.txt.nonselected_scaled_PC1_vs_PC2_assigned_group_boundaries.EUR.50.50.txt
```

This file should contain the following header and line (tab-separated):
```
assignment.group.name   first.PC        second.PC       EUR-CEU_vs_AFR-YRI.EUR-CEU      EUR-CEU_vs_AFR-YRI.AFR-YRI      EUR-CEU_vs_EAS-CHBJPT.EUR-CEU   EUR-CEU_vs_EAS-CHBJPT.EAS-CHBJPT        AFR-YRI_vs_EAS-CHBJPT.AFR-YRI   AFR-YRI_vs_EAS-CHBJPT.EAS-CHBJPT        assigned.eth    assigned.ref.pop.for.refilter
default PC1 PC2 50  NA  50  NA  NA  NA  EUR.assign.50.50    EUR-CEU,AFR-YRI,EAS-CHBJPT
```

Specific details regarding what each of these values means can be found in the comments of the ```rescale_and_filter_pca.R``` script.

We can now run the script to first generate a re-filtered subset of European ancestry individuals:
```
Rscript ~/rescale_and_filter_pca.R \
--workingdir=./ \
--pcaresultspath=pcs.txt \
--persongroupfilepath=gnomad_person_group_nopop.txt \
--genomes1ksamplepoppath=genomes1k_sample_pop.txt \
--genomes1krefgroups='EUR-CEU,AFR-YRI,EAS-CHBJPT' \
--plottinggroups='gnomad' \
--assignedgroups=EUR.assign.50.50 \
--groupassignmentfilepath=pcs.txt.nonselected_scaled_PC1_vs_PC2_assigned_group_boundaries.EUR.50.50.txt \
--rescaledpcaresultspath=pcs.txt.nonselected.scaled.pcs.txt
```

Find the number of remaining European ancestry gnomAD samples after filtering:
```
wc -l pcs.txt.rescaled_pca.resubsetted.EUR.assign.50.50.refilter.no.refs.keep.txt
606 pcs.txt.rescaled_pca.resubsetted.EUR.assign.50.50.refilter.no.refs.keep.txt
```

So 742 individuals would be classified as having primarily European ancestry based on PCs.

**Identifying individuals of African ancestry**

Next, we should make a ```pcs.txt.nonselected_scaled_PC1_vs_PC2_assigned_group_boundaries.AFR.50.50.txt``` file to classify African ancestry individuals and set the datapoints to the following values:
```
nano pcs.txt.nonselected_scaled_PC1_vs_PC2_assigned_group_boundaries.AFR.50.50.txt
```

This file should contain the following header and line (tab-separated):
```
assignment.group.name   first.PC        second.PC       EUR-CEU_vs_AFR-YRI.EUR-CEU      EUR-CEU_vs_AFR-YRI.AFR-YRI      EUR-CEU_vs_EAS-CHBJPT.EUR-CEU   EUR-CEU_vs_EAS-CHBJPT.EAS-CHBJPT        AFR-YRI_vs_EAS-CHBJPT.AFR-YRI   AFR-YRI_vs_EAS-CHBJPT.EAS-CHBJPT        assigned.eth    assigned.ref.pop.for.refilter
default PC1     PC2     NA	50	NA	NA	50	NA	AFR.assign.50.50        EUR-CEU,AFR-YRI,EAS-CHBJPT
```

Specific details regarding what each of these values means can be found in the comments of the ```rescale_and_filter_pca.R``` script.

We can now run the script to first generate a re-filtered subset of European ancestry individuals:
```
Rscript ~/rescale_and_filter_pca.R \
--workingdir=./ \
--pcaresultspath=pcs.txt \
--persongroupfilepath=gnomad_person_group_nopop.txt \
--genomes1ksamplepoppath=genomes1k_sample_pop.txt \
--genomes1krefgroups='EUR-CEU,AFR-YRI,EAS-CHBJPT' \
--plottinggroups='gnomad' \
--assignedgroups=AFR.assign.50.50 \
--groupassignmentfilepath=pcs.txt.nonselected_scaled_PC1_vs_PC2_assigned_group_boundaries.AFR.50.50.txt \
--rescaledpcaresultspath=pcs.txt.nonselected.scaled.pcs.txt
```

Find the number of remaining European ancestry gnomAD samples after filtering:
```
wc -l pcs.txt.rescaled_pca.resubsetted.AFR.assign.50.50.refilter.no.refs.keep.txt
742 pcs.txt.rescaled_pca.resubsetted.AFR.assign.50.50.refilter.no.refs.keep.txt
```

So 742 individuals would be classified as having primarily African ancestry based on PCs.

**Identifying whether any individuals were dropped**

We can use the following R script to identify whether any individuals were dropped (i.e. whether any of these individuals did not meet the criteria to be classified as either European or African ancestry, based on PCs):
```
all_pc_path <- "pcs.txt"
eur_pc_path <- "pcs.txt.rescaled_pca.resubsetted.EUR.assign.50.50.refilter.no.refs.keep.txt"
afr_pc_path <- "pcs.txt.rescaled_pca.resubsetted.AFR.assign.50.50.refilter.no.refs.keep.txt"

all_pc_file <- read.table(all_pc_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_pc_df <- as.data.frame(all_pc_file)

eur_pc_file <- read.table(eur_pc_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
eur_pc_df <- as.data.frame(eur_pc_file)

afr_pc_file <- read.table(afr_pc_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
afr_pc_df <- as.data.frame(afr_pc_file)

classified_df <- rbind(eur_pc_df, afr_pc_df)

length(all_pc_df[!(all_pc_df$IID %in% classified_df$V2),"IID"])
```

We should get the following output:
```
[1] 415
```

So of the original 1,763 gnomAD individuals in our dataset, 415 cannot definitively be classified as European or African ancestry based on PCs.

***k.*** Finally, we can split up they genotype data from our gnomAD individuals into separate sets of European and African ancestry individuals. The genotype data we will use will be the optimal set after variant and individual-level thresholding (generated prior to LD pruning):
```
#European ancestry
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05 \
--keep pcs.txt.rescaled_pca.resubsetted.EUR.assign.50.50.refilter.no.refs.keep.txt \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur

#African ancestry
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05 \
--keep pcs.txt.rescaled_pca.resubsetted.AFR.assign.50.50.refilter.no.refs.keep.txt \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr
```

We should observe the following output:
```
PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05
  --keep pcs.txt.rescaled_pca.resubsetted.EUR.assign.50.50.refilter.no.refs.keep.txt
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur

Start time: Mon Jan  6 14:31:33 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Using up to 72 threads (change this with --threads).
1434 samples (680 females, 753 males, 1 ambiguous; 1434 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.fam.
1246278 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.bim.
Note: No phenotype data present.
--keep: 606 samples remaining.
606 samples (307 females, 298 males, 1 ambiguous; 606 founders) remaining after
main filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur.bed
... done.
End time: Mon Jan  6 14:31:35 2025

PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05
  --keep pcs.txt.rescaled_pca.resubsetted.AFR.assign.50.50.refilter.no.refs.keep.txt
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr

Start time: Mon Jan  6 14:31:36 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Using up to 72 threads (change this with --threads).
1434 samples (680 females, 753 males, 1 ambiguous; 1434 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.fam.
1246278 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.bim.
Note: No phenotype data present.
--keep: 742 samples remaining.
742 samples (357 females, 385 males; 742 founders) remaining after main
filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr.bed
... done.
End time: Mon Jan  6 14:31:37 2025

```

We should now have the following two sets of ancestry-stratified plink bed/bim/fam files:
```
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur

gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr
```
\
**6.** Each of the ancestry-stratified datasets should next undergo additional variant-level filtering based on deviation from Hardy-Weinberg equilibrium (HWE). We will use a HWE p-value threshold of 1 x 10<sup>-8</sup> for European ancestry individuals and 1 x 10<sup>-6</sup> for African ancestry individuals:

**European ancestry**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur \
--hwe 1e-8 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8
```

We should observe the following output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur
  --hwe 1e-8
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 386918 MB successfully, after larger attempt(s) failed.
1246278 variants loaded from .bim file.
606 people (298 males, 307 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.nosex
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 606 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3123 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.999166.
Warning: --hwe observation counts vary by more than 10%.  Consider using
--geno, and/or applying different p-value thresholds to distinct subsets of
your data.
--hwe: 1588 variants removed due to Hardy-Weinberg exact test.
1244690 variants and 606 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.bed
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.bim
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.fam
... done.
```

**African ancestry**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr \
--hwe 1e-6 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6
```

We should observe the following output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr
  --hwe 1e-6
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 386918 MB successfully, after larger attempt(s) failed.
1246278 variants loaded from .bim file.
742 people (385 males, 357 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 742 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 7 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.999324.
Warning: --hwe observation counts vary by more than 10%, due to the X
chromosome.  You may want to use a less stringent --hwe p-value threshold for X
chromosome variants.
--hwe: 2586 variants removed due to Hardy-Weinberg exact test.
1243692 variants and 742 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6.bed
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6.bim
+
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6.fam
... done.
```

So variant filtering using HWE removed an additional 1,588 SNPs from the European ancestry dataset and 2,586 from the African ancestry dataset.
\
**7.** Next, the European and African ancestry individuals should each filtered based on minor allele frequency (MAF), with only SNPs having a MAF > 0.01 in each respective ancestry being retained:
\
**European ancestry**
```
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8 \
--maf 0.01 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01
```

We should observe the following output:
```
PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8
  --maf 0.01
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01

Start time: Mon Jan  6 16:17:47 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Allocated 386918 MiB successfully, after larger attempt(s) failed.
Using up to 72 threads (change this with --threads).
606 samples (307 females, 298 males, 1 ambiguous; 606 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.fam.
1244690 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8.bim.
Note: No phenotype data present.
Calculating allele frequencies... done.
662365 variants removed due to allele frequency threshold(s)
(--maf/--max-maf/--mac/--max-mac).
582325 variants remaining after main filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.bed
... done.
```
\
**African ancestry**
```
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6 \
--maf 0.01 \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01
```

We should observe the following output:
```
PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6
  --maf 0.01
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01

Start time: Mon Jan  6 16:19:01 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Allocated 386918 MiB successfully, after larger attempt(s) failed.
Using up to 72 threads (change this with --threads).
742 samples (357 females, 385 males; 742 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6.fam.
1243692 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6.bim.
Note: No phenotype data present.
Calculating allele frequencies... done.
501586 variants removed due to allele frequency threshold(s)
(--maf/--max-maf/--mac/--max-mac).
742106 variants remaining after main filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.bed
... done.
```

So variant filtering using HWE removed an additional 662,365 SNPs from the European ancestry dataset and 742,106 from the African ancestry dataset.

\
**8.** Next, we should compute the F<sub>het</sub> statistic for each ancestry-stratified dataset and remove individuals with an F coefficient more than 3 SDs from the group mean:

**European ancestry**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01 \
--het \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01
```

We should observe the following output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01
  --het
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 386918 MB successfully, after larger attempt(s) failed.
582325 variants loaded from .bim file.
606 people (298 males, 307 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.nosex
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 606 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3062 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.999143.
582325 variants and 606 people pass filters and QC.
Note: No phenotypes present.
--het: 566991 variants scanned, report written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.het
.
```

We can use the following R script to determine the individuals to prune out:
```
options(scipen = 100)
library("data.table")

path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.het"
file <- fread(path, header = TRUE, sep = " ", quote = "")
df <- as.data.frame(file)

df_mean <- mean(df$F)
df_sd <- sd(df$F)

upper_bound <- df_mean + (3 * df_sd)
lower_bound <- df_mean - (3 * df_sd)

df_pruned <- df[df$F > upper_bound | df$F < lower_bound,]

df_pruned <- df_pruned[,1:2]

write.table(df_pruned, file = "eur_fhet.remove", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

We should identify 10 European ancestry individuals to prune out, which have been output to the ```eur_fhet.remove``` file. We can then prune these individuals using plink:
```
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01 \
--remove eur_fhet.remove \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet
```

We should observe the following output:
```
PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet
  --remove eur_fhet.remove

Start time: Mon Jan  6 16:31:26 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Allocated 386918 MiB successfully, after larger attempt(s) failed.
Using up to 72 threads (change this with --threads).
606 samples (307 females, 298 males, 1 ambiguous; 606 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.fam.
582325 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01.bim.
Note: No phenotype data present.
--remove: 596 samples remaining.
596 samples (302 females, 293 males, 1 ambiguous; 596 founders) remaining after
main filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.bed
... done.
```
\
**African ancestry**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01 \
--het \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01
```

We should observe the following output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01
  --het
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 386918 MB successfully, after larger attempt(s) failed.
742106 variants loaded from .bim file.
742 people (385 males, 357 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 742 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 7 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.hh
); many commands treat these as missing.
Total genotyping rate is 0.999246.
742106 variants and 742 people pass filters and QC.
Note: No phenotypes present.
--het: 720010 variants scanned, report written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.het
.
```

We can use the following R script to determine the individuals to prune out:
```
options(scipen = 100)
library("data.table")

path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.het"
file <- fread(path, header = TRUE, sep = " ", quote = "")
df <- as.data.frame(file)

df_mean <- mean(df$F)
df_sd <- sd(df$F)

upper_bound <- df_mean + (3 * df_sd)
lower_bound <- df_mean - (3 * df_sd)

df_pruned <- df[df$F > upper_bound | df$F < lower_bound,]

df_pruned <- df_pruned[,1:2]

write.table(df_pruned, file = "afr_fhet.remove", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

We should identify 11 European ancestry individuals to prune out, which have been output to the ```afr_fhet.remove``` file. We can then prune these individuals using plink:
```
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01 \
--remove afr_fhet.remove \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet
```

We should observe the following output:
```
PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet
  --remove afr_fhet.remove

Start time: Mon Jan  6 16:44:14 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Allocated 386918 MiB successfully, after larger attempt(s) failed.
Using up to 72 threads (change this with --threads).
742 samples (357 females, 385 males; 742 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.fam.
742106 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01.bim.
Note: No phenotype data present.
--remove: 731 samples remaining.
731 samples (354 females, 377 males; 731 founders) remaining after main
filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.bed
... done.
```
\
**9.** Next, we will perform a sex check to test whether reported sex coincides with genetically inferred sex in each ancestry-stratified dataset. For the purpose of running the sex check, we will implement an additional stringent MAF threshold of 0.2:

**European ancestry**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--maf 0.2 \
--check-sex \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet
```

We should observe the following output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet
  --check-sex
  --maf 0.2
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 386918 MB successfully, after larger attempt(s) failed.
582325 variants loaded from .bim file.
596 people (293 males, 302 females, 1 ambiguous) loaded from .fam.
Ambiguous sex ID written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.nosex
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 596 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3062 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.hh
); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.999148.
409165 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
173160 variants and 596 people pass filters and QC.
Note: No phenotypes present.
--check-sex: 4203 Xchr and 0 Ychr variant(s) scanned, 7 problems detected.
Report written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.sexcheck
.
```

There are 7 problems detected, meaning that 7 individuals have a genetically inferred sex that differs from the reported sex. 

There are two possible solutions for resolving these discrepancies:

1. We can remove these discordant individuals from subsequent analysis
2. For discordant individuals, we can recode sex to that which was genetically inferred.

First, we will walk through an example of removing discordant individuals altogether. We will compile a list of these individuals to remove using the following R script:
```
library("data.table")

path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.sexcheck"
file <- fread(path, header = TRUE, sep = " ", quote = "")
df <- as.data.frame(file)

remove_df <- df[!(df$STATUS == "OK"),1:2]

write.table(remove_df, file = "eur_sexcheck.remove", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

Using this list, we can prune out the 7 discordant individuals using plink:
```
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--remove eur_sexcheck.remove \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet_sexcheck
```

We should observe the following output:
```
PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet_sexcheck.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet_sexcheck
  --remove eur_sexcheck.remove

Start time: Mon Jan  6 17:07:09 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Allocated 386918 MiB successfully, after larger attempt(s) failed.
Using up to 72 threads (change this with --threads).
596 samples (302 females, 293 males, 1 ambiguous; 596 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam.
582325 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.bim.
Note: No phenotype data present.
--remove: 589 samples remaining.
589 samples (297 females, 292 males; 589 founders) remaining after main
filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet_sexcheck.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet_sexcheck.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet_sexcheck.bed
... done.
End time: Mon Jan  6 17:07:09 2025
```

Alternatively, we can recode the sex of discordant individuals in the fam file to the genetically inferred sex using the following R script:
```
library("data.table")

sexcheck_path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.sexcheck"
sexcheck_file <- fread(sexcheck_path, header = TRUE, sep = " ", quote = "")
sexcheck_df <- as.data.frame(sexcheck_file)

fam_path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = "")
fam_df <- as.data.frame(fam_file)

fam_df$V5 <- sexcheck_df$SNPSEX

write.table(fam_df, file = fam_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```
\
**African ancestry**
```
plink \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--maf 0.2 \
--check-sex \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet
```

We should observe the following output:
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet
  --check-sex
  --maf 0.2
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet

1031782 MB RAM detected; reserving 515891 MB for main workspace.
Allocated 386918 MB successfully, after larger attempt(s) failed.
742106 variants loaded from .bim file.
731 people (377 males, 354 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 731 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 7 het. haploid genotypes present (see
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.hh
); many commands treat these as missing.
Total genotyping rate is 0.999265.
550600 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
191506 variants and 731 people pass filters and QC.
Note: No phenotypes present.
--check-sex: 5754 Xchr and 0 Ychr variant(s) scanned, 5 problems detected.
Report written to
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.sexcheck
.
```

There are 5 problems detected, meaning that 5 individuals have a genetically inferred sex that differs from the reported sex. 

There are two possible solutions for resolving these discrepancies:

1. We can remove these discordant individuals from subsequent analysis
2. For discordant individuals, we can recode sex to that which was genetically inferred.

First, we will walk through an example of removing discordant individuals altogether. We will compile a list of these individuals to remove using the following R script:
```
library("data.table")

path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.sexcheck"
file <- fread(path, header = TRUE, sep = " ", quote = "")
df <- as.data.frame(file)

remove_df <- df[!(df$STATUS == "OK"),1:2]

write.table(remove_df, file = "afr_sexcheck.remove", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

Using this list, we can prune out the 5 discordant individuals using plink:
```
plink2 \
--bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--remove afr_sexcheck.remove \
--make-bed \
--out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet_sexcheck
```

We should observe the following output:
```
PLINK v2.00a2.3LM 64-bit Intel (24 Jan 2020)   www.cog-genomics.org/plink/2.0/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet_sexcheck.log.
Options in effect:
  --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet
  --make-bed
  --out gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet_sexcheck
  --remove afr_sexcheck.remove

Start time: Mon Jan  6 17:13:06 2025
1031782 MiB RAM detected; reserving 515891 MiB for main workspace.
Allocated 386918 MiB successfully, after larger attempt(s) failed.
Using up to 72 threads (change this with --threads).
731 samples (354 females, 377 males; 731 founders) loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam.
742106 variants loaded from
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.bim.
Note: No phenotype data present.
--remove: 726 samples remaining.
726 samples (349 females, 377 males; 726 founders) remaining after main
filters.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet_sexcheck.fam
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet_sexcheck.bim
... done.
Writing
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet_sexcheck.bed
... done.
```

Alternatively, we can recode the sex of discordant individuals in the fam file to the genetically inferred sex using the following R script:
```
library("data.table")

sexcheck_path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.sexcheck"
sexcheck_file <- fread(sexcheck_path, header = TRUE, sep = " ", quote = "")
sexcheck_df <- as.data.frame(sexcheck_file)

fam_path <- "gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = "")
fam_df <- as.data.frame(fam_file)

fam_df$V5 <- sexcheck_df$SNPSEX

write.table(fam_df, file = fam_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```
\
Ultimately, the approach you decide to use depends on the specifics of your analysis and the level of stringency you prefer to use.

In this example, we will move forward with the recoded datasets, saved to the following file paths:
```
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet

gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet
```
\
**10.** Finally, before genotype imputation, you also may want to test for batch effects. This is particularly relevant if you have DNA derived from various sample types (i.e. blood vs. saliva vs. mouthwash). We will not test for batch effects here, but sample code for such a scenario can be found in the ```batch_effects``` sub-directory.

## Imputation
The pre-imputed, ancestry-stratified datasets on which we performed QC should be stored at the following paths:

**European ancestry**
```
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet
```

**African ancestry**
```
gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet
```
\
**1.** These variants should be converted to VCF format before being lifted over using the liftoverVCF function from Picard/GATK:
```
mkdir imputation

plink2 --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--recode vcf \
--out imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--snps-only just-acgt \
--autosome

plink2 --bfile gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--recode vcf \
--out imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--snps-only just-acgt \
--autosome
```

The resulting VCF files are saved to the following path:
```
imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.vcf

imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.vcf
```
\
**2.** Many array manifests (i.e. Illumina MEGA) use SNP coordinates based on genome build hg19, while the TOPMed Imputation Server expects coordinates based on hg38. If this is the case, we will need to lift the VCF from hg19 to hg38 prior to imputation. 

Picard has a ```LiftoverVcf``` function that can be used to lift VCF coordinates from one genome build to another.

The following inputs are required to run the LiftoverVcf function:
- liftOver chain file (hg19 to hg38)
- input VCF
- output (desired location of the lifted over VCF)
- reference sequence (FASTA and corresponding .dict file for target genome build hg38)
- reject (file path where rejected variant [i.e. those that cannot be lifted over] will be written)

An hg19 FASTA file can be downloaded from the UCSC Genome Browser (and is also available from the linked Dropbox here):
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip hg38.fa
```

Prior to running the liftover, another Picard function, ```CreateSequenceDictionary```, was used to create the .dict file for the hg38 FASTA:
```
picard CreateSequenceDictionary \
R=hg38.fa \
O=hg38.dict
```

The ```liftOver``` script will expect genomic coordinates to have the prefix ```chr``` in front of chromosome names. The following scripts were written to add the "chr" prefix to the VCF files before lifting over:

**European ancestry VCF**
\
First output the VCF header to a new file:
```
grep '^##' imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.vcf > imputation/eur_correct_header.txt
echo '##INFO=<ID=ReverseComplementedAlleles,Number=1,Type=Integer,Description="Flag for Reverse Complemented Alleles">' >> imputation/eur_correct_header.txt
grep '^#CHROM' imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.vcf >> imputation/eur_correct_header.txt
```

Split up the VCF into individual files for each chromosome:
```
#Split up by chromosome
for i in {1..22}; do
sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=60G \
    --time=1-00:00:00 \
    --job-name="${i}_split_eur" \
    --wrap="plink2 --vcf imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.vcf --chr $i --recode vcf --out imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet"
done
```
Next, we can add the ```chr``` prefix to chromosome identifiers in each VCF:
```
for i in {1..22}; do
    sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=60G \
    --time=1-00:00:00 \
    --job-name="${i}_liftover" \
    --wrap=$'awk \'{if($0 !~ /^#/) print "chr"$0; else print $0}\' imputation/chr'$i'_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.vcf > imputation/chr'$i'_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.vcf'
done
```

Finally, we can lift over each VCF from hg19 to hg38:
```
for i in {1..22}; do
sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=60G \
--time=1-00:00:00 \
--job-name=$i\_liftover \
--wrap="picard -Xmx60g LiftoverVcf \
I=imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.vcf \
O=imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.vcf \
CHAIN=hg19ToHg38.over.chain.gz \
REJECT=imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg19.reject.vcf \
R=hg38.fa"
done
```

After performing the liftover, we should filter the hg38-based VCF files to ensure that there are no alternative chromosomes/contigs introduced:
```
for i in {1..22}; do
sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=60G \
--time=1-00:00:00 \
--job-name=$i\_liftover \
--wrap="plink2 \
--vcf imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.vcf \
--chr $i \
--recode vcf \
--allow-extra-chr \
--out imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt"
done
```

Using the new header we made above, we will replace the header of each new hg38-based VCF:
```
for i in {1..22}; do
    sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=60G \
    --time=1-00:00:00 \
    --job-name=$i \
    --wrap="bcftools reheader -h imputation/eur_correct_header.txt -o imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.vcf"
done
```

Next, we will concatenate all of the individual VCFs into a single file:
```
IN_DIR=imputation

bcftools concat -Ov -o $IN_DIR\/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr1_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr2_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr3_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr4_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr5_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr6_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr7_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr8_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr9_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr10_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr11_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr12_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr13_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr14_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr15_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr16_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr17_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr18_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr19_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr20_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr21_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr22_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf
```

The lifted over VCF will be saved to the following path:
```
imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf
```

**African ancestry VCF**

First output the VCF header to a new file:
```
grep '^##' imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.vcf > imputation/afr_correct_header.txt
echo '##INFO=<ID=ReverseComplementedAlleles,Number=1,Type=Integer,Description="Flag for Reverse Complemented Alleles">' >> imputation/afr_correct_header.txt
grep '^#CHROM' imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.vcf >> imputation/afr_correct_header.txt
```

Split up the VCF into individual files for each chromosome:
```
#Split up by chromosome
for i in {1..22}; do
sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=60G \
    --time=1-00:00:00 \
    --job-name="${i}_split_afr" \
    --wrap="plink2 --vcf imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.vcf --chr $i --recode vcf --out imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet"
done
```
Next, we can add the ```chr``` prefix to chromosome identifiers in each VCF:
```
for i in {1..22}; do
    sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=60G \
    --time=1-00:00:00 \
    --job-name="${i}_liftover" \
    --wrap=$'awk \'{if($0 !~ /^#/) print "chr"$0; else print $0}\' imputation/chr'$i'_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.vcf > imputation/chr'$i'_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.vcf'
done
```

Finally, we can lift over each VCF from hg19 to hg38:
```
for i in {1..22}; do
sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=60G \
--time=1-00:00:00 \
--job-name=$i\_liftover \
--wrap="picard -Xmx60g LiftoverVcf \
I=imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.vcf \
O=imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.vcf \
CHAIN=hg19ToHg38.over.chain.gz \
REJECT=imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg19.reject.vcf \
R=hg38.fa"
done
```

After performing the liftover, we should filter the hg38-based VCF files to ensure that there are no alternative chromosomes/contigs introduced:
```
for i in {1..22}; do
sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=60G \
--time=1-00:00:00 \
--job-name=$i\_liftover \
--wrap="plink2 \
--vcf imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.vcf \
--chr $i \
--recode vcf \
--allow-extra-chr \
--out imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt"
done
```

Using the new header we made above, we will replace the header of each new hg38-based VCF:
```
for i in {1..22}; do
    sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=60G \
    --time=1-00:00:00 \
    --job-name=$i \
    --wrap="bcftools reheader -h imputation/afr_correct_header.txt -o imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf imputation/chr$i\_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.vcf"
done
```

Next, we will concatenate all of the individual VCFs into a single file:
```
IN_DIR=imputation

bcftools concat -Ov -o $IN_DIR\/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr1_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr2_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr3_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr4_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr5_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr6_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr7_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr8_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr9_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr10_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr11_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr12_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr13_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr14_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr15_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr16_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr17_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr18_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr19_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr20_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr21_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf \
$IN_DIR\/chr22_gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf
```

The lifted over VCF will be saved to the following path:
```
imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf
```
\
**3.** After lifting over and reformatting, the European and African ancestry VCFs also need to be bgzipped and indexed:
```
bgzip imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf

bgzip imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf

bcftools index -t imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf.gz

bcftools index -t imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf.gz
```
\
**4.** We can then use ```bcftools``` to merge the ancestry-stratified VCF files into a single file:
```
bcftools reheader -h imputation/eur_correct_header.txt -o imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.revcomp.vcf.gz imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf.gz

bcftools reheader -h imputation/afr_correct_header.txt -o imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.revcomp.vcf.gz imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.vcf.gz

bcftools index -t imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.revcomp.vcf.gz

bcftools index -t imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.revcomp.vcf.gz

bcftools merge \
imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.withchr.hg38.filt.reheader.revcomp.vcf.gz \
imputation/gnomad_afr_eur_mega_snps.hg19.nodup.remove_chm13.biallelic_geno.0.05_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.withchr.hg38.filt.reheader.revcomp.vcf.gz > \
imputation/EUR_AFR_gnomAD.QC.vcf

bgzip imputation/EUR_AFR_gnomAD.QC.vcf

bcftools index -t imputation/EUR_AFR_gnomAD.QC.vcf.gz
```
\
**5.** We should also identify any multi-allelic variants that need to be filtered out prior to imputation:
```
bcftools query -f'%ID\n' -i 'N_ALT>1' imputation/EUR_AFR_gnomAD.QC.vcf.gz > multiallelic_sites.txt
```

We should find that there is a single multi-allelic site:
```
chr21:37444120:C:T;chr21:37444120:C:A
```

Next, we can converted the VCF to plink format, removing the multi-allelic site:
```
plink2 \
--vcf imputation/EUR_AFR_gnomAD.QC.vcf.gz \
--make-bed \
--out imputation/EUR_AFR_gnomAD.QC \
--chr chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
--allow-extra-chr \
--exclude multiallelic_sites.txt
```

We should observe the following output:
```
PLINK v2.00a5.12LM 64-bit Intel (25 Jun 2024)  www.cog-genomics.org/plink/2.0/
(C) 2005-2024 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to imputation/EUR_AFR_gnomAD.QC.log.
Options in effect:
  --allow-extra-chr
  --chr chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
  --exclude multiallelic_sites.txt
  --make-bed
  --out imputation/EUR_AFR_gnomAD.QC
  --vcf imputation/EUR_AFR_gnomAD.QC.vcf.gz

Start time: Thu Jan  9 08:49:03 2025
64224 MiB RAM detected, ~61884 available; reserving 32112 MiB for main
workspace.
Using up to 12 threads (change this with --threads).
--vcf: 865984 variants scanned.
--vcf: imputation/EUR_AFR_gnomAD.QC-temporary.pgen +
imputation/EUR_AFR_gnomAD.QC-temporary.pvar.zst +
imputation/EUR_AFR_gnomAD.QC-temporary.psam written.
1327 samples (0 females, 0 males, 1327 ambiguous; 1327 founders) loaded from
imputation/EUR_AFR_gnomAD.QC-temporary.psam.
865984 variants loaded from imputation/EUR_AFR_gnomAD.QC-temporary.pvar.zst.
Note: No phenotype data present.
--exclude: 865983 variants remaining.
865983 variants remaining after main filters.
Writing imputation/EUR_AFR_gnomAD.QC.fam ... done.
Writing imputation/EUR_AFR_gnomAD.QC.bim ... done.
Writing imputation/EUR_AFR_gnomAD.QC.bed ... done.
```
\
**6.** We will next prepare our dataset for imputation on the TOPMed Imputation Server. As of right now (January 2025), the TOPMed server supports a maximum sample size of 25,000 per dataset. If your dataset exceeds 25,000 samples, you will need to split it into multiple imputation batches.

As our dataset contains only 1,327 individuals, we will not run into this issue. However, a set of sample scripts can be found in the ```imputation``` sub-directory that demonstrates how this could be done. 

The developers of the TOPMed Imputation server provide a Perl script ```TOPMED-check-bim.pl``` that should be run to perform a series of pre-imputation checks on your data to verify that it is in the format expected by the server. This script can be downloaded from https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/ and can also be found in the ```imputation``` sub-directory.

Before running the ```TOPMED-check-bim.pl``` script, we need to create a .freq file of allele frequencies using ```plink```:
```
plink \
--bfile imputation/EUR_AFR_gnomAD.QC \
--freq \
--out imputation/EUR_AFR_gnomAD.QC
```

Next, we can run the ```TOPMED-check-bim.pl``` script. The ```topmed.ref.new.gz``` file is a list of SNPs in the TOPMed imputation panel and can be downloaded here:
```
perl TOPMED-check-bim.pl \
-b imputation/EUR_AFR_gnomAD.QC.bim \
-f imputation/EUR_AFR_gnomAD.QC.frq \
-r topmed.ref.new.gz \
-h \
-n
```

Once this script completes running, it will output a shell script called ```Run-plink.sh```. We can next execute this script to filter our dataset so that it is compatible with the TOPMed reference panel:
```
chmod 777 imputation/Run-plink.sh

./imputation/Run-plink.sh
```
\
**7.** In the resulting VCF files, we need to convert all of the chromosome names so that they have a ```chr``` prefix:
```
for vcf in $(ls imputation/*.vcf); do
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $vcf > $vcf\.withchr.vcf
done
```
\
**8.** Finally, all of the final VCF files can be bgzipped and indexed:
```
for vcf in $(ls imputation/*.withchr.vcf); do
bgzip $vcf
done

for vcf in $(ls imputation/*.withchr.vcf.gz); do
bcftools index -t $vcf
done
```
\
**9.** After compression and indexing is completed, all files should be downloaded locally.

To begin imputation, navigate to the TOPMed Imputation Server at https://imputation.biodatacatalyst.nhlbi.nih.gov and create an account, if you do not have one already.

On the home page, click "Run" from the upper menu bar and select "Genotype Imputation (Minimac4)" from the dropdown menu:
![alt text](https://github.com/mjbetti/genotype_qc/blob/main/topmed_home.png?raw=true)

Next, select the VCF files and indices to upload to the server. We will use the following parameters to run the imputation:
![alt text](https://github.com/mjbetti/genotype_qc/blob/main/topmed_submission.png?raw=true)

\
**10.** Upon completion of imputation, we will download and store it in a new sub-folder in our working directory:
```
mkdir imputation/imputed
```

We will unzip each of the imputed chromosome file, using the password provided by the TOPMed server:
```
for chr in {1..22}; do
unzip chr_$chr\.zip
done
```

## Post-imputation QC

The first post-imputation QC step is to filter out variants with an INFO score < 0.4. The INFO score can range from 0-1 and represents the confidence in each imputed SNP. In this tutorial, we only have a single imputation batch, so we can perform simple filtering based on INFO score.

If you are working with multiple batches, however, you should instead each set of INFO scores from each of the batches, averaged them for each SNP, and then made a list of SNPs that will be retained based on mean INFO score. An example of such a workflow can be found in the ```post-imputation``` sub-directory.

**1.** We will make a new directory for post-imputation QC:
```
mkdir imputation/imputed/qc
```
\
**2.** Next, we will write a shell script to concatenate INFO scores across all chromosomes:
```
nano gnomad_eur_afr_concat.sh
```
The following contents will be added to this script:
```
#!/bin/bash

zcat imputation/imputed/chr1.info.gz > imputation/imputed/qc/gnomad_eur_afr.info.txt

for i in {2..22}; do
    zcat imputation/imputed/batch1/chr$i.info.gz | tail -n +2 >> imputation/imputed/qc/gnomad_eur_afr.info.txt
done
```
This script can be run as a SLURM batch submission:
```
sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=60G \
    --time=1-00:00:00 \
    --job-name=batch1 \
    --wrap="./gnomad_eur_afr_concat.sh"
```
\
**3.** Using an ```R``` script, we can open the set of concatenated SNP INFO scores and compile a set of SNPs to keep based on an INFO score threshold > 0.4:
```
library("data.table")

path_batch1 <- "imputation/imputed/qc/gnomad_eur_afr.info.txt"

out_dir <- "imputation/imputed/qc"

#Open the info files as data frames
df_batch1 <- read_delim(path_batch1, delim = "\t", col_names = TRUE)
df_batch1 <- as.data.frame(df_batch1)
df_batch1 <- df_batch1[,c("SNP", "Rsq")]
names(df_batch1) <- c("SNP", "INFO")

final_df <- df_batch1[(df_batch1$INFO > 0.4),]
final_df <- final_df[,1]

write.table(final_df, file = paste(out_dir, "gnomad_eur_afr_info_cutoff_0.4.extract", sep = "/"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
\
**4.** Next, we can extract these SNPs from the imputed data:
```
for i in {1..22}; do
sbatch \
--job-name=chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=60G \
--time=2-00:00:00 \
--wrap="plink2 \
--vcf imputation/imputed/chr$i\.dose.vcf.gz \
--extract imputation/imputed/qc/gnomad_eur_afr_info_cutoff_0.4.extract \
--make-bed \
--out imputation/imputed/qc/chr$i\.info_0.4_gnomad_eur_afr"
done
```
\
**5.** After filtering SNPs based on INFO score, we will also remove any multi-allelic SNPs:
```
module load PLINK/1.9b_5.2

for i in {1..22}; do
sbatch \
--job-name=chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=60G \
--time=1-00:00:00 \
--wrap="plink \
--bfile imputation/imputed/qc/chr$i\.info_0.4_gnomad_eur_afr \
--snps-only \
--make-bed \
--out imputation/imputed/qc/chr$i\.info_0.4_gnomad_eur_afr_biallelic"
done
```

The resulting plink bed/bim/fam files should be saved to the following path:
```
imputation/imputed/qc/chr$i\.info_0.4_gnomad_eur_afr_biallelic
```
\
**6.** We imputed all gnomAD individuals together in a single batch, irrespective of genetic ancestry. Now, we will split these individuals back into groups of European and African ancestry based on the PCs we computed during the pre-imputation QC.

First, new directories were made to output the European and African ancestry individuals, respectively:
```
plink_dir=imputation/imputed/qc

mkdir $plink_dir\/eur
mkdir $plink_dir\/afr
```
\
**7.** Next, the keep files that we previously used to split individuals into African and European ancestry groups were reformatted to align with the FID and IID formatting in the imputed data:
```
keep_path_eur <- "pcs.txt.rescaled_pca.resubsetted.EUR.assign.refilter.no.refs.keep.txt"
keep_path_afr <- "pcs.txt.rescaled_pca.resubsetted.AFR.assign.refilter.no.refs.keep.txt"

keep_file_eur <- read.table(keep_path_eur, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
keep_df_eur <- as.data.frame(keep_file_eur)

keep_file_afr <- read.table(keep_path_afr, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
keep_df_afr <- as.data.frame(keep_file_afr)

keep_df_eur[,2] <- paste(keep_df_eur[,1], keep_df_eur[,2], keep_df_eur[,2], sep = "_")
keep_df_afr[,2] <- paste(keep_df_afr[,1], keep_df_afr[,2], keep_df_afr[,2], sep = "_")

write.table(keep_df_eur, file = "imputation/imputed/qc/eur/pcs.txt.rescaled_pca.resubsetted.EUR.assign.refilter.no.refs.keep.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(keep_df_afr, file = "imputation/imputed/qc/afr/pcs.txt.rescaled_pca.resubsetted.AFR.assign.refilter.no.refs.keep.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
These new keep files can then be used to split the imputed dataset into separate sets of European ancestry and African ancestry individuals:
### European ancestry
```
for i in {1..22}; do
sbatch \
--job-name=eur_chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink2 \
--bfile imputation/imputed/qc/chr$i\.info_0.4_gnomad_eur_afr_biallelic \
--keep imputation/imputed/qc/eur/pcs.txt.rescaled_pca.resubsetted.EUR.assign.refilter.no.refs.keep.txt \
--make-bed \
--out imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic"
done
```
### African ancestry
```
for i in {1..22}; do
sbatch \
--job-name=afr_chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink2 \
--bfile imputation/imputed/qc/chr$i\.info_0.4_gnomad_eur_afr_biallelic \
--keep imputation/imputed/qc/afr/pcs.txt.rescaled_pca.resubsetted.AFR.assign.refilter.no.refs.keep.txt \
--make-bed \
--out imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic"
done
```
\
**8.** Next, SNPs in each ancestry-stratified dataset with a MAF < 0.01 were removed:

**European ancestry**
```
for i in {1..22}; do
sbatch \
--job-name=afr_chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink \
--bfile imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic \
--maf 0.01 \
--make-bed \
--out imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01"
done
```

**African ancestry**
```
for i in {1..22}; do
sbatch \
--job-name=afr_chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink \
--bfile imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic \
--maf 0.01 \
--make-bed \
--out imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01"
done
```

The outputs were saved to the following paths:
```
imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01

imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01
```
\
**9.** Next, each set of SNPs should be pruned based on deviation from HWE. As we did in the pre-imputation QC, a HWE p-value threshold of 1 x 10<sup>-8</sup> will be used for the European ancestry individuals, and a threshold of 1 x 10<sup>-6</sup> will be used for African ancestry individuals:

**European ancestry**
```
for i in {1..22}; do
sbatch \
--job-name=eur_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink \
--bfile imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01 \
--hwe 1e-8 \
--make-bed \
--out imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-8"
done
```

**African ancestry**
```
for i in {1..22}; do
sbatch \
--job-name=afr_chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink \
--bfile imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01 \
--hwe 1e-6 \
--make-bed \
--out imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-6"
done
```

The outputs were saved to the following paths:
```
imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-8

imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-6
```
\
**10.** In the ```post-imputation``` sub-directory, we provide a ```python``` script called ```compute_filtered_snps_post_imputation.py``` that can compute the number of SNPs remaining after each post-imputation filtering step. This script takes in user arguments and can be run from the command line using the following shell script:
```
python compute_filtered_snps_post_imputation.py \
-d imputation/imputed/qc #the input file directory
-ii imputation/imputed/qc/gnomad_eur_afr.info.txt #the name of the input post-imputation SNP list with info scores \
-is imputation/imputed/qc/gnomad_eur_afr_info_cutoff_0.4.extract #the name of the input SNP list after filtering by info score \
-gi gnomad_eur_afr #the root name of the input genotype files
-o imputation/imputed/qc/gnomad_post_imputation_qc_stats.txt #the path of the output file
-v #print out logging to the terminal
```

We should observe the following output:

**11.** Next, we need to filter SNPs with MAF deviating > 0.1 from the 1000 Genomes MAF in the same target population. A detailed description of how reference allele frequencies from the 1000 Genomes dataset were prepared can be found in the ```reference_individuals_1kg/allele_frequencies_hg38``` sub-directory. 

We will create a new sub-folder in which this QC step will be performed:
```
mkdir imputation/imputed/qc/freq
```

Using the filtered post-imputation gnomAD dataset, we will compute allele frequencies in European and African ancestry individuals:

**European ancestry**
```
for i in {1..22}; do
sbatch \
--job-name=EUR_chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="plink2 --bfile imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-8 --freq --out imputation/imputed/qc/eur/EUR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-8"
done
```

**African ancestry**
```
for i in {1..22}; do
sbatch \
--job-name=AFR_chr$i \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="plink2 --bfile imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-6 --freq --out imputation/imputed/qc/afr/AFR_chr$i\.info_0.4_gnomad_eur_afr_biallelic_maf.0.01_hwe.1e-6"
done
```

We will next write the following script to compute allele frequency difference across all SNPs shared by our gnomAD datasets and European and African ancestry individuals from the 1000 Genomes reference dataset. We will convert any 1000 Genomes allele frequencies greater than 0.5 to the corresponding MAF by calculating (1 - AF):
```
library("data.table")

##EUR
#Open each of the files as a data frame
concat_df <- data.frame()

for (i in seq(1, 22)) {
	print(paste0("Processing chr", i, "..."))
	file_source <- fread(paste0("/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_chr", i, ".r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.afreq"), header = TRUE, sep = "\t", quote = FALSE)
	df_source <- as.data.frame(file_source)
	names(df_source)[5] <- "MAF_SOURCE"

	file_reference <- fread(paste0("/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/ceu/CEU.chr", i, "_GRCh38.genotypes.20170504.afreq"), header = TRUE, sep = "\t", quote = FALSE)
	df_reference <- as.data.frame(file_reference)
	names(df_reference)[5] <- "ALT_FREQS_REFERENCE"
	df_reference$MAF_REFERENCE <- df_reference$ALT_FREQS_REFERENCE
	df_reference <- df_reference[!is.na(df_reference$MAF_REFERENCE),]
	df_reference[df_reference$MAF_REFERENCE > 0.5,"MAF_REFERENCE"] <- (1 - df_reference[df_reference$MAF_REFERENCE > 0.5,"MAF_REFERENCE"])

	#Merge the two together, retaining shared variants
	merged_df <- merge(df_source, df_reference, by = "ID")
	merged_df <- merged_df[,c("ID", "MAF_SOURCE", "MAF_REFERENCE")]

	#Calculate difference in allele frequencies between the two samples
	merged_df$MAF_DIFFERENCE <- abs(merged_df$MAF_SOURCE - 	merged_df$MAF_REFERENCE)

	#Concatenate with genome-wide data frame
	concat_df <- rbind(concat_df, merged_df)
}

#Write out to a .txt file
write.table(concat_df, file = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_ALL_freq_differences.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##AFR
#Open each of the files as a data frame
concat_df <- data.frame()

for (i in seq(1, 22)) {
	print(paste0("Processing chr", i, "..."))
	file_source <- fread(paste0("/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/AFR_chr", i, ".r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.afreq"), header = TRUE, sep = "\t", quote = FALSE)
	df_source <- as.data.frame(file_source)
	names(df_source)[5] <- "MAF_SOURCE"

	file_reference <- fread(paste0("/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/asw/ASW.chr", i, "_GRCh38.genotypes.20170504.afreq"), header = TRUE, sep = "\t", quote = FALSE)
	df_reference <- as.data.frame(file_reference)
	names(df_reference)[5] <- "ALT_FREQS_REFERENCE"
	df_reference$MAF_REFERENCE <- df_reference$ALT_FREQS_REFERENCE
	df_reference <- df_reference[!is.na(df_reference$MAF_REFERENCE),]
	df_reference[df_reference$MAF_REFERENCE > 0.5,"MAF_REFERENCE"] <- (1 - df_reference[df_reference$MAF_REFERENCE > 0.5,"MAF_REFERENCE"])

	#Merge the two together, retaining shared variants
	merged_df <- merge(df_source, df_reference, by = "ID")
	merged_df <- merged_df[,c("ID", "MAF_SOURCE", "MAF_REFERENCE")]

	#Calculate difference in allele frequencies between the two samples
	merged_df$MAF_DIFFERENCE <- abs(merged_df$MAF_SOURCE - 	merged_df$MAF_REFERENCE)

	#Concatenate with genome-wide data frame
	concat_df <- rbind(concat_df, merged_df)
}

#Write out to a .txt file
write.table(concat_df, file = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/AFR_ALL_freq_differences.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
```

Next, the following script was written to output files of those variants whose frequencies deviate by more than 0.1 and 0.2, as well as a plot of the concordance between datasets:


```
library("optparse")
library("data.table")
library("ggplot2")
library("ggpmisc")

option_list = list(
	make_option(c("-s", "--path_source"), type = "character", default = NULL, help = "path of concated genotype frequency file", metavar = "character"),
	make_option(c("-o", "--output_prefix"), type = "character", default = NULL, help = "path of the output files and their prefix file name for all output files")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Open each of the files as a data frame
file_source <- fread(opt$path_source, header = TRUE, sep = "\t", quote = FALSE)
df_source <- as.data.frame(file_source)

#Plot the conccordance between the two sets of allele frequencies
pdf(paste0(opt$output_prefix, ".pdf"))
ggplot(df_source, aes(MAF_REFERENCE, MAF_SOURCE)) +
	geom_point(color='black') +
	geom_smooth(method='lm', formula = y ~ x) +
	stat_poly_eq(formula = y ~ x)
dev.off()

#Find variants with high differences in allele frequency
diff_0.1 <- df_source[(df_source$MAF_DIFFERENCE > 0.1),]
diff_0.2 <- df_source[(df_source$MAF_DIFFERENCE > 0.2),]

#Write these out to .txt files
write.table(diff_0.1, file = paste0(opt$output_prefix, "_diff0.1.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(diff_0.2, file = paste0(opt$output_prefix, "_diff0.2.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
```

This script was run for both the AFR and EUR samples:
```
#AFR
Rscript /home/bettimj/aldrich_rotation/sccs_aalc_inhale_joint_qc/freq/maf_concordance_with_1000_genomes.R \
-s /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/AFR_ALL_freq_differences.txt \
-o /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/AFR_ALL_freq_differences

#EUR
Rscript /home/bettimj/aldrich_rotation/sccs_aalc_inhale_joint_qc/freq/maf_concordance_with_1000_genomes.R \
-s /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_ALL_freq_differences.txt \
-o /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_ALL_freq_differences
```

The following were the counts of variants:
```
wc -l *.txt
     93984 AFR_ALL_freq_differences_diff0.1.txt
      3357 AFR_ALL_freq_differences_diff0.2.txt
  13732969 AFR_ALL_freq_differences.txt
     16483 EUR_ALL_freq_differences_diff0.1.txt
      2387 EUR_ALL_freq_differences_diff0.2.txt
   7979385 EUR_ALL_freq_differences.txt
  21828565 total
```

...and the following were the two plots:
AFR

EUR

Relatedness testing prior to GWAS

Each of the AFR and EUR datasets was first concatenated together (all chromosomes combined):
```
#AFR
cd /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr

nano afr_allfiles.txt

AFR_chr1.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr2.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr3.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr4.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr5.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr6.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr7.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr8.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr9.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr10.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr11.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr12.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr13.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr14.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr15.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr16.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr17.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr18.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr19.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr20.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr21.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6
AFR_chr22.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6

#EUR
cd /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur

nano eur_allfiles.txt

EUR_chr1.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr2.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr3.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr4.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr5.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr6.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr7.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr8.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr9.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr10.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr11.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr12.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr13.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr14.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr15.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr16.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr17.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr18.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr19.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr20.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr21.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
EUR_chr22.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
```

```
module load PLINK/1.9b_5.2

#AFR
cd /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr

plink --merge-list afr_allfiles.txt --make-bed --out AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6

#EUR
cd /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur

plink --merge-list eur_allfiles.txt --make-bed --out EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
```

Reformatted the fam files of each set of merged files:
```
afr_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.fam"
eur_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.fam"

#AFR
afr_file <- read.table(afr_path, header = FALSE, sep = " ", stringsAsFactors = FALSE)
afr_df <- as.data.frame(afr_file)
afr_ids <- strsplit(afr_df[,2], "_")
afr_ids <- unlist(lapply(afr_ids, `[[`, 3))
afr_df[,c(1:2)] <- afr_ids

write.table(afr_df, file = afr_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#EUR
eur_file <- read.table(eur_path, header = FALSE, sep = " ", stringsAsFactors = FALSE)
eur_df <- as.data.frame(eur_file)
eur_ids <- strsplit(eur_df[,2], "_")
eur_ids <- unlist(lapply(eur_ids, `[[`, 3))
eur_df[,c(1:2)] <- eur_ids

write.table(eur_df, file = eur_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

...and then ran KING:
```
mkdir /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king
cd /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king

mkdir afr eur
```

AFR
```
#AFR
cd /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr

king \
-b /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.bed \
--related \
--degree 3 \
--cpus 64
```

Wrote the following script to compile the kinship results with IBS-determined relationships:
```
library("data.table")

kinship_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/king.kin0"

kinship_file <- fread(kinship_path, header = TRUE, sep = "\t", quote = "")
kinship_df <- as.data.frame(kinship_file)

#ibs_file <- fread(ibs_path, header = TRUE, sep = "\t", quote = "")
#ibs_df <- as.data.frame(ibs_file)

#Merge the two files together
#merged_df <- merge(kinship_df, ibs_df, by = c("FID1", "ID1", "FID2", "ID2"))

#Retain only those with a kinship coefficient > 0.088, which would be the geometric mean of third-degree relatives
kinship_df <- kinship_df[(kinship_df$Kinship > 0.088),]

#Add familial relationships based on the parameters that Melinda sent
#Sort by kinship coefficient
kinship_df <- kinship_df[order(kinship_df$Kinship),]
#ibs_df$predicted_relationship <- NA
#ibs_df[((ibs_df$N_IBS1 / ibs_df$N_SNP) > 0.8),ncol(ibs_df)] <- "parent_offspring"
#ibs_df[((ibs_df$N_IBS0 / ibs_df$N_SNP) > 0.1 & (ibs_df$N_IBS0 / #ibs_df$N_SNP) < 0.38),ncol(ibs_df)] <- "full_sibs"
#ibs_df[((ibs_df$N_IBS0 / ibs_df$N_SNP) > 0.38 & (ibs_df$N_IBS0 / ibs_df$N_SNP) < 0.6),ncol(ibs_df)] <- "half_sibs"
#ibs_df[((ibs_df$N_IBS1 / ibs_df$N_SNP) < 0.1 & (ibs_df$N_IBS0 / ibs_df$N_SNP) < 0.1),ncol(ibs_df)] <- "duplicates"

write.table(kinship_df, file = "afr_biovu_king_related_3_degree.sorted.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
```

Print out the individual missingness of each individual in the dataset:
```
module load PLINK/1.9b_5.2
plink \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6 \
--missing
```

Next, based on the number of missingness, the a list was compiled of the duplicate individuals to be removed. Of each of the duplicate individuals, the one with a lower SNP missing rate was retained:
```
library("data.table")

miss_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/plink.imiss"
miss_file <- fread(miss_path, header = TRUE, sep = " ", quote = "")
miss_df <- as.data.frame(miss_file)

rel_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/afr_biovu_king_related_3_degree.sorted.txt"
rel_file <- fread(rel_path, header = TRUE, sep = "\t", quote = "")
rel_df <- as.data.frame(rel_file)

#df_dup <- df[duplicated(df$IID),]

new_df <- data.frame(matrix(NA, nrow = 0, ncol = 2))
names(new_df) <- c("FID", "IID")
#unique_names <- unique(df_dup$IID)

to_remove <- c()

for (row in seq(1, nrow(rel_df))) {
    ind1 <- miss_df[(miss_df$IID == rel_df[row,2]),]
    ind2 <- miss_df[(miss_df$IID == rel_df[row,4]),]
    if (ind1[6] > ind2[6]) {
        missing_ind <- ind1[2]
    } else {
        missing_ind <- ind2[2]
    }
    print(missing_ind) 
    to_remove <- c(to_remove, missing_ind)
}

to_remove <- data.frame(0, unlist(to_remove))
to_remove <- unique(to_remove)

#Write to remove list
write.table(to_remove, file = "AFR_biovu.remove", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
```

Remove the related individuals using the related list:
```
plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6 \
--remove /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/AFR_biovu.remove \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.unrelated
```

So a total of 787 individuals were removed.

Excluded SNPs with MAF deviating from 1000 Genomes ASW by more than 0.2 using plink:
```
library("data.table")

file <- fread("/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/AFR_ALL_freq_differences_diff0.2.txt", sep = "\t", header = TRUE, quote = "")

df <- as.data.frame(file)
df <- df[,1]

write.table(df, file = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_ALL_freq_differences_diff0.2.exclude", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

```
plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.unrelated \
--exclude /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_ALL_freq_differences_diff0.2.exclude \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.unrelated.remove_maf0.2_deviation
```

The resulting files were saved at the following path:
```
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/afr/AFR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.unrelated.remove_maf0.2_deviation
```

EUR
```
#EUR
cd /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur

king \
-b /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.bed \
--relatedness \
--degree 3 \
--cpus 64
```

Wrote the following script to compile the kinship results with IBS-determined relationships:
```
library("data.table")

kinship_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/king.kin0"

kinship_file <- fread(kinship_path, header = TRUE, sep = "\t", quote = "")
kinship_df <- as.data.frame(kinship_file)

#ibs_file <- fread(ibs_path, header = TRUE, sep = "\t", quote = "")
#ibs_df <- as.data.frame(ibs_file)

#Merge the two files together
#merged_df <- merge(kinship_df, ibs_df, by = c("FID1", "ID1", "FID2", "ID2"))

#Retain only those with a kinship coefficient > 0.088, which would be the geometric mean of third-degree relatives
kinship_df <- kinship_df[(kinship_df$Kinship > 0.088),]

#Add familial relationships based on the parameters that Melinda sent
#Sort by kinship coefficient
kinship_df <- kinship_df[order(kinship_df$Kinship),]
#ibs_df$predicted_relationship <- NA
#ibs_df[((ibs_df$N_IBS1 / ibs_df$N_SNP) > 0.8),ncol(ibs_df)] <- "parent_offspring"
#ibs_df[((ibs_df$N_IBS0 / ibs_df$N_SNP) > 0.1 & (ibs_df$N_IBS0 / #ibs_df$N_SNP) < 0.38),ncol(ibs_df)] <- "full_sibs"
#ibs_df[((ibs_df$N_IBS0 / ibs_df$N_SNP) > 0.38 & (ibs_df$N_IBS0 / ibs_df$N_SNP) < 0.6),ncol(ibs_df)] <- "half_sibs"
#ibs_df[((ibs_df$N_IBS1 / ibs_df$N_SNP) < 0.1 & (ibs_df$N_IBS0 / ibs_df$N_SNP) < 0.1),ncol(ibs_df)] <- "duplicates"

write.table(kinship_df, file = "eur_biovu_king_related_3_degree.sorted.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
```

Print out the individual missingness of each individual in the dataset:
```
module load PLINK/1.9b_5.2
plink \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8 \
--missing
```

Next, based on the number of missingness, the a list was compiled of the duplicate individuals to be removed. Of each of the duplicate individuals, the one with a lower SNP missing rate was retained:
```
library("data.table")

miss_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/plink.imiss"
miss_file <- fread(miss_path, header = TRUE, sep = " ", quote = "")
miss_df <- as.data.frame(miss_file)

rel_path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/eur_biovu_king_related_3_degree.sorted.txt"
rel_file <- fread(rel_path, header = TRUE, sep = "\t", quote = "")
rel_df <- as.data.frame(rel_file)

#df_dup <- df[duplicated(df$IID),]

new_df <- data.frame(matrix(NA, nrow = 0, ncol = 2))
names(new_df) <- c("FID", "IID")
#unique_names <- unique(df_dup$IID)

to_remove <- c()

for (row in seq(1, nrow(rel_df))) {
    ind1 <- miss_df[(miss_df$IID == rel_df[row,2]),]
    ind2 <- miss_df[(miss_df$IID == rel_df[row,4]),]
    if (ind1[6] > ind2[6]) {
        missing_ind <- ind1[2]
    } else {
        missing_ind <- ind2[2]
    }
    print(missing_ind) 
    to_remove <- c(to_remove, missing_ind)
}

to_remove <- data.frame(0, unlist(to_remove))
to_remove <- unique(to_remove)

#Write to remove list
write.table(to_remove, file = "EUR_biovu.remove", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
```

Remove the related individuals using the related list:
```
plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8 \
--remove /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/EUR_biovu.remove \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.unrelated
```

So a total of 2,199 individuals were removed.

Excluded SNPs with MAF deviating from 1000 Genomes CEU by more than 0.2 using plink:
```
library("data.table")

file <- fread("/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_ALL_freq_differences_diff0.2.txt", sep = "\t", header = TRUE, quote = "")

df <- as.data.frame(file)
df <- df[,1]

write.table(df, file = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_ALL_freq_differences_diff0.2.exclude", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

```
plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.unrelated \
--exclude /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_ALL_freq_differences_diff0.2.exclude \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.unrelated.remove_maf0.2_deviation
```

The resulting files were saved at the following path:
```
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/king/eur/EUR_biovu.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.unrelated.remove_maf0.2_deviation
```