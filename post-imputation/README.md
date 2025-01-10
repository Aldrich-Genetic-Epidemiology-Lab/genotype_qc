```
mkdir /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc

#Batch 1
nano batch1_concat.sh

#!/bin/bash

zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch1/chr1.info.gz > /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch1.info.txt

for i in {1..22}; do
    zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch1/chr$i.info.gz | tail -n +2 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch1.info.txt
done

sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64G \
    --time=1-00:00:00 \
    --job-name=batch1 \
    --account=aldrich_lab \
    ./batch1_concat.sh

#Batch 2
nano batch2_concat.sh

#!/bin/bash

zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch2/chr1.info.gz > /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch2.info.txt

for i in {1..22}; do
    zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch2/chr$i.info.gz | tail -n +2 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch2.info.txt
done

sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64G \
    --time=1-00:00:00 \
    --job-name=batch2 \
    --account=aldrich_lab \
    ./batch2_concat.sh

#Batch 3
nano batch3_concat.sh

#!/bin/bash

zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch3/chr1.info.gz > /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch3.info.txt

for i in {1..22}; do
    zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch3/chr$i.info.gz | tail -n +2 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch3.info.txt
done

sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64G \
    --time=1-00:00:00 \
    --job-name=batch3 \
    --account=aldrich_lab \
    ./batch3_concat.sh

#Batch 4
nano batch4_concat.sh

#!/bin/bash

zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch4/chr1.info.gz > /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch4.info.txt

for i in {1..22}; do
    zcat /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch4/chr$i.info.gz | tail -n +2 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch4.info.txt
done

sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64G \
    --time=1-00:00:00 \
    --job-name=batch4 \
    --account=aldrich_lab \
    ./batch4_concat.sh
```

```
salloc \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=512G \
--account=aldrich_lab \
--time=2-00:00:00
```

```
library("data.table")

path_batch1 <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch1.info.txt"
path_batch2 <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch2.info.txt"
path_batch3 <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch3.info.txt"
path_batch4 <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/batch4.info.txt"

out_dir <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc"

#Open the info files as data frames
df_batch1 <- read_delim(path_batch1, delim = "\t", col_names = TRUE)
df_batch1 <- as.data.frame(df_batch1)
df_batch1 <- df_batch1[,c("SNP", "Rsq")]
names(df_batch1) <- c("SNP", "R2_batch1")

df_batch2 <- read_delim(path_batch2, delim = "\t", col_names = TRUE)
df_batch2 <- as.data.frame(df_batch2)
df_batch2 <- df_batch2[,c("SNP", "Rsq")]
names(df_batch2) <- c("SNP", "R2_batch2")

df_batch3 <- read_delim(path_batch3, delim = "\t", col_names = TRUE)
df_batch3 <- as.data.frame(df_batch3)
df_batch3 <- df_batch3[,c("SNP", "Rsq")]
names(df_batch3) <- c("SNP", "R2_batch3")

df_batch4 <- read_delim(path_batch4, delim = "\t", col_names = TRUE)
df_batch4 <- as.data.frame(df_batch4)
df_batch4 <- df_batch4[,c("SNP", "Rsq")]
names(df_batch4) <- c("SNP", "R2_batch4")

final_df <- data.frame()

for (i in seq(1:22)) {
print(i)
merged_df <- merge(df_batch1[startsWith(df_batch1$SNP, paste0("chr", i)),], df_batch2[startsWith(df_batch2$SNP, paste0("chr", i)),], by = "SNP")
merged_df <- merge(merged_df, df_batch3[startsWith(df_batch3$SNP, paste0("chr", i)),], by = "SNP")
merged_df <- merge(merged_df, df_batch4[startsWith(df_batch4$SNP, paste0("chr", i)),], by = "SNP")
mean_r2s <- rowMeans(merged_df[,2:ncol(merged_df)])
merged_df <- data.frame(merged_df$SNP, mean_r2s)
merged_df <- merged_df[(merged_df$mean_r2s > 0.4),]
final_df <- rbind(final_df, merged_df)
}

final_df <- final_df[,1]

write.table(final_df, file = "biovu_meanr2_cutoff_0.4.extract", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Next extracted these SNPs from the already-merged sets of imputed data:
```
#Batch 1
for i in {1..22}; do
sbatch \
--job-name=b1_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=64G \
--time=2-00:00:00 \
--wrap="plink2 \
--vcf /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch1/chr$i\.dose.vcf.gz \
--extract /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/biovu_meanr2_cutoff_0.4.extract \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch1"
done

#Batch 2
for i in {1..22}; do
sbatch \
--job-name=b2_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=64G \
--time=2-00:00:00 \
--wrap="plink2 \
--vcf /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch2/chr$i\.dose.vcf.gz \
--extract /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/biovu_meanr2_cutoff_0.4.extract \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch2"
done

#Batch 3
for i in {1..22}; do
sbatch \
--job-name=b3_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=64G \
--time=2-00:00:00 \
--wrap="plink2 \
--vcf /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch3/chr$i\.dose.vcf.gz \
--extract /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/biovu_meanr2_cutoff_0.4.extract \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch3"
done

#Batch 4
for i in {1..22}; do
sbatch \
--job-name=b4_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=64G \
--time=2-00:00:00 \
--wrap="plink2 \
--vcf /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/batch4/chr$i\.dose.vcf.gz \
--extract /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/biovu_meanr2_cutoff_0.4.extract \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch4"
done
```

The by-chromosome VCFs for each batch were next combined together:
```
module load PLINK/1.9b_5.2

for i in {1..22}; do
sbatch \
--job-name=b4_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=256G \
--time=2-00:00:00 \
--wrap="echo /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch1 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.to_merge.txt
echo /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch2 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.to_merge.txt
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch3 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.to_merge.txt
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_batch4 >> /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.to_merge.txt

plink --merge-list /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.to_merge.txt --make-bed --out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_merged_batches"
done
```

The resulting sets of filtered files were saved to the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_merged_batches
```

With the R2-filtered dataset, only bi-allelic SNPs were retained:
```
module load PLINK/1.9b_5.2

for i in {1..22}; do
sbatch \
--job-name=biovu_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=64G \
--time=1-00:00:00 \
--wrap="plink \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_merged_batches \
--snps-only \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_merged_batches_biallelic"
done
```

The outputs were saved to the following path:
```
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr$i\.r2_0.4_merged_batches_biallelic
```

Split up by genetic ancestry:
First reformatted the keep files to match the fam files:
```
plink_dir=/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc

mkdir $plink_dir\/afr
mkdir $plink_dir\/eur
```

```
keep_path_afr <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/script_agnostic_50/pcs.txt.rescaled_pca.resubsetted.AFR.assign.refilter.no.refs.keep.txt"
keep_path_eur <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/script_agnostic_50/pcs.txt.rescaled_pca.resubsetted.EUR.assign.refilter.no.refs.keep.txt"

keep_file_afr <- read.table(keep_path_afr, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
keep_df_afr <- as.data.frame(keep_file_afr)

keep_file_eur <- read.table(keep_path_eur, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
keep_df_eur <- as.data.frame(keep_file_eur)

keep_df_afr[,2] <- paste(keep_df_afr[,1], keep_df_afr[,2], keep_df_afr[,2], sep = "_")
keep_df_eur[,2] <- paste(keep_df_eur[,1], keep_df_eur[,2], keep_df_eur[,2], sep = "_")

write.table(keep_df_afr, file = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/pcs.txt.rescaled_pca.resubsetted.AFR.assign.refilter.no.refs.keep.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(keep_df_eur, file = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/pcs.txt.rescaled_pca.resubsetted.EUR.assign.refilter.no.refs.keep.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

```
#AFR
for i in {1..22}; do
sbatch \
--job-name=afr_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink2 \
--bfile $plink_dir\/chr$i\.r2_0.4_merged_batches_biallelic \
--keep /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/pcs.txt.rescaled_pca.resubsetted.AFR.assign.refilter.no.refs.keep.txt \
--make-bed \
--out $plink_dir\/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic"
done

#EUR
for i in {1..22}; do
sbatch \
--job-name=eur_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink2 \
--bfile $plink_dir\/chr$i\.r2_0.4_merged_batches_biallelic \
--keep /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/pcs.txt.rescaled_pca.resubsetted.EUR.assign.refilter.no.refs.keep.txt \
--make-bed \
--out $plink_dir\/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic"
done
```

Next, SNPs in each dataset with a MAF < 0.01 were removed:
```
module load PLINK/1.9b_5.2

#AFR
for i in {1..22}; do
sbatch \
--job-name=afr_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic \
--maf 0.01 \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01"
done

#EUR
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
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic \
--maf 0.01 \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01"
done
```

The outputs were saved to the following paths:
```
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01

/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01
```

Next, each set of SNPs was pruned based on deviation from HWE:
```
module load PLINK/1.9b_5.2

#AFR
for i in {1..22}; do
sbatch \
--job-name=afr_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=8G \
--time=2-00:00:00 \
--wrap="plink \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01 \
--hwe 1e-6 \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6"
done

#EUR
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
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01 \
--hwe 1e-8 \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8"
done
```

The outputs were saved to the following paths:
```
/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6

/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8
```

The following Python script was written to calculate all of the SNPs removed over the course of these past few filtering steps:


```
from dask import dataframe as dd
import dask.array as da
from dask.diagnostics import ProgressBar
ProgressBar().register()
import pandas as pd

##AFR
print("Importing merged file...")
merged_snps_afr = 0
#Merged
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr" + str(i) + ".r2_0.4_merged_batches.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	merged_snps_chr = len(df.index)
	merged_snps_afr += merged_snps_chr

print("Importing biallelic file...")
biallelic_snps_afr = 0
#Biallelic
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr" + str(i) + ".r2_0.4_merged_batches_biallelic.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	biallelic_snps_chr = len(df.index)
	biallelic_snps_afr += biallelic_snps_chr
	
print("Importing MAF file...")
maf_snps_afr = 0
#MAF
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr" + str(i) + ".r2_0.4_merged_batches_biallelic_maf.0.01.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	maf_snps_chr = len(df.index)
	maf_snps_afr += maf_snps_chr
	
print("Importing HWE file...")
hwe_snps_afr = 0
#HWE
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr" + str(i) + ".r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	hwe_snps_chr = len(df.index)
	hwe_snps_afr += hwe_snps_chr

##EUR
print("Importing merged file...")
merged_snps_eur = 0
#Merged
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr" + str(i) + ".r2_0.4_merged_batches.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	merged_snps_chr = len(df.index)
	merged_snps_eur += merged_snps_chr

print("Importing biallelic file...")
biallelic_snps_eur = 0
#Biallelic
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/chr" + str(i) + ".r2_0.4_merged_batches_biallelic.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	biallelic_snps_chr = len(df.index)
	biallelic_snps_eur += biallelic_snps_chr
	
print("Importing MAF file...")
maf_snps_eur = 0
#MAF
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr" + str(i) + ".r2_0.4_merged_batches_biallelic_maf.0.01.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	maf_snps_chr = len(df.index)
	maf_snps_eur += maf_snps_chr
	
print("Importing HWE file...")
hwe_snps_eur = 0
#HWE
for i in range(1, 23):
	print("Processing chr" + str(i) + "...")
	path = "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr" + str(i) + ".r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8.bim"
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	hwe_snps_chr = len(df.index)
	hwe_snps_eur += hwe_snps_chr
	
##Make output file
pops = ["AFR", "EUR"]
merged_snps = [merged_snps_afr, merged_snps_eur]
biallelic_snps = [biallelic_snps_afr, biallelic_snps_eur]
maf_snps = [maf_snps_afr, maf_snps_eur]
hwe_snps = [hwe_snps_afr, hwe_snps_eur]

out_df = pd.DataFrame(list(zip(pops, merged_snps, biallelic_snps, maf_snps, hwe_snps)),
               columns = ["pop", "merged", "biallelic", "maf", "hwe"])

out_df.to_csv("/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/biovu_merged_n_snps_filtered_via_biallelic_maf_hwe.txt", sep = "\t", header = True, index = False)
```

...and the following was the result:


Next, we need to filter SNPs with MAF deviating > 0.1 from the 1000 Genomes MAF in the same target population.

I previously downloaded a set of 1000 Genomes VCF files based on hg38 at the following path:
```
/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf
```

These were converted to plink format:
```
mkdir /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink

for i in {1..22}; do
sbatch \
--job-name=chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="plink2 \
--vcf /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/ALL.chr$i\_GRCh38.genotypes.20170504.vcf.gz \
--make-bed \
--vcf-half-call missing \
--max-alleles 2 \
--out /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/ALL.chr$i\_GRCh38.genotypes.20170504"
done
```

An Excel spreadsheet containing population information for all of the individuals in the 1000 Genomes database was downloaded from https://www.internationalgenome.org/faq/can-i-get-phenotype-gender-and-family-relationship-information-samples/.

This first sheet within this Excel spreadsheet was converted to a TSV file for easier use:


The following R script was then used to replace the family IDs in the original fam file (originally all written as 0) with each individual's corresponding population identifier, as recorded in the TSV file above:


```
library("optparse")

option_list = list(
	make_option(c("-f", "--fam_path"), type = "character", default = NULL, help = "path of input fam file", metavar = "character"),
	make_option(c("-o", "--output_path"), type = "character", default = NULL, help = "path of the output fam file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

pop_key_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/Sample_Info-Table_1.tsv"

#Open each of the files as a data frame
fam_file <- read.table(opt$fam_path, header = FALSE)
fam_df <- as.data.frame(fam_file)

pop_key_file <- read.delim(pop_key_path, header = TRUE)
pop_key_df <- as.data.frame(pop_key_file)

#For each match between sample IDs in the fam and key files, replace the family ID in the fam file with the 1000 Genomes population identifier from the key file
fam_ind_ids <- fam_df[,2]
pop_key_ind_ids <- pop_key_df[,1]

fam_df[,1] <- pop_key_df$Population[match(fam_ind_ids, pop_key_ind_ids)]

#Write this new fam file out to a fam file
out_file_name <- opt$output_path
write.table(fam_df, file = out_file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

Submitted script to run in parallel:
```
for i in {1..22}; do
sbatch \
--job-name=chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="Rscript /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/annotate_fam_with_pops_1kg_hg38.R \
-f /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/ALL.chr$i\_GRCh38.genotypes.20170504.fam \
-o /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/ALL.chr$i\_GRCh38.genotypes.20170504.fam"
done
```

Symbolic links of the original bed and bim files were created in the same directory as the new population-annotated fam file:
```
for i in {1..22}; do
ln -s /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/ALL.chr$i\_GRCh38.genotypes.20170504.bed /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno
done
```

The following script was written to re-code the variant IDs in the bim file:


```
library("optparse")

option_list = list(
	make_option(c("-b", "--bim_path"), type = "character", default = NULL, help = "path of input fam file", metavar = "character"),
	make_option(c("-o", "--output_path"), type = "character", default = NULL, help = "path of the output fam file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Open each of the files as a data frame
bim_file <- read.table(opt$bim_path, header = FALSE)
bim_df <- as.data.frame(bim_file)

#Replace the rsids in the bim file with the same format variant ID as in TOPMed
varids <- paste(paste0("chr", bim_df[,1]), bim_df[,4], bim_df[,6], bim_df[,5], sep = ":")
bim_df[,2] <- varids

#Write this new fam file out to a fam file
out_file_name <- opt$output_path
write.table(bim_df, file = out_file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

...and the following shell code was written to run this in parallel:
```
for i in {1..22}; do
sbatch \
--job-name=chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="Rscript /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/recode_bim_with_snpid_hg38.R \
-b /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/ALL.chr$i\_GRCh38.genotypes.20170504.bim \
-o /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/ALL.chr$i\_GRCh38.genotypes.20170504.bim"
done
```

So the following sets of 1000 Genomes plink binary files containing only individuals from the CEU or ASW populations was generated using the following script:
```
nano subset_ceu_asw_1kg.sh

CHR=${1}

module load PLINK/1.9b_5.2
#The tutorial found on https://www.gnxp.com/WordPress/2018/07/13/tutorial-to-run-supervised-admixture-analyses/?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+RazibKhansTotalFeed+%28Razib+Khan%27s+total+feed%29 was helpful in generating this workflow.

##CEU
mkdir /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/ceu

#Declare the source directory where the plink binary files (including fam annotated with founder populations)are located
source_dir=/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno
out_dir=/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/ceu

#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from CEU 
grep "CEU" $source_dir\/ALL.chr$CHR\_GRCh38.genotypes.20170504.fam > $out_dir\/keep_$CHR\.keep

#Use plink to generate an entire new set of binary files containing just CEU 1000 Genomes individuals
plink2 --bfile $source_dir\/ALL.chr$CHR\_GRCh38.genotypes.20170504 --keep $out_dir\/keep_$CHR\.keep --make-bed --allow-extra-chr --out $out_dir\/CEU.chr$CHR\_GRCh38.genotypes.20170504

##ASW
mkdir /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/asw

#Declare the source directory where the plink binary files (including fam annotated with founder populations)are located
source_dir=/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno
out_dir=/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/asw

#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from CEU 
grep "ASW" $source_dir\/ALL.chr$CHR\_GRCh38.genotypes.20170504.fam > $out_dir\/keep_$CHR\.keep

#Use plink to generate an entire new set of binary files containing just CEU 1000 Genomes individuals
plink2 --bfile $source_dir\/ALL.chr$CHR\_GRCh38.genotypes.20170504 --keep $out_dir\/keep_$CHR\.keep --make-bed --allow-extra-chr --out $out_dir\/ASW.chr$CHR\_GRCh38.genotypes.20170504
```

...and the script was run in parallel using SLURM:
```
for i in {1..22}; do
sbatch \
--job-name=chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/subset_ceu_asw_1kg.sh $i"
done
```

The following script was next written to calculate allele frequency for the CEU and ASW 1000 Genomes samples:
```
#CEU
for i in {1..22}; do
sbatch \
--job-name=CEU_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="plink2 --bfile /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/ceu/CEU.chr$i\_GRCh38.genotypes.20170504 --freq --out /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/ceu/CEU.chr$i\_GRCh38.genotypes.20170504"
done

#ASW
for i in {1..22}; do
sbatch \
--job-name=ASW_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="plink2 --bfile /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/asw/ASW.chr$i\_GRCh38.genotypes.20170504 --freq --out /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/plink/pop_anno/asw/ASW.chr$i\_GRCh38.genotypes.20170504"
done
```

The same allele frequency calculations were performed with the EUR and AFR imputed datasets:
```
mkdir /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq

#EUR
for i in {1..22}; do
sbatch \
--job-name=EUR_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="plink2 --bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/eur/EUR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8 --freq --out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/EUR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-8"
done

#AFR
for i in {1..22}; do
sbatch \
--job-name=AFR_chr$i \
--account=aldrich_lab \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem=16G \
--time=1-00:00:00 \
--wrap="plink2 --bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/afr/AFR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6 --freq --out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/imputed/qc/freq/AFR_chr$i\.r2_0.4_merged_batches_biallelic_maf.0.01_hwe.1e-6"
done
```

I will also convert the 1000 Genomes allele frequencies greater than 0.5 to MAF by calculating (1 - AF):


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
