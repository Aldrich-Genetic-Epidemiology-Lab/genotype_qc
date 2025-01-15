We need to fit 4 models:
- Model 1: Blood ~ genotype + age + sex + PC1-10
- Model 2: saliva ~ genotype + age + sex + PC1-10
- Model 3: mouthwash ~ genotype + age + sex + PC1-10
- Model 4: Batch ~ genotype + age + sex + PC1-10

Copied sample types files from Box to ACCRE:
```
scp /Users/michaelbetti/Library/CloudStorage/Box-Box/AldrichLab_Shared/Projects/LCDC_U01/SCCS\ Genotyping\ data/SCCS\ MEGA\ Batch1/SCCS\ ID\ list\ and\ matching\ vars/App379Req122bio_AddVar_06May2021.csv bettimj@vgi01.accre.vanderbilt.edu:/home/bettimj/aldrich_rotation

scp /Users/michaelbetti/Library/CloudStorage/Box-Box/AldrichLab_Shared/Projects/LCDC_U01/SCCS\ Genotyping\ data/SCCS\ MEGA\ Batch1/SCCS\ ID\ list\ and\ matching\ vars/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv bettimj@vgi01.accre.vanderbilt.edu:/home/bettimj/aldrich_rotation
```

...and the following path contains a file I previously made containing batch information for the samples:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/pca/batches/sccs_batches.txt
```

Also needed to generate covariate files for the AFR and EUR groups:

AFR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,5])
#df_1[(df_1[,2] == "M"),2] <- 1
#df_1[(df_1[,2] == "F"),2] <- 2
names(df_1) <- c("IDNumber", "Enrollment_Age")

#path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
#file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
#df_2 <- as.data.frame(file_2)
#df_2[,6] <- gsub("_", "-", df_2[,6])
#df_2 <- data.frame(df_2[,6], df_2[,1])
#df_2[(df_2[,2] == "M"),2] <- 1
#df_2[(df_2[,2] == "F"),2] <- 2
#names(df_2) <- c("IDNumber", "Sex")

sccs_ind_path <- "/home/bettimj/aldrich_rotation/data_harmonization/ind_study_data_20211027/379-Biospecimen-122.r2_App379Req122_harmon_07Dec2021_di_EJ.csv"
sccs_ind_file <- fread(sccs_ind_path, header = TRUE, sep = ",", quote = FALSE)
sccs_ind_df <- as.data.frame(sccs_ind_file)
sccs_ind_df[,1] <- gsub("_", "-", sccs_ind_df[,1])
sccs_ind_df <- data.frame(sccs_ind_df[,1], sccs_ind_df[,5])
names(sccs_ind_df) <- c("IDNumber", "Enrollment_Age")

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

df <- rbind(df_1, sccs_ind_df)
df <- unique(df)

df <- merge(fam_df, df, by.x = "V1", by.y = "IDNumber")

pc_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/pca/jackie_params/SCCS_eigenoutput0.pca.evec"
pc_file <- fread(pc_path, header = TRUE, sep = " ", quote = "")
pc_df <- as.data.frame(pc_file)
pc_df <- pc_df[,1:11]
names(pc_df) <- c("IDNumber", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

df <- merge(df, pc_df, by.x = "V1", by.y = "IDNumber")

df <- df[(df[,1] %in% fam_df[,1]),]
df <- data.frame(df[,1], df[,1], df[,c(5, 7:17)])

#Export as a phenotype file
write.table(df, file = "afr.cov", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

EUR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,5])
#df_1[(df_1[,2] == "M"),2] <- 1
#df_1[(df_1[,2] == "F"),2] <- 2
names(df_1) <- c("IDNumber", "Enrollment_Age")

#path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
#file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
#df_2 <- as.data.frame(file_2)
#df_2[,6] <- gsub("_", "-", df_2[,6])
#df_2 <- data.frame(df_2[,6], df_2[,1])
#df_2[(df_2[,2] == "M"),2] <- 1
#df_2[(df_2[,2] == "F"),2] <- 2
#names(df_2) <- c("IDNumber", "Sex")

sccs_ind_path <- "/home/bettimj/aldrich_rotation/data_harmonization/ind_study_data_20211027/379-Biospecimen-122.r2_App379Req122_harmon_07Dec2021_di_EJ.csv"
sccs_ind_file <- fread(sccs_ind_path, header = TRUE, sep = ",", quote = FALSE)
sccs_ind_df <- as.data.frame(sccs_ind_file)
sccs_ind_df[,1] <- gsub("_", "-", sccs_ind_df[,1])
sccs_ind_df <- data.frame(sccs_ind_df[,1], sccs_ind_df[,5])
names(sccs_ind_df) <- c("IDNumber", "Enrollment_Age")

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

df <- rbind(df_1, sccs_ind_df)
df <- unique(df)

df <- merge(fam_df, df, by.x = "V1", by.y = "IDNumber")

pc_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/pca/jackie_params/SCCS_eigenoutput0.pca.evec"
pc_file <- fread(pc_path, header = TRUE, sep = " ", quote = "")
pc_df <- as.data.frame(pc_file)
pc_df <- pc_df[,1:11]
names(pc_df) <- c("IDNumber", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

df <- merge(df, pc_df, by.x = "V1", by.y = "IDNumber")

df <- df[(df[,1] %in% fam_df[,1]),]
df <- data.frame(df[,1], df[,1], df[,c(5, 7:17)])

#Export as a phenotype file
write.table(df, file = "eur.cov", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

So first fit the model for blood:
The following script was written to generate a phenotype file that could be fed into plink:
Blood:
AFR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,2])
names(df_1) <- c("IDNumber", "SampleType")

path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
df_2 <- as.data.frame(file_2)
df_2[,6] <- gsub("_", "-", df_2[,6])
df_2 <- data.frame(df_2[,6], df_2[,7])
names(df_2) <- c("IDNumber", "SampleType")

df <- rbind(df_1, df_2)
df <- unique(df)

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

pheno_col <- df[,2]
pheno_col[pheno_col == "B"] <- 2
pheno_col[!(pheno_col == 2)] <- 1

phenotype_df <- data.frame(df[,1], df[,1], pheno_col)
phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "blood_afr.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/blood_afr.phe
```

Next run association analysis:
789 cases, 443 controls
```
plink2 \
--threads 6 \
--out sccs_afr_blood_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/blood_afr.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/afr.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_afr_blood_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 73 associated SNPs.
Using 1e-6, there are 3.
Using 5e-8, there are 3.

EUR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,2])
names(df_1) <- c("IDNumber", "SampleType")

path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
df_2 <- as.data.frame(file_2)
df_2[,6] <- gsub("_", "-", df_2[,6])
df_2 <- data.frame(df_2[,6], df_2[,7])
names(df_2) <- c("IDNumber", "SampleType")

df <- rbind(df_1, df_2)
df <- unique(df)

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

pheno_col <- df[,2]
pheno_col[pheno_col == "B"] <- 2
pheno_col[!(pheno_col == 2)] <- 1

phenotype_df <- data.frame(df[,1], df[,1], pheno_col)
phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "blood_eur.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/blood_eur.phe
```

Next run association analysis:
844 cases, 442 controls
```
plink2 \
--threads 6 \
--out sccs_eur_blood_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/blood_eur.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/eur.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_eur_blood_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 64 associated SNPs.
Using 1e-6, there are 1.
Using 5e-8, there are 0.

Mouthwash:
AFR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,2])
names(df_1) <- c("IDNumber", "SampleType")

path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
df_2 <- as.data.frame(file_2)
df_2[,6] <- gsub("_", "-", df_2[,6])
df_2 <- data.frame(df_2[,6], df_2[,7])
names(df_2) <- c("IDNumber", "SampleType")

df <- rbind(df_1, df_2)
df <- unique(df)

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

pheno_col <- df[,2]
pheno_col[pheno_col == "M"] <- 2
pheno_col[!(pheno_col == 2)] <- 1

phenotype_df <- data.frame(df[,1], df[,1], pheno_col)
phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "mouthwash_afr.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/mouthwash_afr.phe
```

Next run association analysis:
420 cases, 812 controls
```
plink2 \
--threads 6 \
--out sccs_afr_mouthwash_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/mouthwash_afr.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/afr.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_afr_mouthwash_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 90 associated SNPs.
Using 1e-6, there are 4.
Using 5e-8, there are 3.

EUR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,2])
names(df_1) <- c("IDNumber", "SampleType")

path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
df_2 <- as.data.frame(file_2)
df_2[,6] <- gsub("_", "-", df_2[,6])
df_2 <- data.frame(df_2[,6], df_2[,7])
names(df_2) <- c("IDNumber", "SampleType")

df <- rbind(df_1, df_2)
df <- unique(df)

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

pheno_col <- df[,2]
pheno_col[pheno_col == "M"] <- 2
pheno_col[!(pheno_col == 2)] <- 1

phenotype_df <- data.frame(df[,1], df[,1], pheno_col)
phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "mouthwash_eur.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/mouthwash_eur.phe
```

Next run association analysis:
388 cases, 898 controls
```
plink2 \
--threads 6 \
--out sccs_eur_mouthwash_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/mouthwash_eur.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/eur.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_eur_mouthwash_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 65 associated SNPs.
Using 1e-6, there are 2.
Using 5e-8, there are 1.

Saliva:
AFR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,2])
names(df_1) <- c("IDNumber", "SampleType")

path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
df_2 <- as.data.frame(file_2)
df_2[,6] <- gsub("_", "-", df_2[,6])
df_2 <- data.frame(df_2[,6], df_2[,7])
names(df_2) <- c("IDNumber", "SampleType")

df <- rbind(df_1, df_2)
df <- unique(df)

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

pheno_col <- df[,2]
pheno_col[pheno_col == "O"] <- 2
pheno_col[!(pheno_col == 2)] <- 1

phenotype_df <- data.frame(df[,1], df[,1], pheno_col)
phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "saliva_afr.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/saliva_afr.phe
```

Next run association analysis:
23 cases, 1209 controls
```
plink2 \
--threads 6 \
--out sccs_afr_saliva_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/saliva_afr.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/afr.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_afr_saliva_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 177 associated SNPs.
Using 1e-6, there are 2.
Using 5e-8, there are 0.

EUR:
```
library("data.table")

#Import the phenotype data
path_1 <- "/home/bettimj/aldrich_rotation/App379Req122bio_AddVar_06May2021.csv"
file_1 <- fread(path_1, header = TRUE, sep = ",", quote = FALSE)
df_1 <- as.data.frame(file_1)
df_1[,1] <- gsub("_", "-", df_1[,1])
df_1 <- data.frame(df_1[,1], df_1[,2])
names(df_1) <- c("IDNumber", "SampleType")

path_2 <- "/home/bettimj/aldrich_rotation/379-Biospecimen-122.r2_app379req122bio_forj_sampletype.csv"
file_2 <- fread(path_2, header = TRUE, sep = ",", quote = FALSE)
df_2 <- as.data.frame(file_2)
df_2[,6] <- gsub("_", "-", df_2[,6])
df_2 <- data.frame(df_2[,6], df_2[,7])
names(df_2) <- c("IDNumber", "SampleType")

df <- rbind(df_1, df_2)
df <- unique(df)

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = FALSE)
fam_df <- as.data.frame(fam_file)

pheno_col <- df[,2]
pheno_col[pheno_col == "O"] <- 2
pheno_col[!(pheno_col == 2)] <- 1

phenotype_df <- data.frame(df[,1], df[,1], pheno_col)
phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "saliva_eur.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/saliva_eur.phe
```

Next run association analysis:
54 cases, 1232 controls
```
plink2 \
--threads 6 \
--out sccs_eur_saliva_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/saliva_eur.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/eur.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_eur_saliva_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 126 associated SNPs.
Using 1e-6, there are 0.
Using 5e-8, there are 0.

Batch:
The following script was written to generate a phenotype file that could be fed into plink:
AFR:
```
library("data.table")

#Import the phenotype data
#Import the phenotype data
path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/pca/batches/sccs_batches.txt"
file <- fread(path, header = TRUE, sep = "\t", quote = FALSE)
df <- as.data.frame(file)
df[,1] <- gsub("_", "-", df[,1])

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = "")
fam_df <- as.data.frame(fam_file)

phenotype_df <- data.frame(df[,1], df[,1], df[,2])
phenotype_df[(phenotype_df[,3] == "batch1"),3] <- 1
phenotype_df[(phenotype_df[,3] == "batch2"),3] <- 2

phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "batch_afr.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/batch_afr.phe
```

Next run association analysis:
393 cases, 839 controls
```
plink2 \
--threads 6 \
--out sccs_afr_batch_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.afr_hwe.1e-6_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/batch_afr.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/afr.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_afr_batch_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 76 associated SNPs.
Using 1e-6, there are 3.
Using 5e-8, there are 1.

EUR:
```
library("data.table")

#Import the phenotype data
#Import the phenotype data
path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/pca/batches/sccs_batches.txt"
file <- fread(path, header = TRUE, sep = "\t", quote = FALSE)
df <- as.data.frame(file)
df[,1] <- gsub("_", "-", df[,1])

fam_path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet.fam"
fam_file <- fread(fam_path, header = FALSE, sep = "\t", quote = "")
fam_df <- as.data.frame(fam_file)

phenotype_df <- data.frame(df[,1], df[,1], df[,2])
phenotype_df[(phenotype_df[,3] == "batch1"),3] <- 1
phenotype_df[(phenotype_df[,3] == "batch2"),3] <- 2

phenotype_df <- phenotype_df[(phenotype_df[,1] %in% fam_df[,1]),]

#Export as a phenotype file
write.table(phenotype_df, file = "batch_eur.phe", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
```

The resulting file was stored at the following path on ACCRE:
```
/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/batch_eur.phe
```

Next run association analysis:
540 cases, 746 controls
```
plink2 \
--threads 6 \
--out sccs_eur_batch_gwas_sumstats \
--bfile /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/deliverable_MA-5702_BestPerformingSNPs_nodup_remove_hapmap_biallelic_geno.0.02_mind.0.05.eur_hwe.1e-8_maf.0.01_fhet \
--pheno /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/batch_eur.phe \
--covar /home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/eur.cov \
--mpheno 1 \
--chr 1-22 \
--logistic \
--covar-variance-standardize
```

Evaluated for associations:
```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/sccs_batches1and2_qc/post_pca_qc/sccs_eur_batch_gwas_sumstats.PHENO1.glm.logistic"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)
df <- df[(df$TEST == "ADD"),]
df <- df[!is.na(df$P),]
df_sig <- df[(df$P < 1e-4),]
nrow(df_sig)
df_sig <- df[(df$P < 1e-6),]
nrow(df_sig)
df_sig <- df[(df$P < 5e-8),]
nrow(df_sig)
```

Using 1e-4, there are 41 associated SNPs.
Using 1e-6, there are 0.
Using 5e-8, there are 0.
