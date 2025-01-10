```
library("data.table")

path <- "/home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC.fam"
file <- fread(path, header = FALSE, sep = "\t", quote = "")
df <- as.data.frame(file)

# Shuffle the dataframe
set.seed(123) # setting seed for reproducibility
df_shuffled <- df[sample(nrow(df)), ]

# Split the shuffled dataframe into batches
batch_size <- 25000

keep1 <- df_shuffled[c(1:batch_size), c("V1", "V2")]
keep2 <- df_shuffled[(batch_size+1):(2*batch_size), c("V1", "V2")]
keep3 <- df_shuffled[(2*batch_size+1):(3*batch_size), c("V1", "V2")]
keep4 <- df_shuffled[(3*batch_size+1):nrow(df_shuffled), c("V1", "V2")]

# Write each batch to a file
write.table(keep1, file = "batch1.keep", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(keep2, file = "batch2.keep", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(keep3, file = "batch3.keep", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(keep4, file = "batch4.keep", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Split using plink:
```
plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC \
--keep batch1.keep \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC.batch1

plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC \
--keep batch2.keep \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC.batch2

plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC \
--keep batch3.keep \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC.batch3

plink2 \
--bfile /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC \
--keep batch4.keep \
--make-bed \
--out /home/bettimj/aldrich_rotation/biovu_2023_pull_qc/imputation/AFR_EUR_BioVU.QC.batch4
```
