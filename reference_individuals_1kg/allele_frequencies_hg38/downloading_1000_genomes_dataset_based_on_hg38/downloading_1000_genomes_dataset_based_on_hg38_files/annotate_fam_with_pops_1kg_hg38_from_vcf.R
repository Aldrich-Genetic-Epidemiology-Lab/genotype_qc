#Declare the path of the 1KG fam file, the reference TSV file with the corresponding populations of individuals in the cohort, and the target output directory
fam_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.fam"
pop_key_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/Sample_Info-Table_1.tsv"
out_dir <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno"

#Open each of the files as a data frame
fam_file <- read.table(fam_path, header = FALSE)
fam_df <- as.data.frame(fam_file)

pop_key_file <- read.delim(pop_key_path, header = TRUE)
pop_key_df <- as.data.frame(pop_key_file)

#For each match between sample IDs in the fam and key files, replace the family ID in the fam file with the 1000 Genomes population identifier from the key file
fam_ind_ids <- fam_df[,2]
pop_key_ind_ids <- pop_key_df[,1]

fam_df[,1] <- pop_key_df$Population[match(fam_ind_ids, pop_key_ind_ids)]

#Write this new fam file out to a fam file
out_file_name <- "ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.fam"
write.table(fam_df, file = paste(out_dir, out_file_name, sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)