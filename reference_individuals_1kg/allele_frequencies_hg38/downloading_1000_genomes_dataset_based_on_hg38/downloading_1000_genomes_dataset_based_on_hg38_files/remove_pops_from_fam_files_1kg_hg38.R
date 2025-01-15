#Declare the paths of each super-population fam file
eur_fam_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/super_pop_files/eur_ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.fam"
afr_fam_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/super_pop_files/afr_ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.fam"
eas_fam_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/super_pop_files/eas_ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.fam"
sas_fam_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/super_pop_files/sas_ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.fam"
amr_fam_path <- "/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno/super_pop_files/amr_ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.fam"

#Open each of the files as a data frame
eur_fam_file <- read.table(eur_fam_path, header = FALSE)
eur_fam_df <- data.frame(eur_fam_file)

afr_fam_file <- read.table(afr_fam_path, header = FALSE)
afr_fam_df <- data.frame(afr_fam_file)

eas_fam_file <- read.table(eas_fam_path, header = FALSE)
eas_fam_df <- data.frame(eas_fam_file)

sas_fam_file <- read.table(sas_fam_path, header = FALSE)
sas_fam_df <- data.frame(sas_fam_file)

amr_fam_file <- read.table(amr_fam_path, header = FALSE)
amr_fam_df <- data.frame(amr_fam_file)

#Make new fam data frames, in which the first column containing the 1000 Genomes population identifiers is replaced with the individual identifier (column 2)
new_eur_fam_df <- data.frame(eur_fam_df[,2], eur_fam_df[,2], eur_fam_df[,3], eur_fam_df[,4], eur_fam_df[,5], eur_fam_df[,6])
write.table(new_eur_fam_df, file = eur_fam_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

new_afr_fam_df <- data.frame(afr_fam_df[,2], afr_fam_df[,2], afr_fam_df[,3], afr_fam_df[,4], afr_fam_df[,5], afr_fam_df[,6])
write.table(new_afr_fam_df, file = afr_fam_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

new_eas_fam_df <- data.frame(eas_fam_df[,2], eas_fam_df[,2], eas_fam_df[,3], eas_fam_df[,4], eas_fam_df[,5], eas_fam_df[,6])
write.table(new_eas_fam_df, file = eas_fam_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

new_sas_fam_df <- data.frame(sas_fam_df[,2], sas_fam_df[,2], sas_fam_df[,3], sas_fam_df[,4], sas_fam_df[,5], sas_fam_df[,6])
write.table(new_sas_fam_df, file = sas_fam_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

new_amr_fam_df <- data.frame(amr_fam_df[,2], amr_fam_df[,2], amr_fam_df[,3], amr_fam_df[,4], amr_fam_df[,5], amr_fam_df[,6])
write.table(new_amr_fam_df, file = amr_fam_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

