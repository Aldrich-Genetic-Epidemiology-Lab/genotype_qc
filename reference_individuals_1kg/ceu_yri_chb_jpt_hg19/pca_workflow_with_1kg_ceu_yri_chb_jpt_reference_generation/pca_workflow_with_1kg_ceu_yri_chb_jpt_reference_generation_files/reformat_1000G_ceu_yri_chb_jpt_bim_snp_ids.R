library("data.table")

#Declare the path of the bim file
bim_path <- "/home/bettimj/gamazon_rotation/1000_genomes/b37/pop_anno/ceu_yri_chb_jpt/ceu_yri_phase3_ceu_yri_chb_jpt.bim"

#Open the bim file as a data frame
bim_file <- fread(bim_path, header = FALSE, sep = "\t", quote = "")
bim_df <- as.data.frame(bim_file)

#Generate new variant IDs for each SNP using the coordinate and allele information in the other columns
bim_df[,2] <- paste(paste0("chr", bim_df[,1]), bim_df[,4], bim_df[,6], bim_df[,5], "b37", sep = "_")

#Write the modified bim out
write.table(bim_df, file = bim_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)