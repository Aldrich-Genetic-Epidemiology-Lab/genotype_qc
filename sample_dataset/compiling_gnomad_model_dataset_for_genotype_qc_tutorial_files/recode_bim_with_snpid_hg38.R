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
varids <- paste(paste0("chr", bim_df[,1]), bim_df[,4], bim_df[,5], sep = "_")
bim_df[,2] <- varids

#Write this new fam file out to a fam file
out_file_name <- opt$output_path
write.table(bim_df, file = out_file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)