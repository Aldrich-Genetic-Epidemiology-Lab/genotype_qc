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
