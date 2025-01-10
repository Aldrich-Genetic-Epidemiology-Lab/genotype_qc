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
