#The tutorial found on https://www.gnxp.com/WordPress/2018/07/13/tutorial-to-run-supervised-admixture-analyses/?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+RazibKhansTotalFeed+%28Razib+Khan%27s+total+feed%29 was helpful in generating this workflow.

#Declare the source directory where the plink binary files (including fam annotated with founder populations)are located
SOURCE_DIR=/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno
TARGET_DIR=/home/bettimj/gamazon_rotation/1000_genomes/compile_plink/pop_anno
FAM_NAME=ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic

#European
#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from each of the 5 super-populations
grep "CEU\|TSI\|FIN\|GBR\|IBS" $SOURCE_DIR\/$FAM_NAME\.fam > $TARGET_DIR\/eur_keep.keep

#Use plink to generate an entire new set of binary files containing just the EUR 1000 Genomes individuals
plink2 --bfile $SOURCE_DIR\/$FAM_NAME \
--keep $TARGET_DIR\/eur_keep.keep \
--make-bed \
--out $TARGET_DIR\/eur_$FAM_NAME

#African
#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from each of the 5 super-populations
grep "YRI\|LWK\|GWD\|MSL\|ESN\|ASW\|ACB" $SOURCE_DIR\/$FAM_NAME.fam > $TARGET_DIR\/afr_keep.keep

#Use plink to generate an entire new set of binary files containing just the AFR 1000 Genomes individuals
plink2 --bfile $SOURCE_DIR\/$FAM_NAME \
--keep $TARGET_DIR\/afr_keep.keep \
--make-bed \
--out $TARGET_DIR\/afr_$FAM_NAME

#East Asian
#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from each of the 5 super-populations
grep "CHB\|JPT\|CHS\|CDX\|KHV" $SOURCE_DIR\/$FAM_NAME\.fam > $TARGET_DIR\/eas_keep.keep

#Use plink to generate an entire new set of binary files containing just the EAS 1000 Genomes individuals
plink2 --bfile $SOURCE_DIR\/$FAM_NAME \
--keep $TARGET_DIR\/eas_keep.keep \
--make-bed \
--out $TARGET_DIR\/eas_$FAM_NAME

#South Asian
#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from each of the 5 super-populations
grep "GIH\|PJL\|BEB\|STU\|ITU" $SOURCE_DIR\/$FAM_NAME\.fam > $TARGET_DIR\/sas_keep.keep

#Use plink to generate an entire new set of binary files containing just the SAS 1000 Genomes individuals
plink2 --bfile $SOURCE_DIR\/$FAM_NAME \
--keep $TARGET_DIR\/sas_keep.keep \
--make-bed \
--out $TARGET_DIR\/sas_$FAM_NAME

#Admixed American
#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from each of the 5 super-populations
grep "MXL\|PUR\|CLM\|PEL" $SOURCE_DIR\/$FAM_NAME\.fam > $TARGET_DIR\/amr_keep.keep

#Use plink to generate an entire new set of binary files containing just the SAS 1000 Genomes individuals
plink2 --bfile $SOURCE_DIR\/$FAM_NAME \
--keep $TARGET_DIR\/amr_keep.keep \
--make-bed \
--out $TARGET_DIR\/amr_$FAM_NAME




