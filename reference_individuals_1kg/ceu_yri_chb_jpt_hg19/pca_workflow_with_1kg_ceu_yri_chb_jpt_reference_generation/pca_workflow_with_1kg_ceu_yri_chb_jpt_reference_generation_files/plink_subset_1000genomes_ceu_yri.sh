#The tutorial found on https://www.gnxp.com/WordPress/2018/07/13/tutorial-to-run-supervised-admixture-analyses/?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+RazibKhansTotalFeed+%28Razib+Khan%27s+total+feed%29 was helpful in generating this workflow.

#Declare the source directory where the plink binary files (including fam annotated with founder populations)are located
source_dir=/home/mbetti/aldrich_rotation/final_pop_fam
out_dir=/home/mbetti/aldrich_rotation/final_pop_fam/ceu_yri

#Using regular expressions, subset the individuals in the 1000 Genomes fam file that are from either the ASW or CEU 
grep "CEU\|YRI" $source_dir\/all_phase3.fam > $out_dir\/keep.keep

#Use plink to generate an entire new set of binary files containing just the ASW and CEU 1000 Genomes individuals
plink --bfile $source_dir\/all_phase3 --keep $out_dir\/keep.keep --make-bed --allow-extra-chr --out $out_dir\/ceu_yri_phase3
