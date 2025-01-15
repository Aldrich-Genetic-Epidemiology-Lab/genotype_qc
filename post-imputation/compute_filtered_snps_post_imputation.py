import argparse
import os
#from dask import dataframe as dd
#import dask.array as da
#from dask.diagnostics import ProgressBar
#ProgressBar().register()
import pandas as pd

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-d", "--input_dir", type = str, required = True, help = "the input file directory (required)")
parser.add_argument("-ip", "--input_preimputation_file", type = str, required = True, help = "the root name of the input pre-imputation bim files (required)")
parser.add_argument("-ii", "--input_imputation_file", type = str, required = True, help = "the name of the input post-imputation SNP list with info scores (required)")
parser.add_argument("-is", "--input_info_score_file", type = str, required = True, help = "the name of the input SNP list after filtering by info score (required)")
parser.add_argument("-gi", "--input_genotype_file", type = str, required = True, help = "the root name of the input genotype files (required)")
parser.add_argument("-o", "--output", type = str, required = True, help = "the path of the output file")
parser.add_argument("-v", "--verbose", required = False, help = "return logging as terminal output", action = "store_true")

args = parser.parse_args()

##Pre-imputation datasets
#EUR
eur_pre_imputation_snps = 0

if args.verbose:
	print("Importing EUR pre-imputed SNPs...")
path = os.path.join(args.input_preimputation_file + ".eur_hwe.1e-8_maf.0.01_fhet.bim")

# Count lines that do not start with '#'
with open(path, 'r') as file:
    for line in file:
        if not line.startswith('#'):  # Skip lines starting with '#'
            eur_pre_imputation_snps += 1

if args.verbose:
    print(f"Total pre-imputation EUR SNPs: {eur_pre_imputation_snps}")
    print("\n")

#AFR
afr_pre_imputation_snps = 0

if args.verbose:
	print("Importing AFR pre-imputed SNPs...")
path = os.path.join(args.input_preimputation_file + ".afr_hwe.1e-6_maf.0.01_fhet.bim")

# Count lines that do not start with '#'
with open(path, 'r') as file:
    for line in file:
        if not line.startswith('#'):  # Skip lines starting with '#'
            afr_pre_imputation_snps += 1

if args.verbose:
    print(f"Total pre-imputation AFR SNPs: {afr_pre_imputation_snps}")
    print("\n")

##Combined dataset
#Post-imputation
post_imputation_snps = 0

if args.verbose:
	print("Importing imputed SNPs...")
path = os.path.join(args.input_dir, args.input_imputation_file)

# Count lines that do not start with '#'
with open(path, 'r') as file:
    for line in file:
        if not line.startswith('#'):  # Skip lines starting with '#'
            post_imputation_snps += 1

if args.verbose:
    print(f"Total post-imputation SNPs: {post_imputation_snps}")
    print("\n")

#Imputation R2 filtering
r2_snps = 0

if args.verbose:
	print("Importing remaining SNPs after imputation R2 filtering...")
path = os.path.join(args.input_dir, args.input_info_score_file)

# Count lines that do not start with '#'
with open(path, 'r') as file:
    for line in file:
        if not line.startswith('#'):  # Skip lines starting with '#'
            r2_snps += 1

if args.verbose:
    print(f"Total SNPs after imputation R2 filtering: {r2_snps}")
    print("\n")

#Biallelic SNPs
biallelic_snps = 0

if args.verbose:
    print("Importing biallelic file...")

for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
    path = os.path.join(args.input_dir, "chr" + str(i) + "." + args.input_genotype_file + ".r2_0.4_biallelic.bim")
    
    # Count lines that do not start with '#'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Skip lines starting with '#'
                biallelic_snps += 1

if args.verbose:
    print(f"Total SNPs after removing multi-allelic variants: {biallelic_snps}")
    print("\n")

##EUR
#Split into EUR ancestry
eur_biallelic_snps = 0

if args.verbose:
    print("Importing EUR biallelic file...")

for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
    path = os.path.join(args.input_dir, "eur", "EUR_chr" + str(i) + "." + args.input_genotype_file + ".r2_0.4_biallelic.bim")
    
    # Count lines that do not start with '#'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Skip lines starting with '#'
                eur_biallelic_snps += 1

if args.verbose:
    print(f"Total SNPs after splitting into European ancestry only: {eur_biallelic_snps}")
    print("\n")

#MAF
maf_snps_eur = 0

if args.verbose:
    print("Importing EUR MAF file...")

for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
    path = os.path.join(args.input_dir, "eur", "EUR_chr" + str(i) + "." + args.input_genotype_file + ".r2_0.4_biallelic_maf.0.01.bim")
    
    # Count lines that do not start with '#'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Skip lines starting with '#'
                maf_snps_eur += 1

if args.verbose:
    print(f"Total EUR SNPs after filtering by MAF: {maf_snps_eur}")
    print("\n")

#HWE
hwe_snps_eur = 0

if args.verbose:
    print("Importing EUR HWE file...")

for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
    path = os.path.join(args.input_dir, "eur", "EUR_chr" + str(i) + "." + args.input_genotype_file + ".r2_0.4_biallelic_maf.0.01_hwe.1e-8.bim")
    
    # Count lines that do not start with '#'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Skip lines starting with '#'
                hwe_snps_eur += 1

if args.verbose:
    print(f"Total EUR SNPs after HWE pruning: {maf_snps_eur}")
    print("\n")

##AFR
#Split into AFR ancestry
afr_biallelic_snps = 0

if args.verbose:
    print("Importing AFR biallelic file...")

for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
    path = os.path.join(args.input_dir, "afr", "AFR_chr" + str(i) + "." + args.input_genotype_file + ".r2_0.4_biallelic.bim")
    
    # Count lines that do not start with '#'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Skip lines starting with '#'
                afr_biallelic_snps += 1

if args.verbose:
    print(f"Total SNPs after splitting into African ancestry only: {afr_biallelic_snps}")
    print("\n")

#MAF
maf_snps_afr = 0

if args.verbose:
    print("Importing AFR MAF file...")

for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
    path = os.path.join(args.input_dir, "afr", "AFR_chr" + str(i) + "." + args.input_genotype_file + ".r2_0.4_biallelic_maf.0.01.bim")
    
    # Count lines that do not start with '#'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Skip lines starting with '#'
                maf_snps_afr += 1

if args.verbose:
    print(f"Total AFR SNPs after filtering by MAF: {maf_snps_afr}")
    print("\n")

#HWE
hwe_snps_afr = 0

if args.verbose:
    print("Importing AFR HWE file...")

for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
    path = os.path.join(args.input_dir, "afr", "AFR_chr" + str(i) + "." + args.input_genotype_file + ".r2_0.4_biallelic_maf.0.01_hwe.1e-6.bim")
    
    # Count lines that do not start with '#'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):  # Skip lines starting with '#'
                hwe_snps_afr += 1

if args.verbose:
    print(f"Total AFR SNPs after HWE pruning: {maf_snps_afr}")
    print("\n")
	
##Make output file
pops = ["EUR", "AFR"]
post_imputation_snps_array = [post_imputation_snps, post_imputation_snps]
r2_snps_array = [r2_snps, r2_snps]
biallelic_snps_array = [biallelic_snps, biallelic_snps]
ancestry_biallelic_snps_array = [eur_biallelic_snps, afr_biallelic_snps]
maf_snps = [maf_snps_eur, maf_snps_afr]
hwe_snps = [hwe_snps_eur, hwe_snps_afr]

out_df = pd.DataFrame(list(zip(pops, post_imputation_snps_array, r2_snps_array, biallelic_snps_array, ancestry_biallelic_snps_array, maf_snps, hwe_snps)),
               columns = ["pop", "post_imputation", "imputation_r2", "biallelic", "biallelic_ancestry_split", "maf", "hwe"])

out_df.to_csv(os.path.join(args.output), sep = "\t", header = True, index = False)