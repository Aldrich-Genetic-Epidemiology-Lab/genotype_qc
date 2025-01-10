import argparse
import os
from dask import dataframe as dd
import dask.array as da
from dask.diagnostics import ProgressBar
ProgressBar().register()
import pandas as pd

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-d", "--input_dir", type = str, required = True, help = "the input file directory (required)")
parser.add_argument("-ii", "--input_imputation_file", type = str, required = True, help = "the name of the input post-imputation SNP list with info scores (required)")
parser.add_argument("-is", "--input_info_score_file", type = str, required = True, help = "the name of the input SNP list after filtering by info score (required)")
parser.add_argument("-gi", "--input_genotype_file", type = str, required = True, help = "the root name of the input genotype files (required)")
parser.add_argument("-o", "--output", type = str, required = True, help = "the path of the output file")
parser.add_argument("-v", "--verbose", required = False, help = "return logging as terminal output", action = "store_true")

args = parser.parse_args()

##Combined dataset
#Post-imputation
post_imputation_snps = 0

if args.verbose:
	print("Importing imputed SNPs...")
	print("\n")
path = os.file.path(args.input_dir, args.input_imputation_file)
file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
df = file.compute()
post_imputation_snps = len(df.index)

#INFO score thresholding
info_snps = 0

if args.verbose:
	print("Importing remaining SNPs after INFO score filtering...")
	print("\n")
path = os.file.path(args.input_dir, args.input_info_score_file)
file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
df = file.compute()
info_snps = len(df.index)

#Biallelic SNPs
biallelic_snps = 0

if args.verbose:
    print("Importing biallelic file...")
    print("\n")
for i in range(1, 23):
	if args.verbose:
		print("Processing chr" + str(i) + "...")
		print("\n")
	path = os.file.path(args.input_dir, "chr" + i + ".info_0.4_" + args.input_genotype_file + "_biallelic.bim")
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	biallelic_snps_chr = len(df.index)
	biallelic_snps += biallelic_snps_chr

##EUR
#Split into EUR ancestry
eur_biallelic_snps = 0

if args.verbose:
    print("Importing EUR biallelic file...")
    print("\n")
for i in range(1, 23):
	if args.verbose:
		print("Processing chr" + str(i) + "...")
		print("\n")
	path = os.file.path(args.input_dir, "EUR_chr" + i + ".info_0.4_" + args.input_genotype_file + "_biallelic.bim")
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	eur_biallelic_snps_chr = len(df.index)
	eur_biallelic_snps += eur_biallelic_snps_chr

#MAF
maf_snps_eur = 0

if args.verbose:
    print("Importing MAF file...")
    print("\n")
for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
        print("\n")
    path = os.file.path(args.input_dir, "EUR_chr" + i + ".info_0.4_" + args.input_genotype_file + "_biallelic_maf.0.01.bim")
    file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
    df = file.compute()
    maf_snps_chr_eur = len(df.index)
    maf_snps_eur += maf_snps_chr_eur

#HWE
hwe_snps_eur = 0

if args.verbose:
    print("Importing HWE file...")
    print("\n")
for i in range(1, 23):
	if args.verbose:
	    print("Processing chr" + str(i) + "...")
		print("\n")
    path = os.file.path(args.input_dir, "EUR_chr" + i + ".info_0.4_" + args.input_genotype_file + "_biallelic_maf.0.01_hwe.1e-8.bim")
    file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
    df = file.compute()
    hwe_snps_chr = len(df.index)
    hwe_snps_eur += hwe_snps_chr

##AFR
#Split into AFR ancestry
afr_biallelic_snps = 0

if args.verbose:
    print("Importing AFR biallelic file...")
    print("\n")
for i in range(1, 23):
	if args.verbose:
		print("Processing chr" + str(i) + "...")
		print("\n")
	path = os.file.path(args.input_dir, "AFR_chr" + i + ".info_0.4_" + args.input_genotype_file + "_biallelic.bim")
	file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
	df = file.compute()
	afr_biallelic_snps_chr = len(df.index)
	afr_biallelic_snps += afr_biallelic_snps_chr

#MAF
maf_snps_afr = 0

if args.verbose:
    print("Importing MAF file...")
    print("\n")
for i in range(1, 23):
    if args.verbose:
        print("Processing chr" + str(i) + "...")
        print("\n")
    path = os.file.path(args.input_dir, "AFR_chr" + i + ".info_0.4_" + args.input_genotype_file + "_biallelic_maf.0.01.bim")
    file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
    df = file.compute()
    maf_snps_chr_afr = len(df.index)
    maf_snps_afr += maf_snps_chr_afr

#HWE
hwe_snps_afr = 0

if args.verbose:
    print("Importing HWE file...")
    print("\n")
for i in range(1, 23):
	if args.verbose:
	    print("Processing chr" + str(i) + "...")
		print("\n")
    path = os.file.path(args.input_dir, "AFR_chr" + i + ".info_0.4_" + args.input_genotype_file + "_biallelic_maf.0.01_hwe.1e-6.bim")
    file = dd.read_csv(path, sep = "\t", low_memory = False, dtype = {'0': 'float64', '1': 'object'})
    df = file.compute()
    hwe_snps_chr = len(df.index)
    hwe_snps_afr += hwe_snps_chr
	
##Make output file
pops = ["EUR", "AFR"]
post_imputation_snps_array = [post_imputation_snps, post_imputation_snps]
info_snps_array = [info_snps, info_snps]
biallelic_snps_array = [biallelic_snps, biallelic_snps]
ancestry_biallelic_snps_array = [eur_biallelic_snps, afr_biallelic_snps]
maf_snps = [maf_snps_eur, maf_snps_afr]
hwe_snps = [hwe_snps_eur, hwe_snps_afr]

out_df = pd.DataFrame(list(zip(pops, post_imputation_snps_array, info_snps_array, biallelic_snps_array, maf_snps, hwe_snps)),
               columns = ["pop", "post_imputation", "info_score", "biallelic", "biallelic_ancestry_split", "maf", "hwe"])

out_df.to_csv(os.file.path(args.output), sep = "\t", header = True, index = False)