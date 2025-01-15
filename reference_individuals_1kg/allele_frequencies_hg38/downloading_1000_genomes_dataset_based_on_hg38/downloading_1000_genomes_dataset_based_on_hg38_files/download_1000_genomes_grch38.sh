#Download each of the individual (by chr) 1000 Genomes genotype VCF files and their corresponding indices:
sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr1_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr1_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr1_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr2_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr2_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr2_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr3_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr3_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr3_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr4_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr4_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr4_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr5_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr5_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr5_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr6_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr6_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr6_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr7_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr7_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr7_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr8_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr8_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr8_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr9_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr9_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr9_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr10_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr10_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr10_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr11_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr11_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr11_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr12_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr12_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr12_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr13_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr13_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr13_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr14_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr14_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr14_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr15_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr15_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr15_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr16_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr16_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr16_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr17_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr17_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr17_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr18_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr18_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr18_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr19_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr19_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr19_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr20_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr20_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr20_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr21_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr21_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr21_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chr22_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chr22_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chrX_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrX_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chrX_GRCh38.genotypes.20170504.vcf.gz"

sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=8G \
--time=0-02:00:00 \
--job-name=chrY_download \
--wrap="wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrY_GRCh38.genotypes.20170504.vcf.gz
tabix -f -p vcf ALL.chrY_GRCh38.genotypes.20170504.vcf.gz"



