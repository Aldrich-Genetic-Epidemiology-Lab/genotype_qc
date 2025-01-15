sbatch \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=1 \
--mem-per-cpu=128G \
--time=1-00:00:00 \
--job-name=rm_multiallelic \
--wrap="bcftools view \
--max-alleles 2 \
/home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/ALL.GRCh38.genotypes.20170504.sorted.vcf.gz \
> /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.vcf

bgzip /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.vcf

tabix -f -p vcf /home/bettimj/gamazon_rotation/1000_genomes/compile_all_grch38_vcf/ALL.GRCh38.genotypes.20170504.sorted.no.multiallelic.vcf.gz"