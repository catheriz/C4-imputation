module load bcftools
bcftools annotate -x ^INFO/DR2,^INFO/AF,^FORMAT/DS DMPM_C4_imputation_hg19.vcf.gz -Oz -o DMPM_C4_imputation_hg19_DS.vcf.gz

(echo -e "ALT\tDR2\tAF\t$(bcftools query -l DMPM_C4_imputation_hg19_DS.vcf.gz | tr '\n' '\t' | sed 's/\t$//')"; \
bcftools view -i 'ID=="C4"' DMPM_C4_imputation_hg19_DS.vcf.gz | bcftools query -f '%ALT\t%INFO/DR2\t%INFO/AF[\t%DS]\n') > dosages.txt


