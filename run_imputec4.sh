vcf="..."
out="..."
build=38 # build=37
declare -A reg=( ["37"]="6:24894177-33890574" ["38"]="chr6:24893949-33922797" )

bcftools view --no-version "$vcf" -r ${reg[$build]} | \
  java -Xmx8g -jar /mount/ictr1/Users/catheriz/imputec4/beagle.25Nov19.28d.jar gt=/dev/stdin gp=true \
  ref=/mount/ictr1/Users/catheriz/imputec4/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf.gz out="$out" \
  map=<(bcftools query -f "%CHROM\t%POS\n" /mount/ictr1/Users/catheriz/imputec4/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf.gz | \
  awk '{print $1"\t.\t"$2/1e7"\t"$2}')
