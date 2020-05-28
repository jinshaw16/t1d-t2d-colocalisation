#convert_hrc_to_plink_illu.sh

#converting the VCF files from Anubha to PLINK files that I am more familiar wit hto use in SnpStats
#Also going to filter these genotypes to keep only SNPs imputed with high certainty (r2>0.3) and maf>0.1


for i in {1..22}
do
#generate the info list for each SNP
zcat /well/todd/users/jinshaw/t1d_risk/T1D.Illumina.chr$i.vcf.gz | grep -v '#' | awk -F  '\t' '{print $3, $8}' | tee /well/todd/users/jinshaw/t1d_risk/processed/look1 

#run the R script to keep variant names with an imputation r2 of >0.3
Rscript filterr2_illu.R $i


/apps/well/vcftools/16102015/bin/vcftools --gzvcf /well/todd/users/jinshaw/t1d_risk/T1D.Illumina.chr$i.vcf.gz --recode --maf 0.01 \
--exclude /well/todd/users/jinshaw/t1d_risk/processed/exclude_imp_qual_illu  --out /well/todd/users/jinshaw/t1d_risk/processed/t1.vcf
/apps/well/plink/1.90b3/plink --vcf /well/todd/users/jinshaw/t1d_risk/processed/t1.vcf.recode.vcf --const-fid 0 --make-bed --out /well/todd/users/jinshaw/t1d_risk/processed/illu_chr$i
done
