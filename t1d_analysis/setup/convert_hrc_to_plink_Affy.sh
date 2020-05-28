#convert_hrc_to_plink.sh

#converting the VCF files from Anubha to PLINK files that I am more familiar wit hto use in SnpStats
#Also going to filter these genotypes to keep only SNPs imputed with high certainty (r2>0.3) and maf>0.01


for i in {1..22}
do
#generate the info list for each SNP
zcat /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.vcf.gz | grep -v '#' | awk -F  '\t' '{print $3, $8}' | tee /well/todd/users/jinshaw/t1d_risk/processed/look 

#run the R script to keep variant names with an imputation r2 of >0.3
Rscript filterr2.R $i


/apps/well/vcftools/16102015/bin/vcftools --gzvcf /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.vcf.gz --min-alleles 2 --max-alleles 2 --recode --maf 0.01 \
--exclude /well/todd/users/jinshaw/t1d_risk/processed/exclude_imp_qual  --out /well/todd/users/jinshaw/t1d_risk/processed/s1.vcf
/apps/well/plink/1.90b3/plink --vcf /well/todd/users/jinshaw/t1d_risk/processed/s1.vcf.recode.vcf --make-bed --out /well/todd/users/jinshaw/t1d_risk/processed/affy_chr$i
done

