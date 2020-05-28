#rename_variants_Affy.sh

#ensuring each snp is names chrom:pos:alleleA:alleleB as opposed to chrom:pos

for i in {1..22}
do
tabix -H /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.dose.vcf.gz > /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.vcf
zgrep -v '#' /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.dose.vcf.gz  | awk '$3=$3":"$4":"$5' OFS="\t" >> /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.vcf
bgzip /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.vcf
tabix /well/todd/users/jinshaw/t1d_risk/T1D.Affy.chr$i.vcf.gz
done
