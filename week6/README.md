1.1
/usr/local/bin/plink --vcf genotypes.vcf --pca -out gtPCA

2.1
/usr/local/bin/plink --freq --vcf genotypes.vcf

3.1
/usr/local/bin/plink --vcf genotypes.vcf --linear --pheno GS451_IC50.txt --covar gtPCA.eigenvec --allow-no-sex --out GS451_IC50_GWAS

/usr/local/bin/plink --vcf genotypes.vcf --linear --pheno CB1908_IC50.txt --covar gtPCA.eigenvec --allow-no-sex --out CB1908_IC50_GWAS

3.4
The rs17113501 SNP is located in the ACTR1A gene at chromosome 10 base pair 104230536. This gene encodes a subunit of dynactin, which transports molecules in the cell. The site of this SNP is conserved in mammal species, and is an important gene in a convserved cellular pathway. This may be why I did not see any variants with TT, only homozygotes for A or heterozygotes. 
