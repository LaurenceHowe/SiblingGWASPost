# Generate allele frequencies of HM3 SNPs 
plink \
--bfile eur \
--extract hm3.snplist \
--make-bed \
--out eur_hm3


plink \
--bfile eur_hm3 \
--freq \
--out eur_hm3

