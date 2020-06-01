phen=${1}

plink \
--bfile ../eur \
--clump ${phen}_imputed.txt \
--clump-snp-field SNP \
--clump-field P_BOLT_LMM_INF \
--clump-p1 0.00000005 \
--clump-p2 0.00000005 \
--clump-r2 0.001 \
--clump-kb 10000 \
--out ${phen}_GWS

plink \
--bfile ../eur \
--clump ${phen}_imputed.txt \
--clump-snp-field SNP \
--clump-field P_BOLT_LMM_INF \
--clump-p1 0.00001 \
--clump-p2 0.00001 \
--clump-r2 0.001 \
--clump-kb 10000 \
--out ${phen}_05

module add languages/r/3.6.0

Rscript \
shrinkage-part1b.R \
${phen}
