# Simulate effect sizes under LDSC model

require(devtools)
require(data.table)
require(simulateGP)
#devtools::install_github("explodecomputer/simulateGP")


freq <- fread("eur_hm3.frq")
bim <- fread("eur_hm3.bim")
bim2 <- bim[, c(2, 4)]
names(bim2) <- c("SNP", "BP")

for (i in 1:10)
{ 

#Outfile
outfile <- paste("/ReferenceSets/","Simulation",i, sep = "")

#Set seed
k <- i + 134
set.seed(k)

#Extract random subset of variants
freq$norm <- rnorm(nrow(freq), 0, 1)
freq2 <- freq[order(freq$norm), ]
freq3 <- freq2[1:1000]

#Generate effect sizes
effectsizes <- generate_gwas_params(freq3$MAF, 0.3)

#Tidy up
output1 <- cbind(freq3, effectsizes)
output2 <- merge(output1, bim2, by = "SNP")
output2$NCHROBS <- NULL
output2$maf <- NULL
output2$norm <- NULL

#Write to file
write.table(output2, outfile, quote = F, row.names = F)
}
