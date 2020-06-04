##Script to generate parental genotypes based on sibling genotypes
##Assumes laws of Mendelian inheritance, random mating and HWE

require(data.table)

rawfile <- fread("Simulation1-tweaked.txt.raw")
rawfile2 <- rawfile[order(rawfile$FID), ]
maf <- fread("Simulation1-tweaked.txt.frq")

#Define variables
maf$off00P2 <- NA 
maf$off00P1 <- NA
maf$off00P0 <- NA

maf$off01P2 <- NA
maf$off01P1 <- NA

maf$off02P2 <- NA

maf$off11P3 <- NA
maf$off11P2 <- NA
maf$off11P1 <- NA

maf$off12P3 <- NA
maf$off12P2 <- NA

maf$off22P4 <- NA
maf$off22P3 <- NA
maf$off22P2 <- NA

##


for (i in 1:nrow(maf)) 
{
freq <- maf$MAF[i]

freq2 <- freq ^2
freq1 <- 2 * freq * (1-freq)
freq0 <- (1-freq) ^2

freq22 <- freq2 ^2
freq21 <- 2 * freq2 * freq1
freq20 <- 2 * freq2 * freq0
freq11 <- freq1 ^2 
freq10 <- 2 * freq1 * freq0
freq00 <- freq0 ^2

freq22 + freq21 + freq20 + freq11 + freq10 + freq00


total00 <- (freq11 / 16) + (freq10 / 4) + freq00
maf$off00P2[i] <- (freq11 / 16) / total00
maf$off00P1[i] <- (freq10 / 4) / total00
maf$off00P0[i] <- freq00 / total00

total01 <- (freq11 / 4) + (freq10 / 2)
maf$off01P2[i] <- (freq11 / 4) / total01
maf$off01P1[i] <- (freq10 / 2) / total01

total02 <- (freq11 / 8)
maf$off02P2[i] <- 1

total11 <- (freq11 / 4) + freq20 + (freq10 / 4) + (freq21 / 4)
maf$off11P3[i] <- (freq21 / 4) / total11
maf$off11P2[i] <- ((freq11 / 4) + freq20) / total11
maf$off11P1[i] <- (freq10 / 4) / total11

total12 <- (freq11 / 4) + (freq21 / 2)
maf$off12P3[i] <- (freq21 / 2) / total12
maf$off12P2[i] <- (freq11 / 4) / total12

total22 <- (freq11 / 16) + (freq21 / 4) + freq22
maf$off22P4[i] <- freq22 / total22
maf$off22P3[i] <- (freq21 / 4) / total22
maf$off22P2[i] <- (freq11 / 16) / total22
}

test <- rawfile2[duplicated(rawfile2$FID), ]
families <- unique(test$FID)
output <- data.frame(maf$SNP)

for (j in 1:length(families)) 

for (j in 1:100) 
{
restrict <- families[j]

temp <- rawfile2[which(rawfile2$FID == restrict), ]
temp2 <- temp[1:2, ]


parent <- data.table(SNP = maf$SNP, Parent = NA)
print(j)

for (i in 1:nrow(maf)) 
{

snpvector <- NULL
k <- i + 6
snpvector <- c(as.numeric(temp2[1,..k]), as.numeric(temp2[2,..k]))

random <- runif(1, 0, 1)

#Parental probabilities

if (snpvector[1] == 0 & snpvector[2] == 0) 
{
thresh1 <- (1- maf$off00P2[i])
thresh2 <- (thresh1 - maf$off00P1[i])
parent$Parent[i] <- 0
parent$Parent[i][random > thresh2] <- 1
parent$Parent[i][random > thresh1] <- 2
}

if (snpvector[1] == 1 & snpvector[2] == 1) 
{
thresh1 <- (1- maf$off11P3[i])
thresh2 <- (thresh1 - maf$off11P2[i])
parent$Parent[i] <- 1
parent$Parent[i][random > thresh2] <- 2
parent$Parent[i][random > thresh1] <- 3
}

if (snpvector[1] == 2 & snpvector[2] == 2) 
{
thresh1 <- (1- maf$off22P4[i])
thresh2 <- (thresh1 - maf$off22P3[i])
parent$Parent[i] <- 2
parent$Parent[i][random > thresh2] <- 3
parent$Parent[i][random > thresh1] <- 4
}

if (snpvector[1] == 0 & snpvector[2] == 2 | snpvector[1] == 2 & snpvector[2] == 0) 
{
parent$Parent[i] <- 2
}


if (snpvector[1] == 0 & snpvector[2] == 1 | snpvector[1] == 1 & snpvector[2] == 0) 
{
thresh1 <- (1 - maf$off01P2[i])
parent$Parent[i] <- 1
parent$Parent[i][random > thresh1] <- 2
}


if (snpvector[1] == 1 & snpvector[2] == 2 | snpvector[1] == 2 & snpvector[2] == 1) 
{
thresh1 <- (1 - maf$off12P3[i])
parent$Parent[i] <- 2
parent$Parent[i][random > thresh1] <- 3
}
}
parent$SNP <- NULL
output <- cbind(output, parent)
}

names(output) <- c("SNP", families)
write.table(output, "parent-genotypes.txt", quote = F, row.names = F)
