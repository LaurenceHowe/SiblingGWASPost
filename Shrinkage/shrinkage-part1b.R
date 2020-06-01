require(data.table)

# import argument
arguments <- commandArgs(trailingOnly = T)
phenotype <- arguments[1]

#Summary statistics
datafile <- paste(phenotype,"_imputed.txt", sep = "")
data <- fread(datafile)

#Clumped SNPs
clumped1file <- paste(phenotype,"_GWS.clumped", sep = "")
clumped1 <- fread(clumped1file)

clumped2file <- paste(phenotype,"_05.clumped", sep = "")
clumped2 <- fread(clumped2file)

#Restrict summary stats to clumped SNPs
subset1 <- data[which(data$SNP %in% clumped1$SNP), ]
subset2 <- data[which(data$SNP %in% clumped2$SNP), ]

#Subset and output
outfile1 <- paste(phenotype,"_GWS.txt", sep = "")
outfile2 <- paste(phenotype,"_05.txt", sep = "")
write.table(subset1, outfile1, quote = F, row.names = F)
write.table(subset2, outfile2, quote = F, row.names = F)
