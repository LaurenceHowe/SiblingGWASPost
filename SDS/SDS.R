require(data.table)

sds <- fread("SDS_QC.txt")

arguments <- commandArgs(trailingOnly = T)
phenotype <- arguments[1]
filename <- paste(",phenotype,"/",phenotype,"WFTotal.txt", sep = "")

data <- fread(filename)
data2 <- data[which(data$SE_BETA_WF > 0 & data$SE_BETA_TOTAL > 0), ]
data3 <- data2[which(data2$N_REG_TOTAL >= max(data2$N_REG_TOTAL)/2), ]

merge <- merge(data3, sds, by = c("SNP", "CHR"))

final <- NULL
for (i in 1:2) {

if (i == 1) {
merge$BETA <- merge$BETA_TOTAL
merge$PV <- merge$P_BETA_TOTAL
}

if (i == 2) {
merge$BETA <- merge$BETA_WF
merge$PV <- merge$P_BETA_WF
}

merge$SDS_tweak[merge$A1 == merge$DA & merge$BETA >0] <- merge$SDS_STD[merge$A1 == merge$DA & merge$BETA >0]
merge$SDS_tweak[merge$A1 == merge$DA & merge$BETA <0] <- (-1)*merge$SDS_STD[merge$A1 == merge$DA & merge$BETA <0]
merge$SDS_tweak[merge$A2 == merge$DA & merge$BETA <0] <- merge$SDS_STD[merge$A2 == merge$DA & merge$BETA <0]
merge$SDS_tweak[merge$A2 == merge$DA & merge$BETA >0] <- (-1)*merge$SDS_STD[merge$A2 == merge$DA & merge$BETA >0]

merge <- merge[!which(is.na(merge$SDS_tweak)), ]

#Spearmans Rank on whole sample
srtest <- cor.test(merge$PV, merge$SDS_tweak,  method = "spearman")
print("Whole sample")
print(srtest)

#Jackknife blocks
merge2 <- merge[order(merge$CHR, merge$BP), ]
merge2$Order <- 1:nrow(merge2)

#Block size to exclude
nb <- 100
k <- round(nrow(merge2)/nb)

#Loop

output <- NULL

for (i in 1:nb) {
l <- k * (i-1)
m <- k * i

if (i == nb) {
m <- nrow(merge2)
}

block <- merge2[! which(merge2$Order > l & merge2$Order <= m), ]
srtest <- cor.test(block$PV, block$SDS_tweak,  method = "spearman")

temp <- as.numeric(srtest$estimate)
output <- rbind(output, temp)
}

output2 <- as.data.frame(output)
output2$V2 <- (mean(output2$V1) - output2$V1) ^2
print(mean(output2$V1))
print(sqrt(sum(output2$V2)* ((nb-1)/nb)))

set <- rbind(c(phenotype, mean(output2$V1), sqrt(sum(output2$V2)* ((nb-1)/nb))))
final <- rbind(final, set)
}

outname <- paste(phenotype,".SDS.txt", sep = "")
write.table(final, outname, quote = F, row.names = F, col.names = F)
