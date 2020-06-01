require(data.table)

sds <- fread("temp")
names(sds) <- c("CHR", "POS", "ID", "AA", "DA", "DAF", "SDS", "SNP")

#Remove MHC and lactase: chr6: 25,892,529-33,436,144 and chr2: 134,608,646-138,608,646
sds2 <- sds[! which(sds$CHR == 6 & sds$POS >= 25892529 & sds$POS <= 33436144 | sds$CHR ==2 & sds$POS >= 134608646 & sds$POS <= 138608646), ]

#Normalise within each DAF bin of 1%

output <- NULL

for (i in 1:100) {
j <- (i-1)/100
k <- i/100

bin <- sds2[which(sds2$DAF > j & sds2$DAF <= k), ]
bin$SDS_STD <- scale(bin$SDS)

output <- rbind(output, bin)
}

write.table(output, "SDS_QC.txt", quote = F, row.names = F)
