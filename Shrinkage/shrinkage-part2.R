require(data.table)
require(metafor)
require(msm)

arguments <- commandArgs(trailingOnly = T)
sumstatsfile <- arguments[1]
phenotype <- arguments[2]
set <- as.numeric(arguments[3])

sumstats <- fread(sumstatsfile)
print(sumstatsfile)

if (set == 1) {
filename <- paste(,phenotype,"_GWS.txt", sep = "")
snps <- fread(filename, h = T)
print("GWS")
}

if (set == 2) {
filename <- paste(,phenotype,"_05.txt", sep = "")
snps <- fread(filename, h = T)
print("10-5")
}


merge <- merge(sumstats, snps,by = c("CHR", "BP"))
print(nrow(merge))

merge$Match[merge$A1 == merge$ALLELE1] <- 1
merge$Match[merge$A1 == merge$ALLELE0] <- 0

merge <- merge[! is.na(merge$Match), ]

merge$BETA_TOTAL_H[merge$A1 == merge$ALLELE1] <- merge$BETA_TOTAL[merge$A1 == merge$ALLELE1]
merge$BETA_TOTAL_H[merge$A1 == merge$ALLELE0] <- (-1)*merge$BETA_TOTAL[merge$A1 == merge$ALLELE0]
merge$TOTAL_TOP <- merge$BETA * merge$BETA_TOTAL_H * (1/merge$SE_BETA_TOTAL^2)
merge$TOTAL_BOT <- merge$BETA^2 * (1/merge$SE_BETA_TOTAL^2)


IVWBeta_TOTAL <- sum(merge$TOTAL_TOP) / sum(merge$TOTAL_BOT)
IVWSE_TOTAL <- 1 / sum(merge$TOTAL_BOT)
print(c(IVWBeta_TOTAL, IVWSE_TOTAL))

merge$BETA_WF_H[merge$A1 == merge$ALLELE1] <- merge$BETA_WF[merge$A1 == merge$ALLELE1]
merge$BETA_WF_H[merge$A1 == merge$ALLELE0] <- (-1)*merge$BETA_WF[merge$A1 == merge$ALLELE0]
merge$WF_TOP <- merge$BETA * merge$BETA_WF_H * (1/merge$SE_BETA_WF^2)
merge$WF_BOT <- merge$BETA^2 * (1/merge$SE_BETA_WF^2)

IVWBeta_WF <- sum(merge$WF_TOP) / sum(merge$WF_BOT)
IVWSE_WF <- sqrt(1 / sum(merge$WF_BOT))
print(c(IVWBeta_WF, IVWSE_WF))

Ratio <- IVWBeta_WF / IVWBeta_TOTAL

cov_matrix <- matrix(c(IVWSE_WF^2, 0, 0, IVWSE_TOTAL^2), nrow=2, ncol=2)
RatioSE <- deltamethod(g=~(x1)/x2, mean=c(IVWBeta_WF, IVWBeta_TOTAL), cov=cov_matrix)

#Jackknife blocks for each SNP

k <- nrow(merge)
merge$Order <- 1:k

#Loop

output <- NULL

for (i in 1:k) {

block <- merge[! which(merge$Order == i), ]

block$BETA_TOTAL_H[block$A1 == block$ALLELE1] <- block$BETA_TOTAL[block$A1 == block$ALLELE1]
block$BETA_TOTAL_H[block$A1 == block$ALLELE0] <- (-1)*block$BETA_TOTAL[block$A1 == block$ALLELE0]
block$TOTAL_TOP <- block$BETA * block$BETA_TOTAL_H * (1/block$SE_BETA_TOTAL^2)
block$TOTAL_BOT <- block$BETA^2 * (1/block$SE_BETA_TOTAL^2)


temp1 <- sum(block$TOTAL_TOP) / sum(block$TOTAL_BOT)


block$BETA_WF_H[block$A1 == block$ALLELE1] <- block$BETA_WF[block$A1 == block$ALLELE1]
block$BETA_WF_H[block$A1 == block$ALLELE0] <- (-1)*block$BETA_WF[block$A1 == block$ALLELE0]
block$WF_TOP <- block$BETA * block$BETA_WF_H * (1/block$SE_BETA_WF^2)
block$WF_BOT <- block$BETA^2 * (1/block$SE_BETA_WF^2)

temp2 <- sum(block$WF_TOP) / sum(block$WF_BOT)

temp3 <- cbind(temp1, temp2)

output <- rbind(output, temp3)
}


output2 <- as.data.frame(output)
names(output2) <- c("Total", "WF")
output2$Ratio <- output2$WF / output2$Total

output2$SS <- (mean(output2$Ratio) - output2$Ratio) ^2
print(mean(output2$Ratio))
BootSE <- sqrt(sum(output2$SS)* ((k-1)/k))


#Original SE from delta method
print(c(1-Ratio, 1-(Ratio + 1.96*RatioSE), 1-(Ratio - 1.96*RatioSE)))

#Bootstrapped SEs
print(c(1-Ratio, 1-(Ratio + 1.96*BootSE), 1-(Ratio - 1.96*BootSE)))


####IVW/Heterogeneity code (adapted from Gibran Hemani code)

merge$W <- (merge$BETA_TOTAL ^2) / merge$SE_BETA_WF^2
merge$Wj <- sqrt(merge$W)

merge$SNPRatios <- merge$BETA_WF / merge$BETA_TOTAL
merge$BetaWj <- merge$SNPRatios * merge$Wj

mod <- lm(BetaWj ~ -1 + Wj, data = merge)$coef[1]
merge$Wj <- ((merge$SE_BETA_WF^2 + (mod^2*merge$SE_BETA_TOTAL^2)) / merge$BETA_TOTAL^2)^-1 %>% sqrt
merge$BetaWj <- merge$SNPRatios * merge$Wj

mod2 <- summary(lm(BetaWj ~ -1 + Wj, data = merge))
ivw <- coefficients(mod2)[1,1]
ivw_se <- coefficients(mod2)[1,2]

merge$Qj <- merge$Wj^2 * (merge$SNPRatios - ivw)^2
Q <- sum(merge$Qj)
Qpval <- pchisq(Q, nrow(merge)-1, lower.tail=FALSE)
