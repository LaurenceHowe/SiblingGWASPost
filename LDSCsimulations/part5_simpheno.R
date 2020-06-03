require(data.table)

#Read in sibling genotypes in raw format and order
rawfile <- fread("Simulation1-tweaked.txt.raw")
rawfile2 <- rawfile[order(rawfile$FID), ]

#Read in allele frequencies
maf <- fread("Simulation1-tweaked.txt.frq")
maf$order <- 1:nrow(maf)

#Read in simulated parental genotypes
parent <- fread("parent-genotypes.txt", h = T)

#Read in SNP weightings for parental score
score <- fread("Simulation1-tweaked.txt")
names(score) <- c("rs", "SNP", "A1", "A2", "MAF", "Beta")

#Identify families with 2+ individuals
test <- rawfile2[duplicated(rawfile2$FID), ]
families <- unique(test$FID)


#Merge and harmonise alleles between score file/genotype file
merge <- merge(maf, score, by ="SNP")
merge$TweakedBeta[merge$A1.y == merge$A1.x] <- merge$Beta[merge$A1.y == merge$A1.x]
merge$TweakedBeta[merge$A1.y == merge$A2.x] <- (-1 * merge$Beta)[merge$A1.y == merge$A2.x]

merge2 <- merge[order(merge$order),]

#Run loop to generate parental genetic score for each family

output <- NULL

for (i in 1:nrow(merge2)) {
j <- merge2$TweakedBeta[i]
test <- as.vector(as.numeric(parent[i, 2:19525]))
test2 <- j * test
output <- cbind(output, test2)
}

final <- NULL

for (i in 1:19524) {
temp <- output[i,]
temp2 <- sum(temp)
print(temp2)
final <- rbind(final, temp2)
}

parentscores <- data.table(FID = families, SCORE = final) 

#Read in and merge with the individual sibling genetic scores
sibscores <- fread("Simulation1-tweaked.txt.profile")
newmerge <-  merge(sibscores, parentscores, by = "FID")

#Order by FID and generate scaled individual scores, parental scores and mean sibling scores.
newmerge2 <- newmerge[order(newmerge$FID), ]
newmerge2$SCORE_STD <- scale(newmerge2$SCORE)
newmerge2$PARENT_STD <- scale(newmerge2$SCORE.V1)
newmerge2$FAM_STD <- scale(ave(as.numeric(newmerge2$SCORE_STD), newmerge2$FID, FUN=mean))

#Estimate the correlations between the scores to demonstrate that they are correlated
cor(newmerge2$SCORE_STD, newmerge2$PARENT_STD)
#[1] 0.6998169
cor(newmerge2$SCORE_STD, newmerge2$FAM_STD)
#[1] 0.8625404
cor(newmerge2$PARENT_STD, newmerge2$FAM_STD)
#[1] 0.811344


#Simulate the phenotypes
data <- newmerge2

#Model A: Direct genetic effects explaining 30% of the phenotypic variance
data$ModelA <- scale(data$SCORE_STD + sqrt(7/3)*rnorm(nrow(data), 0, 1))

#Model B: Parental genetic effects explaining 30% of the phenotypic variance
data$ModelB <- scale(data$PARENT_STD + sqrt(7/3)*rnorm(nrow(data), 0, 1))

#Model C: DGE explaining 30% and PGE explaining an additional 30% of the phenotypic variance.
data$ModelC <- scale(scale(data$SCORE_STD + data$PARENT_STD) + sqrt(4/6)*rnorm(nrow(data), 0, 1))

