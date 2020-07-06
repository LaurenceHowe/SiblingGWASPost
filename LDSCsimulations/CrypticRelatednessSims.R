#Install packages
#devtools::install_github("explodecomputer/simulateGP")
require(simulateGP)
require(data.table)
require(devtools)
require(sandwich)
require(lmtest)

#Set seed
set.seed(57)

#Set allele frequency
freq <- list()
freq <- 0.33

#Set sample size
n <- 10000

#Set output parameters to null
totalout1 <- NULL
wfout1 <- NULL
totalout2 <- NULL
wfout2 <- NULL

noclust1 <- NULL
noclust2 <- NULL

#Loop
for (i in 1:10000)
{
# Simulate common allele for men and women in Generation 1 (G1)
men1 <- make_geno(nid = n, nsnp = 1, af = freq)
women1 <- make_geno(nid = n, nsnp = 1, af = freq)

#Randomly pair up men/women in G1
data <- data.frame(cbind(men1, women1))

#First generation meiosis
#Male transmission to Child 1
data$MT1[data$X1 == 2] <- 1  
data$MT1[data$X1 == 0] <- 0

#Chance vector if male is heterozygous, 50% of transmitting each
data$MC1 <- rbinom(n, 1, 0.5) 
data$MT1[data$X1 == 1] <- data$MC1[data$X1 == 1]

#Male transmission to Child 2
data$MT2[data$X1 == 2] <- 1  
data$MT2[data$X1 == 0] <- 0

#Chance vector if male is heterozygous, 50% of transmitting each
data$MC2 <- rbinom(n, 1, 0.5) 
data$MT2[data$X1 == 1] <- data$MC2[data$X1 == 1]

#Female transmission to Child 1
data$WT1[data$X2 == 2] <- 1  
data$WT1[data$X2 == 0] <- 0
data$WC1 <- rbinom(n, 1, 0.5) 
data$WT1[data$X2 == 1] <- data$WC1[data$X2 == 1]

#Female transmission to Child 2
data$WT2[data$X2 == 2] <- 1  
data$WT2[data$X2 == 0] <- 0
data$WC2 <- rbinom(n, 1, 0.5) 
data$WT2[data$X2 == 1] <- data$WC2[data$X2 == 1]

#Genotypes for Generation 2 Males and Females
data$G2M <- data$MT1 + data$WT1
data$G2F <- data$MT2 + data$WT2

gen2 <- data.table(M1 = data$G2M, F1 = data$G2F)
gen2$FID <- 1:nrow(gen2)

#Simulate common environment term for Generation 2 families (this term is the same for siblings in G2)
gen2$CE <- rnorm(n, 0, 1)


#Extract relevant phenotypes and sort out random mating for Generation 2
gen2male <- gen2[, c(1, 3, 4)]
names(gen2male) <- c("M1", "MFID", "MCE")
gen2male$Order <- rnorm(nrow(gen2male), 0, 1)
temp <- gen2male[order(gen2male$Order), ]

gen2female <- gen2[, c(2, 3, 4)]
names(gen2female) <- c("F1", "FFID", "FCE")

gen2pairs <- cbind(temp, gen2female)
gen2pairs$Order <- NULL

#Second generation meiosis to create Generation 3
#Male transmission to Offspring 1
gen2pairs$MT1[gen2pairs$M1 == 2] <- 1  
gen2pairs$MT1[gen2pairs$M1 == 0] <- 0
gen2pairs$MC1 <- rbinom(n, 1, 0.5) 
gen2pairs$MT1[gen2pairs$M1 == 1] <- gen2pairs$MC1[gen2pairs$M1 == 1]

#Male transmission to Offspring 2
gen2pairs$MT2[gen2pairs$M1 == 2] <- 1  
gen2pairs$MT2[gen2pairs$M1 == 0] <- 0
gen2pairs$MC2 <- rbinom(n, 1, 0.5) 
gen2pairs$MT2[gen2pairs$M1 == 1] <- gen2pairs$MC2[gen2pairs$M1 == 1]

#Female transmission to Offspring 1
gen2pairs$WT1[gen2pairs$F1 == 2] <- 1  
gen2pairs$WT1[gen2pairs$F1 == 0] <- 0
gen2pairs$WC1 <- rbinom(n, 1, 0.5) 
gen2pairs$WT1[gen2pairs$F1 == 1] <- gen2pairs$WC1[gen2pairs$F1 == 1]

#Female transmission to Offspring 2
gen2pairs$WT2[gen2pairs$F1 == 2] <- 1  
gen2pairs$WT2[gen2pairs$F1 == 0] <- 0
gen2pairs$WC2 <- rbinom(n, 1, 0.5) 
gen2pairs$WT2[gen2pairs$F1 == 1] <- gen2pairs$WC2[gen2pairs$F1 == 1]

#Genotypes for Generation 3
gen2pairs$G3M <- gen2pairs$MT1 + gen2pairs$WT1
gen2pairs$G3F <- gen2pairs$MT2 + gen2pairs$WT2

#Create a parental effect using Generation 2 genotypes
gen2pairs$Parent <- gen2pairs$M1 + gen2pairs$F1

#Extract relevant phenotypes for Generation 3
gen3 <- data.table(M1 = gen2pairs$G3M, F1 = gen2pairs$G3F, FID1 = gen2pairs$MFID, FID2 = gen2pairs$FFID, 
                   CE1 = gen2pairs$MCE, CE2 = gen2pairs$FCE, Parent = gen2pairs$Parent)

#Generate CE term for Generation 3 as the average of their parents CE terms. 
#Will be correlated between cousins ~ 0.5 because sibling parents had same CE terms.
gen3$CE <- scale(0.5 * (gen3$CE1 + gen3$CE2))

#New FID for siblings
gen3$NewFID <- 1:n

#Simulate phenotype influenced by CE for Offspring 1 and Offspring 2
gen3$Phen1 <- scale(gen3$CE + sqrt(7/3)*rnorm(n, 0, 1))
gen3$Phen2 <- scale(gen3$CE + sqrt(7/3)*rnorm(n, 0, 1))


#Simulate phenotype influenced by parental genotype for Offspring 1 and Offspring 2
gen3$Phen3 <- scale(scale(gen3$Parent) + sqrt(10)*rnorm(n, 0, 1))
gen3$Phen4 <- scale(scale(gen3$Parent) + sqrt(10)*rnorm(n, 0, 1))


#Convert format from wide to long
gen3A <- gen3[, c(9, 10, 12, 1)]
names(gen3A) <- c("FID", "PhenCE", "PhenP", "G")
gen3B <- gen3[, c(9, 11, 13, 2)]
names(gen3B) <- c("FID", "PhenCE", "PhenP", "G")

gen3C <- rbind(gen3A, gen3B)

#Generate Family Mean genotype + Centred Genotypes
gen3C$G_MEAN <- ave(gen3C$G, gen3C$FID, FUN=mean)
gen3C$G_CENTRED <- gen3C$G - gen3C$G_MEAN


#Fit regression models
#Total model with CE phenotype
lm1 <- lm(PhenCE ~ G, data = gen3C)
#WF model with CE phenotype
lm2 <- lm(PhenCE ~ G_MEAN + G_CENTRED, data = gen3C)

#Total model with dynastic effect phenotype
lm3 <- lm(PhenP ~ G, data = gen3C)
#WF model with dynastic effect phenotype
lm4 <- lm(PhenP ~ G_MEAN + G_CENTRED, data = gen3C)

#Extract non-clustered P-values
noclust1 <- rbind(noclust1, summary(lm1)$coefficients[2,4])
noclust2 <- rbind(noclust2, summary(lm2)$coefficients[3,4])

#Clustering of standard errors across sibships
vcv_matrix1 <- vcovCL(lm1, cluster = gen3C$FID)
vcv_matrix2 <- vcovCL(lm2, cluster = gen3C$FID)
vcv_matrix3 <- vcovCL(lm3, cluster = gen3C$FID)
vcv_matrix4 <- vcovCL(lm4, cluster = gen3C$FID)

test_matrix1 <- coeftest(lm1, vcov.=vcv_matrix1)
test_matrix2 <- coeftest(lm2, vcov.=vcv_matrix2)
test_matrix3 <- coeftest(lm3, vcov.=vcv_matrix3)
test_matrix4 <- coeftest(lm4, vcov.=vcv_matrix4)

temp1 <- test_matrix1[2,4] 
temp2 <- test_matrix2[3,4]
temp3 <- test_matrix3[2,4] 
temp4 <- test_matrix4[3,4]

#Extract clustered P-values
totalout1 <- rbind(totalout1, temp1)
totalout2 <- rbind(totalout2, temp3)

wfout1 <- rbind(wfout1, temp2)
wfout2 <- rbind(wfout2, temp4)
}

#Total model, CE: no clustering
summary(noclust1)
table(noclust1 < 0.05)

#Total model, CE: clustering on sibs
summary(totalout1)
table(totalout1 < 0.05)
#Within-family model, CE: no clustering
summary(noclust2)

#Within-family model, CE: clustering on sibs
summary(wfout1)

#Within-family model: Dynastic Effects: clustering on sibs
summary(wfout2)
table(wfout2 < 0.05)
