#Packages
require(data.table)
require(forestplot)

#Read in data
data <- fread("Table-h2.txt", h = F , sep = " ")

#Define population/within-sibship estimates
data$WF <- rep(1:3, times=1, nrow(data))

Total <- data[which(data$WF == 1), ]
WF <- data[which(data$WF == 2), ]

#Cross-trait intercepts (used as estimate of covariance) when testing heterogeneity.
CTInt <- data[which(data$WF == 3), ]

#Generate confidence intervals
data$LCI <- data$V2 - (1.96* data$V3)
data$UCI <- data$V2 + (1.96* data$V3)

#Set label to phenotype names
label <- WF$V1

#Define parameters in figure
mean <- WF$V2
lci <- WF$LCI
uci <- WF$UCI

mean2 <- Total$V2
lci2 <- Total$LCI
uci2 <- Total$UCI

#Generate figure
tiff(filename="/mnt/storage/scratch/lh14833/Tables/Figure-h2.tiff", width=20, height=20, units="cm", res=600)

forestplot(label, 
           legend = c("Population", "Within-sibship"),
           line.margin = 0.3, # We need to add this to avoid crowding
           mean = cbind(mean2, mean),
           lower = cbind(lci2, lci),
           upper = cbind(uci2, uci),
           col=fpColors(box=c("blue", "darkred")),
           xlab="SNP h2 estimate (95% C.I.)",
           boxsize = 0.1,
          zero=0,
xticks=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

dev.off()
