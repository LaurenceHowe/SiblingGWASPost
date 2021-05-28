#Packages
require(data.table)
require(forestplot)

#Read in estimates
data <- fread("Table-h2_SumHer.txt")

#Define within-sibship and population (total) estimates
data$WF <- rep(1:2, times=1, nrow(data))
Total <- data[which(data$WF == 1), ]
WF <- data[which(data$WF == 2), ]

#Generate LCI and UCI from Beta/SE
data$LCI <- data$V3 - (1.96 * data$V4)
data$UCI <- data$V3 + (1.96 * data$V4)

#Label is phenotype name
label <- WF$V1

#Within-sibship estimates
mean <- WF$V3
lci <- WF$LCI
uci <- WF$UCI

#Population (Total) estimates
mean2 <- Total$V3
lci2 <- Total$LCI
uci2 <- Total$UCI

#Generate figure
tiff(filename="/mnt/storage/scratch/lh14833/Tables/Figure-h2_SumHer.tiff", width=20, height=20, units="cm", res=600)

forestplot(label, 
           legend = c("Population", "Within-sibship"),
           line.margin = 0.3, # We need to add this to avoid crowding
           mean = cbind(mean2, mean),
           lower = cbind(lci2, lci),
           upper = cbind(uci2, uci),
           col=fpColors(box=c("blue", "darkred")),
           xlab="SNP h2 estimate (95% C.I.)",
           boxsize = 0.1,
	   clip = c(-0.1, 1),
          zero=0,
xticks=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

dev.off()
