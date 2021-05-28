#Packages
require(data.table)
require(forestplot)
require(metafor)

#Read in data
data <- fread("Height_SDS.Studies.txt", h = F , sep = " ")

#Define population (Total) and within-sibship (WF) estimates
data$WF <- rep(1:2, times=1, nrow(data))
Total <- data[which(data$WF == 1), ]
WF <- data[which(data$WF == 2), ]

#Confidence intervals
data$LCI <- data$V2 - (1.96* data$V3)
data$UCI <- data$V2 + (1.96* data$V3)

#Subset to individual studies only for pooling (i.e. remove the meta-analysis) 
Total2 <- Total[c(1:8)]
WF2 <- WF[c(1:8)]

#Pooled meta-analysis of the 8 studies
res.fe.total <- metafor::rma(yi = Total2$V2, sei = Total2$V3, method="FE")
res.fe.wf <- metafor::rma(yi = WF2$V2, sei = WF2$V3, method ="FE")

#Extract metafor estimates
CombinedTotal <- data.frame(V1 = "Height", V2 = as.numeric(res.fe.total[2]), V3 = as.numeric(res.fe.total[3]), WF = 1, LCI = as.numeric(res.fe.total[6]), UCI = as.numeric(res.fe.total[7]))
CombinedWF <- data.frame(V1 = "Height", V2 = as.numeric(res.fe.wf[2]), V3 = as.numeric(res.fe.wf[3]), WF = 2, LCI = as.numeric(res.fe.wf[6]), UCI = as.numeric(res.fe.wf[7]))

#Merge with the rest of the data
Total3 <- rbind(Total2, CombinedTotal, Total[9])
WF3 <- rbind(WF2, CombinedWF, WF[9])

#Study labels
label <- c("UKBiobank", "HUNT", "Generation Scotland", "DiscovEHR", "QIMR", "FinnTwin", "Netherlands Twin Registry", "China Kadoorie Biobank", "Pooled study-level estimate", "European meta-analysis GWAS")

#Multiply estimates by -1 to get in the direction of selection pressure on trait increasing alleles
mean <- WF3$V2 * (-1)
lci <- WF3$UCI * (-1)
uci <- WF3$LCI * (-1)

mean2 <- Total3$V2 * (-1)
lci2 <- Total3$UCI * (-1)
uci2 <- Total3$LCI * (-1)

#Generate figure
tiff(filename="/mnt/storage/scratch/lh14833/Tables/Height_SDS_Studies.tiff", width=20, height=10, units="cm", res=600)

forestplot(label, 
           legend = c("Population", "Within-sibship"),
           line.margin = 0.3, # We need to add this to avoid crowding
           mean = cbind(mean2, mean),
           lower = cbind(lci2, lci),
           upper = cbind(uci2, uci),
           col=fpColors(box=c("blue", "darkred")),
           xlab="Spearman rank correlation for height (95% C.I.)",
            boxsize = 0.1,
          zero=0,
xticks=c(-0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04))

dev.off()
