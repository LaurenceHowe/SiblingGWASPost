#Packages
require(data.table)
require(forestplot)

#Read in data
data <- fread("Education-rg.txt")

#Confidence intervals
data$LCI <- data$V3 - (1.96* data$V4)
data$UCI <- data$V3 + (1.96* data$V4)

#Sort out formatting for within-sibship/population estimates
WF <- data[which(data$V1 =="../WF" | data$V1 == "WF"), ]
Total <- data[which(data$V1 =="../Total" | data$V1 == "Total"), ]

#Label as phenotype names
label <- c("Height", "BMI", "EverSmk",  "SBP", "WHR", "Alcohol", "Menarche",
"Children", "Menopause", "Cognition", "Depressive",
"Neuroticism", "CPD", "LDL", "HDL", "TG", "CRP", "eGFR", "FEV1", "HbA1c", "Wellbeing")

#Specify parameters for the figure
mean <- WF$V3
lci <- WF$LCI
uci <- WF$UCI

mean2 <- Total$V3
lci2 <- Total$LCI
uci2 <- Total$UCI

#Generate figure
jpeg(filename="/mnt/storage/scratch/lh14833/Tables/Figure-education-rg.jpg", width=20, height=20, units="cm", res=600)
forestplot(label, 
           legend = c("Population", "Within-sibship"),
           line.margin = 0.3, # We need to add this to avoid crowding
           mean = cbind(mean2, mean),
           lower = cbind(lci2, lci),
           upper = cbind(uci2, uci),
           col=fpColors(box=c("blue", "darkred")),
           xlab="Genetic correlations with Educational attainment (95% C.I.)",
           boxsize = 0.1,
          clip=c(-1, 1),
          zero=0,
xticks=c(-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

dev.off()
