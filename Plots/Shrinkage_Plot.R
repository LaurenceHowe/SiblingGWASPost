#Packages
require(data.table)
require(forestplot)

#Read in data
data <- fread("Shrinkage_results.txt")

#Two genetic variant inclusion thresholds: genome-wide significant and more liberal threshold
GWS <- data[which(data$V2 == "GWS"), ]
NonGWS <- data[! which(data$V2 == "GWS"), ]

#Set label to phenotype name
label <- GWS$V1

#Multiply estimates by 100 to get % shrinkage for GWS/nonGWS
mean <- GWS$V7 * 100
lci <- GWS$V8 * 100
uci <- GWS$V9 * 100

mean2 <- NonGWS$V7 * 100
lci2 <- NonGWS$V8 * 100
uci2 <- NonGWS$V9 * 100

#Generate figure

jpeg("/mnt/storage/scratch/lh14833/Tables/Figure2_shrinkage.jpg", width = 20, height = 20, units = "cm", res = 600)

forestplot(label, 
           legend = c("S_G", "S_L"),
           line.margin = 0.3, # We need to add this to avoid crowding
           mean = cbind(mean, mean2),
           lower = cbind(lci, lci2),
           upper = cbind(uci, uci2),
           col=fpColors(box = c("darkblue", "orange")),
	         lwd.zero = 2,
           xlab = "Percentage shrinkage in effect sizes (95% C.I.)",
	         clip = c(-50, 100),
           boxsize = 0.1,
          zero = 0,
xticks = c(-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

dev.off()
