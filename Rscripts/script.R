#### set working directory
### Tms ###
### Tge ###
setwd("/home/susana/Dropbox/Timema_cryptic_geneflow/genevievae/FB_scaf")


### change scaffold name to cr name
### Tms ###
scaf <- read.csv2(file="/home/susana/Dropbox/Timema_cryptic_geneflow/monikensis/trimmed/Tms.anchored.scafs", sep = "\t")
tms <- read.csv2(file="/home/susana/Dropbox/Timema_cryptic_geneflow/monikensis/trimmed/Tms.fb.bial3.table", sep = "\t")
tms2 <- read.csv2(file="/home/susana/Dropbox/Timema_cryptic_geneflow/monikensis/trimmed/Tms.trim2.bial3.table", sep = "\t")
tms3 <- read.csv2(file="/home/susana/Dropbox/Timema_cryptic_geneflow/monikensis/trimmed/Tms.trim3.bial3.table", sep = "\t")
colnames(tms)
### Tge ###
scaf <- read.csv2(file="/home/susana/Dropbox/Timema_cryptic_geneflow/genevievae/FB_scaf/5_Tge_anchored_scfs.tsv", sep = "\t")
tge <- read.csv2(file="/home/susana/Dropbox/Timema_cryptic_geneflow/genevievae/FB_scaf/Tge.trim3.test.bial3.table", sep = "\t")

## change only one col name
### Tms ###
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Tms9_P2	Tms9_P1	Tms8_P1	Tms7_P2	Tms6_P1	Tms4_P2	Tms22_P1	Tms8_P2	Tms4_P1	Tms22_P2	Tms20_P1	Tms2_P1	Tms14_P1	Tms23_P2	Tms12_P1	Tms14_P2	Tms12_P2	Tms11_P1	Tms3_P1	Tms21_P2	Tms10_P1	Tms5_P1	Tms23_P1	Tms13_P1	Tms10_P2	Tms7_P1	Tms15_P2	Tms21_P1Tms13_P2	Tms1_P1	Tms11_P2	Tms15_P1	Tms5_P2	Tms16_P1	Tms17_P1	Tms6_P2	Tms17_P2	Tms19_P1	Tms1_P2	Tms18_P1	Tms2_P2	Tms18_P2	Tms3_P2	Tms20_P2	Tms16_P2	Tms19_P2

### Tge ###
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Gepop2.7	Gepop2.9	Gepop2.4	Gepop2.2	Gepop2.18	Gepop2.17	Gepop2.6	Gepop2.16	Gepop2.15	Gepop2.11	Gepop2.10	Gepop2.14	Gepop1.12	Gepop2.12	Gepop1.19	Gepop1.17	Gepop2.22	Gepop1.21	Gepop2.1	Gepop1.15	Gepop2.23	Gepop2.20	Gepop2.13	Gepop1.14	Gepop1.9	Gepop1.20	Gepop1.5	Gepop1.7	Gepop1.18	Gepop1.2	Gepop2.21	Gepop1.10	Gepop1.1	Gepop2.5	Gepop1.11	Gepop2.8	Gepop1.13	Gepop1.22	Gepop2.19	Gepop1.16	Gepop1.23	Gepop2.3	Gepop1.3	Gepop1.4	Gepop1.6	Gepop1.8

### Tms ###
colnames(tms)[1] <- "scafname"
colnames(tms2)[1] <- "scafname"
colnames(tms3)[1] <- "scafname"
colnames(scaf)[1] <- "scafname"
### Tge ###
colnames(tge)[1] <- "scafname"
colnames(scaf)[1] <- "scafname"



## merge and keep in the table file only the scaffolds with cr (lg) correspondance
### Tms ###
test <- merge(tms, scaf, by='scafname')
test2 <- merge(tms2, scaf, by='scafname')
test3 <- merge(tms3, scaf, by='scafname')
colnames(test)
### Tge ###
test <- merge(tge, scaf, by='scafname')

#test$paste <- paste(test$POS, '_', test$len, sep='')
#test2$paste <- paste(test2$POS, '_', test2$len, sep='')
#test3$paste <- paste(test3$POS, '_', test3$len, sep='')

### reorder columns
### Tms ###
test.srt <- test[with(test, order(lg, scafname, POS)), ]
test2.srt <- test2[with(test2, order(lg, scafname, POS)), ]
test3.srt <- test3[with(test3, order(lg, scafname, POS)), ]
head(test.srt)
colnames(test.srt)
### Tge ###
test.srt <- test[with(test, order(lg, scafname, POS)), ]

# prepare file to print
### Tms ###
tmscr <- test.srt[,c(58,2,56,4:55)]
tmscr2 <- test2.srt[,c(58,2,56,4:55)]
tmscr3 <- test3.srt[,c(58,2,56,4:55)]
colnames(tmscr)
colnames(tmscr)[1] <- "#CHROM"
colnames(tmscr)[2] <- "POS"
colnames(tmscr)[3] <- "ID"
colnames(tmscr2)[1] <- "#CHROM"
colnames(tmscr2)[2] <- "POS"
colnames(tmscr2)[3] <- "ID"
colnames(tmscr3)[1] <- "#CHROM"
colnames(tmscr3)[2] <- "POS"
colnames(tmscr3)[3] <- "ID"
### Tge ###
tgecr <- test.srt[,c(58,2,3,4:55)]
colnames(tgecr)
colnames(tgecr)[1] <- "#CHROM"
colnames(tgecr)[2] <- "POS"
colnames(tgecr)[3] <- "ID"

## print file to convert to vcf
### Tms ###
write.table(tmscr, file="Tms.fb.lg.txt", sep ="\t", quote=FALSE, row.names=FALSE)
write.table(tmscr2, file="Tms.trim2.lg.txt", sep ="\t", quote=FALSE, row.names=FALSE)
write.table(tmscr3, file="Tms.trim3.lg.txt", sep ="\t", quote=FALSE, row.names=FALSE)
### Tge ###
write.table(tgecr, file="Tge.lg.txt", sep ="\t", quote=FALSE, row.names=FALSE)

### ouside R

## open file "Tge.lg.txt" on excel and run this function on a new column
# after pasting first value per lg group on that column.
# =IF(B3>B2,(C2+(B3-B2)), (C2 + B3))
# delete old POS column and make the new column, POS
# save it and make this new corrected vcf (TEMP until Kamil gives me corrected lg coordinates)

### Tms ###
grep -v '^##' Tms.fb.bial3.vcf > Tms.fb.bial3.table
grep -v '^##' Tms.trim2.bial3.vcf > Tms.trim2.bial3.table
grep -v '^##' Tms.trim3.bial3.vcf > Tms.trim3.bial3.table
grep '^##' Tms.fb.bial3.vcf > fb.header
cat fb.header Tms.fb.lg.txt > Tms.fb.lg.vcf
grep '^##' Tms.trim2.bial3.vcf > trim2.header
cat trim2.header Tms.trim2.lg.txt > Tms.trim2.lg.vcf
grep '^##' Tms.trim3.bial3.vcf > trim3.header
cat trim3.header Tms.trim3.lg.txt > Tms.trim3.lg.vcf
### Tge ###
grep '^##' Tge.trim3.test.bial3.vcf > header
cat header Tge.lg.txt > Tge.lg.vcf


### inside R again
# make het plots per individual
hetfb <- read.csv2(file = "Tms.fb.het", header=T, sep='\t')
het2 <- read.csv2(file = "Tms.trim2.het", header=T, sep='\t')
het3 <- read.csv2(file = "Tms.trim3.het", header=T, sep='\t')

colnames(hetfb) <- c("INDV", "O.HOM", "E.HOM", "N_SITES", "F", "indvno", "POP")
colnames(het2) <- c("INDV", "O.HOM", "E.HOM", "N_SITES", "F", "indvno", "POP")
colnames(het3) <- c("INDV", "O.HOM", "E.HOM", "N_SITES", "F", "indvno", "POP")
hetfb$ratio <- hetfb$O.HOM/hetfb$N_SITES
het2$ratio <- het2$O.HOM/het2$N_SITES
het3$ratio <- het3$O.HOM/het3$N_SITES


## plot by indv order
hetfb$INDV <- factor(hetfb$INDV, levels = hetfb$INDV)
het2$INDV <- factor(het2$INDV, levels = het2$INDV)
het3$INDV <- factor(het3$INDV, levels = het3$INDV)

## rbind all datasets to plot in ggplot
all <- rbind(hetfb, het2, het3)
all$trim <- c(rep("fb", 46), rep("t2", 46), rep("t3", 46))

library(ggplot2)
pdf(file="het-all-DPcorrected.pdf", width=10)
ggplot(data = all, mapping = aes(x = INDV, y = ratio, colour=POP)) +
	geom_point() +
	facet_grid(trim ~ .) +
	theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8))
dev.off()



#### calculate hom1 hom2 and het relationships

hwefb <- read.csv2(file = "Tms.fb.srt.hwe", header=T, sep='\t')
hwe2 <- read.csv2(file = "Tms.trim2.srt.hwe", header=T, sep='\t')
hwe3 <- read.csv2(file = "Tms.trim3.srt.hwe", header=T, sep='\t')




##### plot hom-het ratio

allhet <- read.csv2(file = "all.het", header=T, sep='\t')

# "sample" - sample ID
# "X1.0" -> "het1.0" - no of heterozygous positions GT = "1/0"
# "X0.1" -> "het0.1" - no of heterozygous positions GT = "1/0"
# "X1.0.1" -> "het1..0" - no of heterozygous positions GT = "1/0"
# "X0.1.1" -> "het0..1" - no of heterozygous positions GT = "1/0"
# "trim" -> "trim" - the original vcf file trimming treatment (fb, trim2 or trim3)
# "sum.het." -> "sum.het" - sum of total no of heterozygous positions
#  "X1..1." -> "hom1.1" - no of homozygous positions GT = "1/0"
#  "X0..0" -> "hom0.0" - no of homozygous positions GT = "1/0"
# "X1.1" -> "hom1..1" - no of homozygous positions GT = "1/0"
# "X0.0" ->  "hom0..0" - no of homozygous positions GT = "1/0"
#  "sum.hom." - sum of total no of homozygous positions
# "sum.homnonref." - sum of total no of homozygous positions EXCEPT the reference homozygous (0/0 and 0|0)
# "het.hom" - ratio of sum het and sum hom
# "het.homnonref" - ratio of sum het and sum hom (without reference hom)
# "Vcflib.hethomration" - ratio het hom calculated with vcflib
# "het..het.hom." -> "het.total" - ratio het positions and sum (hom + het)
# "het..het.homnonref." -> "het.totalnonref" - - ratio het positions and sum (hom without ref hom + het)
# "hom..hom.het." -> "hom.total" - ratio hom and total (hom + het)
# "X" ->  "hom.totalnonref" - - ratio hom without reference hom and total (hom without ref hom + het)
# "sampleno" - no of the sample
# "POP" - population where the sample came from



colnames(allhet) <- c("sample", "het1.0", "het0.1", "het1..0", "het0..1", "trim", "sum.het",
			"hom1.1", "hom0.0", "hom1..1", "hom0..0", "sum.hom", "sum.homnonref",
			"het.hom", "het.homnonref", "Vcflib.hethomratio", "het.total",
			"het.totalnonref", "hom.total", "hom.totalnonref", "sampleno", "POP")
str(allhet)
## some variables were called as factor, and they should be called as numerical:
allhet$het.hom <- as.numeric(as.character(allhet$het.hom))
allhet$het.homnonref <- as.numeric(as.character(allhet$het.homnonref))
allhet$Vcflib.hethomratio <- as.numeric(as.character(allhet$Vcflib.hethomratio))
allhet$het.total <- as.numeric(as.character(allhet$het.total))
allhet$het.totalnonref <- as.numeric(as.character(allhet$het.totalnonref))
allhet$hom.total <- as.numeric(as.character(allhet$hom.total))
allhet$hom.totalnonref <- as.numeric(as.character(allhet$hom.totalnonref))

## keep in the plot the same order as in the dataset
allhet$sample <-factor(allhet$sample, levels=unique(allhet$sample))

pdf(file="tms-het-vcflibratio.pdf", width = 10)
ggplot(data = allhet, mapping = aes(x = sample, y = Vcflib.hethomratio, colour=POP)) +
	geom_point() +
	facet_grid(trim ~ .) +
	theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8)) #+
#	ylim(0.0,0.10)
dev.off()

pdf(file="tms-het-total.homnonref.ratio.pdf", width = 10)
ggplot(data = allhet, mapping = aes(x = sample, y = het.totalnonref, colour=POP)) +
	geom_point() +
	facet_grid(trim ~ .) +
	theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8)) #+
#	ylim(0.0,0.10)
dev.off()

pdf(file="tms-het-total.ratio-small.pdf", width = 10)
ggplot(data = allhet, mapping = aes(x = sample, y = het.total, colour=POP)) +
	geom_point() +
	facet_grid(trim ~ .) +
	theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8)) +
	ylim(0.0,0.03)
dev.off()





######## LOOKING FOR POLYPLOIDS #######
######## READ AO AND RO values ########
ratio <- read.csv2("tms.trim3.lg.DPratio.txt", sep=",", header=T)

library(ggplot2)
library(reshape2)

melted <- melt(ratio, id.vars = c('X.CHROM','POS'))
## how to remove all NAs
ratioclean <- na.omit(melted)


## read column value as numeric
ratioclean$value <- as.numeric(as.character(ratioclean$value))

## set the order of samples
ratioclean$variable <- factor(ratioclean$variable, levels = c("Tms1_P1","Tms2_P1","Tms3_P1","Tms4_P1","Tms5_P1","Tms6_P1","Tms7_P1",
			"Tms8_P1","Tms9_P1","Tms10_P1","Tms11_P1","Tms12_P1","Tms13_P1","Tms14_P1","Tms15_P1","Tms16_P1","Tms17_P1",
			"Tms18_P1","Tms19_P1","Tms20_P1","Tms21_P1","Tms22_P1","Tms23_P1","Tms1_P2","Tms2_P2","Tms3_P2","Tms4_P2",
			"Tms5_P2","Tms6_P2","Tms7_P2","Tms8_P2","Tms9_P2","Tms10_P2","Tms11_P2","Tms12_P2","Tms13_P2","Tms14_P2",
			"Tms15_P2","Tms16_P2","Tms17_P2","Tms18_P2","Tms19_P2","Tms20_P2","Tms21_P2","Tms22_P2","Tms23_P2"))

### make a barplot of counts: count the number of positions per ratio class
# Instead of faceting with a variable in the horizontal or vertical direction,
# facets can be placed next to each other, wrapping with a certain number of columns or rows.
# The label for each plot will be at the top of the plot.
pdf(file="tms.trim3.plot_ratio_DP.pdf", width=10)
ggplot(ratioclean, aes(x=value)) +
	geom_histogram(binwidth=.05, stat = "count") +
	facet_wrap(~variable, ncol=6) +  
	xlim(c(0.1,1)) + ylim(c(0, 100)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
dev.off()


## subset data and plot only a few samples
### WEIRD SAMPLES
tms19_P2 <- ratioclean[ which(ratioclean$variable=='Tms19_P2'), ]
tms1_P2 <- ratioclean[ which(ratioclean$variable=='Tms1_P2'), ]
tms22_P2 <- ratioclean[ which(ratioclean$variable=='Tms22_P2'), ]
tms17_P2 <- ratioclean[ which(ratioclean$variable=='Tms17_P2'), ]

## set order of lg
tms19_P2$X.CHROM <- factor(tms19_P2$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))
tms1_P2$X.CHROM <- factor(tms1_P2$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))
tms22_P2$X.CHROM <- factor(tms22_P2$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))
tms17_P2$X.CHROM <- factor(tms17_P2$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))
tms19_P1$X.CHROM <- factor(tms19_P1$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))
tms1_P1$X.CHROM <- factor(tms1_P1$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))
tms16_P2$X.CHROM <- factor(tms16_P2$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))
tms21_P2$X.CHROM <- factor(tms21_P2$X.CHROM, levels = c("lg1","lg2","lg3","lg4","lg5","lg6","lg7", "lg8", "lg9",
							 "lg10", "lg11", "lg12", "lgX"))


### WEIRD SAMPLE: Tms19_P2
pdf(file="tms.trim3.tms19_P2_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms19_P2_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms19_P2_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms19_P2_ratio_DP-density3.pdf", width=10)
ggplot(tms19_P2, aes(x=value)) +
	geom_histogram(binwidth=.1, stat = "count") +
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	xlim(c(0.1,1)) + 
#	ylim(c(0, 70)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms19_P2) + 
#	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.1, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + xlim(c(0.1,1)) #+ 
#	facet_wrap(~X.CHROM, ncol=5, scales="free_y") 
dev.off()

### WEIRD SAMPLE: Tms1_P2
pdf(file="tms.trim3.tms1_P2_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms1_P2_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms1_P2_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms1_P2_ratio_DP-density3.pdf", width=10)
ggplot(tms1_P2, aes(x=value)) +
	geom_histogram(binwidth=.1, stat = "count") +
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	xlim(c(0.1,1)) + 
#	ylim(c(0, 70)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms1_P2) + 
	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.1, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + xlim(c(0.1,1)) + 
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") 
dev.off()

### WEIRD SAMPLE: Tms22_P2
pdf(file="tms.trim3.tms22_P2_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms22_P2_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms22_P2_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms22_P2_ratio_DP-density3.pdf", width=10)
ggplot(tms22_P2, aes(x=value)) +
	geom_histogram(binwidth=.1, stat = "count") +
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	xlim(c(0.1,1)) + 
#	ylim(c(0, 70)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms22_P2) + 
	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.1, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + xlim(c(0.1,1)) + 
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") 
dev.off()

### WEIRD SAMPLE: Tms17_P2
pdf(file="tms.trim3.tms17_P2_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms17_P2_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms17_P2_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms17_P2_ratio_DP-density3.pdf", width=10)
ggplot(tms17_P2, aes(x=value)) +
	geom_histogram(binwidth=.1, stat = "count") +
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	xlim(c(0.1,1)) + 
#	ylim(c(0, 70)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms17_P2) + 
	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.1, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + xlim(c(0.1,1)) + 
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") 
dev.off()



### NORMAL SAMPLES (as examples)
tms19_P1 <- ratioclean[ which(ratioclean$variable=='Tms19_P1'), ]
tms1_P1 <- ratioclean[ which(ratioclean$variable=='Tms1_P1'), ]
tms21_P2 <- ratioclean[ which(ratioclean$variable=='Tms21_P2'), ]
tms16_P2 <- ratioclean[ which(ratioclean$variable=='Tms16_P2'), ]
tms19_P1$value <- as.numeric(as.character(tms19_P1$value))
### NORMAL SAMPLE: Tms19_P1
pdf(file="tms.trim3.tms19_P1_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms19_P1_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms19_P1_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms19_P1_ratio_DP-whylg12.pdf", width=10)
pdf(file="tms.trim3.tms19_P1_ratio_DP-whylg12-2.pdf", width=10)
ggplot(tms19_P1, aes(x=value)) +
	geom_histogram(aes(fill = tms19_P1$X.CHROM), binwidth=.1, stat = "bin") +
#	stat_bin(mapping = NULL, data = NULL, geom = "bar",
#    position = "stack", width = 0.9, drop = FALSE,
#    right = FALSE, binwidth = NULL, origin = NULL,
#    breaks = NULL, ...)
#	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	facet_wrap(~X.CHROM, ncol=5, scales="free") + 
#	xlim(c(0,1)) + 
#	ylim(c(0, 100)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms19_P1) + 
#	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.2, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + ylim(c(0,5)) + 
	facet_wrap(~X.CHROM, ncol=5, scales="free_x") 
dev.off()

### NORMAL SAMPLE: Tms1_P1
pdf(file="tms.trim3.tms1_P1_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms1_P1_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms1_P1_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms1_P1_ratio_DP-whylg12.pdf", width=10)
ggplot(tms1_P1, aes(x=value)) +
	geom_histogram(aes(fill = tms1_P1$X.CHROM), binwidth=.1, stat = "bin") +
#	stat_bin(mapping = NULL, data = NULL, geom = "bar",
#    position = "stack", width = 0.9, drop = FALSE,
#    right = FALSE, binwidth = NULL, origin = NULL,
#    breaks = NULL, ...)
#	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	facet_wrap(~X.CHROM, ncol=5, scales="free") + 
#	xlim(c(0,1)) + 
#	ylim(c(0, 100)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms1_P1) + 
	geom_histogram(aes(x = value, fill = tms1_P1$X.CHROM), binwidth=.2, stat = "count") + 
#	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) +
#	xlim(c(0,1)) + 
	ylim(c(0,20)) +
	facet_wrap(~X.CHROM, ncol=5, scales="free_x") 
dev.off()

### NORMAL SAMPLE: Tms21_P2
pdf(file="tms.trim3.tms21_P2_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms21_P2_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms21_P2_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms21_P2_ratio_DP-density3.pdf", width=10)
ggplot(tms21_P2, aes(x=value)) +
	geom_histogram(binwidth=.1, stat = "count") +
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	xlim(c(0.1,1)) + 
#	ylim(c(0, 70)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms21_P2) + 
	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.1, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + xlim(c(0.1,1)) + 
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") 
dev.off()

### NORMAL SAMPLE: Tms16_P2
pdf(file="tms.trim3.tms16_P2_ratio_DP.pdf", width=10)
pdf(file="tms.trim3.tms16_P2_ratio_DP-density1.pdf", width=10)
pdf(file="tms.trim3.tms16_P2_ratio_DP-density2.pdf", width=10)
pdf(file="tms.trim3.tms16_P2_ratio_DP-density3.pdf", width=10)
ggplot(tms16_P2, aes(x=value)) +
	geom_histogram(binwidth=.1, stat = "count") +
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") +  
	xlim(c(0.1,1)) + 
#	ylim(c(0, 70)) +
	theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0),axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"))
#
ggplot(tms16_P2) + 
	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.1, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + xlim(c(0.1,1)) + 
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") 
dev.off()







df <- with(Galton, data.frame(x = x, y = dnorm(x, mean(parent), sd(parent))))
plot(density(ratioclean$value))
ggplot(tms17_P2) + 
	geom_histogram(aes(x = value, fill = X.CHROM), binwidth=.1, stat = "count") + 
	geom_density(aes(x = value, fill = X.CHROM), alpha = 0.2) + xlim(c(0.1,1)) + 
	facet_wrap(~X.CHROM, ncol=5, scales="free_y") 
### calculate mean sd etc of distribution per lg
summary(tms1_P1)












newdata=data.frame(matrix(NA, nrow = nrow(AO), ncol = ncol(AO)))
for (line in nrow(AO)) {
	for (column in length(AO)) {
		if (AO[line, column] >= RO[line,column]) { ratio=AO[line,column]/RO[line,column] }
		else if (AO[line, column] < RO[line,column]) { ratio=RO[line,column]/AO[line,column] }
		newdata[line,column]=ratio
		}
	}				
	



######################################################################################
########################### plot heterozygosities ####################################
######################################################################################

setwd("C:/Users/User/Dropbox/Timema_cryptic_geneflow/mapping_scaf_coordinates_new/test_plink/allinds/heterozygosity")

ls()


colnames(lg1) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", "hom11", "hom00", 
  "het11", "het00",	"sum(hom)",	"sum(homnonref)",	"het/hom", "het/homnonref", "Vcflib-hethomration", "het/(het+hom)",
    "het/(het+homnonref)", "hom/(hom+het)", "homnonref/(homnonref+het)", "sampleno", "POP")


#### Tge ####
lg1 <- read.csv2(file = "Tms_lg1_coord.csv", header=T, sep='\t')
lg2 <- read.csv2(file = "Tms_lg2_coord.csv", header=T, sep='\t')
lg3 <- read.csv2(file = "Tms_lg3_coord.csv", header=T, sep='\t')
lg4 <- read.csv2(file = "Tms_lg4_coord.csv", header=T, sep='\t')
lg5 <- read.csv2(file = "Tms_lg5_coord.csv", header=T, sep='\t')
lg6 <- read.csv2(file = "Tms_lg6_coord.csv", header=T, sep='\t')
lg7 <- read.csv2(file = "Tms_lg7_coord.csv", header=T, sep='\t')
lg8 <- read.csv2(file = "Tms_lg8_coord.csv", header=T, sep='\t')
lg9 <- read.csv2(file = "Tms_lg9_coord.csv", header=T, sep='\t')
lg10 <- read.csv2(file = "Tms_lg10_coord.csv", header=T, sep='\t')
lg11 <- read.csv2(file = "Tms_lg11_coord.csv", header=T, sep='\t')
lg12 <- read.csv2(file = "Tms_lg12_coord.csv", header=T, sep='\t')
lgX <- read.csv2(file = "Tms_lgX_coord.csv", header=T, sep='\t')


### change colnames to something readable and writeable
colnames(lg1) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")


colnames(lg2) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg3) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg4) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg5) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg6) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg7) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg8) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg9) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")

colnames(lg10) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                    "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                    "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                    "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg11) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                    "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                    "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                    "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lg12) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                    "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                    "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                    "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")
colnames(lgX) <- c("sample", "het10", "het01", "het10phs",	"het01phs", "trim", "sum(het)", 
                   "hom11", "hom00", "het11", "het00",	"sumhom",	"sumhomnonref",	"het_hom",
                   "het_homnonref", "Vcflib-hethomration", "het_hethom", "het_hethomnonref",
                   "hom_homhet", "homnonref_homnonrefhet", "sampleno", "POP")


## add new column with lg no
lg1$lg <- rep("lg1",times = 46)
lg2$lg <- rep("lg2",times = 46)
lg3$lg <- rep("lg3",times = 46)
lg4$lg <- rep("lg4",times = 46)
lg5$lg <- rep("lg5",times = 46)
lg6$lg <- rep("lg6",times = 46)
lg7$lg <- rep("lg7",times = 46)
lg8$lg <- rep("lg8",times = 46)
lg9$lg <- rep("lg9",times = 46)
lg10$lg <- rep("lg10",times = 46)
lg11$lg <- rep("lg11",times = 46)
lg12$lg <- rep("lg12",times = 46)
lgX$lg <- rep("lgX",times = 46)

### concatenate all veritcally
Tgecat <- do.call("rbind", list(lg1, lg2, lg3, lg4, lg5, lg6, lg7, lg8, lg9, lg10, lg11, lg12, lgX))

## plot by indv order
Tgecat$sampleno <- factor(Tgecat$sampleno, levels = c("Gepop1.1", "Gepop1.2", "Gepop1.3", "Gepop1.4",
                                                  "Gepop1.5", "Gepop1.6", "Gepop1.7", "Gepop1.8" ,
                                                   "Gepop1.9", "Gepop1.10", "Gepop1.11", "Gepop1.12",
                                                 "Gepop1.13", "Gepop1.14", "Gepop1.15", "Gepop1.16",
                                                  "Gepop1.17", "Gepop1.18", "Gepop1.19", "Gepop1.20",
                                                  "Gepop1.21", "Gepop1.22", "Gepop1.23"))
Tgecat$sampleno <- factor(Tgecat$sampleno, levels = c("1", "2", "3", "4", "5", "6", "7", 
                                                      "8", "9", "10", "11", "12", "13", "14", "15", 
                                                      "16", "17", "18", "19", "20", "21", "22", "23"))
## convert several columns as numeric
cols = c(2:5, 7:15, 17:20);    
Tgecat[,cols] = apply(Tgecat[,cols], 2, function(x) as.numeric(as.character(x)))



library(ggplot2)
#install.packages("ggplot2")
pdf(file="het-all-DPcorrected.pdf", width=10)
ggplot(data = Tgecat, mapping = aes(x = lg, y = Tgecat$het_hethom, colour = POP)) +
  geom_point() +
#  facet_grid(sample ~ .) +
  theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8))
dev.off()

ggplot(data = Tgecat, mapping = aes(x = lg, y = het_hethom, colour=POP)) +
  geom_boxplot() +
  #  facet_grid(sample ~ .) +
  theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8))

pdf(file="Tge_ratio_het_all.pdf", width=10)
ggplot(data = Tgecat, mapping = aes(x = lg, y = het_hethom, colour=POP)) +
  geom_boxplot() +
  #  facet_grid(sample ~ .) +
  theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8)) +
  ylim(0, 0.05)
dev.off()



ggplot(data = Tgecat, mapping = aes(x = lg, y = hom_homhet, colour=POP)) +
  geom_boxplot() +
  #  facet_grid(sample ~ .) +
  theme(axis.text.x = element_text(angle = 90,  vjust = -0.01, size=8)) +
  ylim(0.9, 1.0)


#######################################################
################ ANALYSE LD RESULTS ###################
#######################################################
setwd("C:/Users/User/Dropbox/Timema_cryptic_geneflow/mapping_scaf_coordinates_new/test_plink/allinds/heterozygosity")
setwd("C:/Users/User/Dropbox/Timema_cryptic_geneflow/simulations/files_plot_R")

ls()
LD <- read.delim(file = "test.geno.ld", header = TRUE)
ggplot(LD, aes(x = ((LD$POS2 + LD$POS1)/2), y = LD$R^2)) + geom_point()

##### in bash


# working file: Tge_alignedcoords_P1.maf.inter.ld
sed -i 's/ /+/g' Tge_alignedcoords_P2.maf.inter.ld
# remove the last character
sed -i 's/.$//' Tge_alignedcoords_P2.maf.inter.ld

sed -i 's/+++++++++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld 
sed -i 's/+++++++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld 
sed -i 's/++++++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld

sed -i 's/+++++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/+++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/+++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld

sed -i 's/+++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/++++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/+++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/++/\t/g' Tge_alignedcoords_P2.maf.inter.ld
sed -i 's/+/\t/g' Tge_alignedcoords_P2.maf.inter.ld

# some '.' were left without tabs between them and the next column. To correct that:
sed -i 's/\.l/\.\tl/g' Tge_alignedcoords_P2.maf.inter.ld

# remove \t in the beginning of the line
sed -i 's/^\t//g' Tge_alignedcoords_P2.maf.inter.ld




# make violin plot
#install.packages("ggplot2", lib="/users/sfreitas1/R")
#install.packages("crayon", lib="/users/sfreitas1/R")
#install.packages("vctrs", lib="/users/sfreitas1/R")
#install.packages("backports", lib="/users/sfreitas1/R")
#install.packages("withr", lib="/users/sfreitas1/R")
#install.packages("labeling", lib="/users/sfreitas1/R")
#install.packages("digest", lib="/users/sfreitas1/R")
#install.packages("dplyr", lib="/users/sfreitas1/R")
#install.packages("fansi", lib="/users/sfreitas1/R")
#install.packages("utf8", lib="/users/sfreitas1/R")
#install.packages("cli", lib="/users/sfreitas1/R")

## to be used in wally
#library("cli", lib.loc="/users/sfreitas1/R")
#library("utf8", lib.loc="/users/sfreitas1/R")
#library("fansi", lib.loc="/users/sfreitas1/R")
#library("digest", lib.loc="/users/sfreitas1/R")
#library("labeling", lib.loc="/users/sfreitas1/R")
#library("withr", lib.loc="/users/sfreitas1/R")
#library("backports", lib.loc="/users/sfreitas1/R")
#library("vctrs", lib.loc="/users/sfreitas1/R")
#library("crayon", lib.loc="/users/sfreitas1/R")
#library("ggplot2", lib.loc="/users/sfreitas1/R")
#library("dplyr", lib.loc="/users/sfreitas1/R")

library("cli")
library("utf8")
library("fansi")
library("digest")
library("labeling")
library("withr")
library("backports")
library("vctrs")
library("crayon")
library("ggplot2")
library("dplyr")


# read in the readable file (after sed/bash tweaking)
# asex pops
p1in <- read.delim("Tge_alignedcoords_P1.maf.inter.ld", header=T)
p2in <- read.delim("Tge_alignedcoords_P2.maf.inter.ld", header=T)
p1in <- read.delim("Tms_alignedcoords_P1.maf.inter.ld", header=T)
p2in <- read.delim("Tms_alignedcoords_P2.maf.inter.ld", header=T)
# changing the collumns names
colnames(p1in) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
colnames(p2in) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")


# sex pop
tce <- read.delim("Tce_alignedcoords.lg.maf.inter.ld", header=T)
colnames(tce) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")



# simulations
# cloning - no recombination
clono <- read.delim("cloning-norecomb.lg.maf.inter.ld")
colnames(clono) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
# cloning - low recombination
clolo <- read.delim("cloning-lowrecomb.lg.maf.inter.ld")
colnames(clolo) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")
# cloning - high recombination
clohi <- read.delim("cloning-highrecomb.lg.maf.inter.ld")
colnames(clohi) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")

# selfing - no recombination
selno <- read.delim("selfing-norecomb.lg.maf.inter.ld")
# selfing - low recombination
sellow <- read.delim("selfing-lowrecomb.lg.maf.inter.ld")
# selfing - high recombination
selhi <- read.delim("selfing-highrecomb.lg.maf.inter.ld")


# make column with chr-chr
p1in$inter <- paste(p1in$CHR_A,'-',p1in$CHR_B, sep='')
p2in$inter <- paste(p2in$CHR_A,'-',p2in$CHR_B, sep='')

clono$inter <- paste(clono$CHR_A,'-',clono$CHR_B, sep='')
clolo$inter <- paste(clolo$CHR_A,'-',clolo$CHR_B, sep='')
clohi$inter <- paste(clohi$CHR_A,'-',clohi$CHR_B, sep='')

selno$inter <- paste(selno$CHR_A,'-',selno$CHR_B, sep='')
sello$inter <- paste(sello$CHR_A,'-',sello$CHR_B, sep='')
selhi$inter <- paste(selhi$CHR_A,'-',selhi$CHR_B, sep='')

tce$inter <- paste(tce$CHR_A,'-',tce$CHR_B, sep='')


# confirm if we are doing it correctly
chr_unq <- unique(p1in$inter)
length(chr_unq)
# 91 
### GOOD!! continue :)

# to make subset of data (example)
lg1 <- subset(p1in, CHR_A %in% c("lg1"))
lg1_P2 <- subset(p2in, CHR_A %in% c("lg1"))
## how many cr interactions we have
"lg1-lg1", "lg1-lg2", "lg1-lg3", "lg1-lg4", "lg1-lg5", "lg1-lg6",
"lg1-lg7", "lg1-lg8", "lg1-lg9", "lg1-lg10", "lg1-lg11", "lg1-lg12", "lg1-lgX"

## example of the violin plot
pdf(file="lg1-all.pdf")
myplot <- ggplot(lg1) + 
	geom_violin(aes(x= inter, y = R2, colour=inter)) +
	facet_wrap(~inter, scales = "free")
myplot <- ggplot(lg1_P2) + 
	geom_violin(aes(x= inter, y = R2, colour=inter)) +
	facet_wrap(~inter, scales = "free")
print(myplot)
dev.off()

## reorder x axis labels
clohi$inter2 <- factor(clohi$inter, levels = c("lg1-lg1", "lg2-lg2", "lg3-lg3", "lg4-lg4", "lg5-lg5", "lg6-lg6", "lg7-lg7", "lg8-lg8",
				"lg9-lg9", "lg10-lg10", "lg11-lg11", "lg12-lg12", "lgX-lgX",
				"lg1-lg2", "lg1-lg3", "lg1-lg4", "lg1-lg5", "lg1-lg6", "lg1-lg7", "lg1-lg8", "lg1-lg9",
				"lg1-lg10", "lg1-lg11", "lg1-lg12", "lg1-lgX",
				"lg2-lg3", "lg2-lg4", "lg2-lg5", "lg2-lg6", "lg2-lg7", "lg2-lg8", "lg2-lg9", "lg10-lg2",
				"lg11-lg2", "lg12-lg2", "lg2-lgX",
				"lg3-lg4", "lg3-lg5", "lg3-lg6", "lg3-lg7", "lg3-lg8", "lg3-lg9", "lg10-lg3", "lg11-lg3",
				"lg12-lg3", "lg3-lgX",
				"lg4-lg5", "lg4-lg6", "lg4-lg7", "lg4-lg8", "lg4-lg9", "lg10-lg4", "lg11-lg4", "lg12-lg4",
				"lg4-lgX",
				"lg5-lg6", "lg5-lg7", "lg5-lg8", "lg5-lg9", "lg10-lg5", "lg11-lg5", "lg12-lg5", "lg5-lgX",
				"lg6-lg7", "lg6-lg8", "lg6-lg9", "lg10-lg6", "lg11-lg6", "lg12-lg6", "lg6-lgX",
				"lg7-lg8", "lg7-lg9", "lg10-lg7", "lg11-lg7", "lg12-lg7", "lg7-lgX",
				"lg8-lg9", "lg10-lg8", "lg11-lg8", "lg12-lg8", "lg8-lgX",
				"lg10-lg9", "lg11-lg9", "lg12-lg9", "lg9-lgX",
				"lg10-lg11", "lg10-lg12", "lg10-lgX",
				"lg11-lg12", "lg11-lgX",
				"lg12-lgX"))
p2in$inter2 <- factor(p2in$inter, levels = c("lg1-lg1", "lg2-lg2", "lg3-lg3", "lg4-lg4", "lg5-lg5", "lg6-lg6", "lg7-lg7", "lg8-lg8",
				"lg9-lg9", "lg10-lg10", "lg11-lg11", "lg12-lg12", "lgX-lgX",
				"lg1-lg2", "lg1-lg3", "lg1-lg4", "lg1-lg5", "lg1-lg6", "lg1-lg7", "lg1-lg8", "lg1-lg9",
				"lg1-lg10", "lg1-lg11", "lg1-lg12", "lg1-lgX",
				"lg2-lg3", "lg2-lg4", "lg2-lg5", "lg2-lg6", "lg2-lg7", "lg2-lg8", "lg2-lg9", "lg10-lg2",
				"lg11-lg2", "lg12-lg2", "lg2-lgX",
				"lg3-lg4", "lg3-lg5", "lg3-lg6", "lg3-lg7", "lg3-lg8", "lg3-lg9", "lg10-lg3", "lg11-lg3",
				"lg12-lg3", "lg3-lgX",
				"lg4-lg5", "lg4-lg6", "lg4-lg7", "lg4-lg8", "lg4-lg9", "lg10-lg4", "lg11-lg4", "lg12-lg4",
				"lg4-lgX",
				"lg5-lg6", "lg5-lg7", "lg5-lg8", "lg5-lg9", "lg10-lg5", "lg11-lg5", "lg12-lg5", "lg5-lgX",
				"lg6-lg7", "lg6-lg8", "lg6-lg9", "lg10-lg6", "lg11-lg6", "lg12-lg6", "lg6-lgX",
				"lg7-lg8", "lg7-lg9", "lg10-lg7", "lg11-lg7", "lg12-lg7", "lg7-lgX",
				"lg8-lg9", "lg10-lg8", "lg11-lg8", "lg12-lg8", "lg8-lgX",
				"lg10-lg9", "lg11-lg9", "lg12-lg9", "lg9-lgX",
				"lg10-lg11", "lg10-lg12", "lg10-lgX",
				"lg11-lg12", "lg11-lgX",
				"lg12-lgX"))

# eliminate all values of R^2 > 0.9
#sub_p1in <- p1in[!(p1in$R2>=1),]
#sub_p2in <- p1in[!(p2in$R2>=1),]

# make error bars
errbar_lims = group_by(clohi, inter2) %>%
summarize(mean=mean(R2), se=sd(R2)/sqrt(n()),
upper=mean+(2*se), lower=mean-(2*se))

errbar_lims = group_by(p2in, inter2) %>%
summarize(mea.ptn=mean(R2), se=sd(R2)/sqrt(n()),
upper=mean+(2*se), lower=mean-(2*se))


##### final commands to make the main plot: all pairwise comparisons
mean_se_violin =  ggplot() + 
	geom_violin(data=clohi, aes(x= inter2, y = R2, colour=inter2)) +
	geom_point(data=clohi, aes(x= inter2, y=R2), stat="summary", fun.y=mean, fun.ymax=mean, fun.ymin=mean, size=1) +
	geom_errorbar(aes(x=errbar_lims$inter2, ymax=errbar_lims$upper, ymin=errbar_lims$lower), stat='identity', width=.25) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, size=4,  hjust=1), legend.position = "none")

mean_se_violin =  ggplot() + 
	geom_violin(data=p2in, aes(x= inter2, y = R2, colour=inter2)) +
	geom_point(data=p2in, aes(x= inter2, y=R2), stat="summary", fun.y=mean, fun.ymax=mean, fun.ymin=mean, size=1) +
	geom_errorbar(aes(x=errbar_lims$inter2, ymax=errbar_lims$upper, ymin=errbar_lims$lower), stat='identity', width=.25) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, size=4,  hjust=1), legend.position = "none")
	
	
png("clohi.ld_coords.png", width = 480, height = 480)
print(mean_se_violin)
dev.off()

############################
###### Plot LD decay #######
############################

ggplot(p1in, aes(dist, r2)) +
	geom_point() +
		geom_smooth(formula=dist ~ r2, se=F)
## we will estimate a smooth conditional mean
# Loess smoothing is a process by which many statistical softwares do smoothing.
# In ggplot2 this should be done when you have less than 1000 points,
# otherwise it can be time consuming.





# read in the readable file (after sed/bash tweaking)
p1dec <- read.delim("Tge.P1.coord.decay.ld", header=T)
p2dec <- read.delim("Tge.P2.coord.decay.ld", header=T)

p1dec <- read.delim("Tms.P1.coord.decay.ld", header=T)
p2dec <- read.delim("Tms.P2.coord.decay.ld", header=T)

### simulations
# cloning - no recombination
clono <- read.delim("cloning-norecomb.lg.decay.ld")
# cloning - low recombination
clolow <- read.delim("cloning-lowrecomb.lg.decay.ld")
# cloning - high recombination
tce <- read.delim("cloning-highrecomb.lg.decay.ld")

# selfing - no recombination
selno <- read.delim("selfing-norecomb.lg.decay.ld")
# selfing - low recombination
sellow <- read.delim("selfing-lowrecomb.lg.decay.ld")
# selfing - high recombination
selhi <- read.delim("selfing-highrecomb.lg.decay.ld")



### sexual species
tce <- read.delim("Tce_alignedcoords.lg.maf.decay.ld")

## convert several columns as numeric
cols = c(2, 5, 7);    
p1dec[,cols] = apply(p1dec[,cols], 2, function(x) as.numeric(as.character(x)))
p2dec[,cols] = apply(p2dec[,cols], 2, function(x) as.numeric(as.character(x)))

clono[,cols] = apply(clono[,cols], 2, function(x) as.numeric(as.character(x)))
clolow[,cols] = apply(clolow[,cols], 2, function(x) as.numeric(as.character(x)))
clohi[,cols] = apply(clohi[,cols], 2, function(x) as.numeric(as.character(x)))

tce[,cols] = apply(tce[,cols], 2, function(x) as.numeric(as.character(x)))

# make new column : distance between SNPs
p1dec$distance <- p1dec$BP_B - p1dec$BP_A
p2dec$distance <- p2dec$BP_B - p2dec$BP_A

clono$distance <- clono$BP_B - clono$BP_A
clolow$distance <- clolow$BP_B - clolow$BP_A
clohi$distance <- clohi$BP_B - clohi$BP_A

tce$distance <- tce$BP_B - tce$BP_A

## plot
#install.packages("gridExtra")
library("gridExtra")


sub1 <- p1dec[(p1dec$CHR_A=="lg1"),]
plot1<-ggplot(sub1, aes(distance, R2)) +
  geom_point() +
#  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg1")+
  theme(axis.text = element_text(size=6))

sub2 <- p1dec[(p1dec$CHR_A=="lg2"),]
plot2<-ggplot(sub2, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg2")+
  theme(axis.text = element_text(size=6))

sub3 <- p1dec[(p1dec$CHR_A=="lg3"),]
plot3<-ggplot(sub3, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg3")+
  theme(axis.text = element_text(size=6))

sub4 <- p1dec[(p1dec$CHR_A=="lg4"),]
plot4<-ggplot(sub4, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg4")+
  theme(axis.text = element_text(size=6))

sub5 <- p1dec[(p1dec$CHR_A=="lg5"),]
plot5<-ggplot(sub5, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg5")+
  theme(axis.text = element_text(size=6))

sub6 <- p1dec[(p1dec$CHR_A=="lg6"),]
plot6<-ggplot(sub6, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg6")+
  theme(axis.text = element_text(size=6))

sub7 <- p1dec[(p1dec$CHR_A=="lg7"),]
plot7<-ggplot(sub7, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg7")+
  theme(axis.text = element_text(size=6))

sub8 <- p1dec[(p1dec$CHR_A=="lg8"),]
plot8<-ggplot(sub8, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg8")+
  theme(axis.text = element_text(size=6))

sub9 <- p1dec[(p1dec$CHR_A=="lg9"),]
plot9<-ggplot(sub9, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg9")+
  theme(axis.text = element_text(size=6))

sub10 <- p1dec[(p1dec$CHR_A=="lg10"),]
plot10<-ggplot(sub10, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg10")+
  theme(axis.text = element_text(size=6))

sub11 <- p1dec[(p1dec$CHR_A=="lg11"),]
plot11<-ggplot(sub11, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg11")+
  theme(axis.text = element_text(size=6))

sub12 <- p1dec[(p1dec$CHR_A=="lg12"),]
plot12<-ggplot(sub12, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lg12")+
  theme(axis.text = element_text(size=6))

subX <- p1dec[(p1dec$CHR_A=="lgX"),]
plotX<-ggplot(subX, aes(distance, R2)) +
  geom_point() +
  #  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)  +
  ggtitle("lgX") +
  theme(axis.text = element_text(size=6))

pdf("Tms_pop1_LDdecay.pdf", width = 10, paper = 'a4')
png("Tms_pop1_LDdecay.png", width = 480, height = 480)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8,
             plot9, plot10, plot11, plot12, plotX, ncol=4)

dev.off()


ggplot(p2dec, aes(distance, R2)) +
  geom_point() +
  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)

subs <- p2dec[(p2dec$CHR_A=="lg8"),]
subs <- clono[(clono$CHR_A=="lg4"),]
ggplot(subs, aes(distance, R2)) +
  geom_point() +
  facet_grid(CHR_A ~ .) +
  geom_smooth(formula=distance ~ R2, se=F)


############################
##### Plot LD heatmap ######
############################
## for managing the data
library(dplyr)
library(tidyr)
## for the plots and colours
library(ggplot2)
library(hrbrthemes)

# read in the readable file (after sed/bash tweaking)
# asex pops
p1in <- read.delim("Tge_alignedcoords_P1.maf.inter.ld", header=T)
p2in <- read.delim("Tge_alignedcoords_P2.maf.inter.ld", header=T)
p1in <- read.delim("Tms_alignedcoords_P1.maf.inter.ld", header=T)
p2in <- read.delim("Tms_alignedcoords_P2.maf.inter.ld", header=T)

## To calculate mean per lg-lg interaction we will make a new column
# make column with chr-chr
p1in$inter <- paste(p1in$CHR_A,'-',p1in$CHR_B, sep='')
p2in$inter <- paste(p2in$CHR_A,'-',p2in$CHR_B, sep='')

## to calculate the mean per group (lg - lg groups)
p1in_mean <- aggregate(p1in[, 7], list(p1in$inter), mean)
p2in_mean <- aggregate(p2in[, 7], list(p2in$inter), mean)

## set the colnames
colnames(p1in_mean) <- c("inter", "R2")
colnames(p2in_mean) <- c("inter", "R2")


### pop1
dat1 <- data.frame(t(matrix(
                     unlist(strsplit(as.vector(p1in_mean$inter), split = "-")), 
                     ncol = length(p1in_mean$inter), nrow = 2)))
names(dat1) <- c("lga", "lgb")
dat3 <- cbind(dat1, p1in_mean$R2)


### pop2
dat1 <- data.frame(t(matrix(
                     unlist(strsplit(as.vector(p2in_mean$inter), split = "-")), 
                     ncol = length(p2in_mean$inter), nrow = 2)))
names(dat1) <- c("lga", "lgb")
dat3 <- cbind(dat1, p2in_mean$R2)



## common part for Pop1 and Pop2
dat4 <- dat3[,c(2,1,3)]
colnames(dat3) <- c("lga", "lgb", "R2")
colnames(dat4) <- c("lga", "lgb", "R2")
dat5 <- rbind(dat3, dat4)
dat5$lga <- factor(dat5$lga, levels = c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8",
				"lg9", "lg10", "lg11", "lg12", "lgX"))
dat5$lgb <- factor(dat5$lgb, levels = c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8",
				"lg9", "lg10", "lg11", "lg12", "lgX"))

## test before printing
ggplot(data = dat5, aes(x = lga, y = lgb)) +
  geom_tile(aes(fill = R2)) 

png("Tms_pop1_heatmap.png", width = 580, height = 480)
png("Tms_pop2_heatmap.png", width = 580, height = 480)
ggplot(data = dat5, aes(x = lga, y = lgb)) +
  geom_tile(aes(fill = R2)) +
  ggtitle(label = "Tms - pop 2 - LD per lg") +
    #edit legends : guide = guide_legend(reverse = TRUE) will reverse the order in the legend
  scale_fill_distiller(palette="RdPu", trans = "reverse", guide = guide_legend(reverse = TRUE)) +
  theme_ipsum()
dev.off()












# working file: Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/ /+/g' Tms_alignedcoords_P2.maf.inter.ld
# remove the last character
sed -i 's/.$//' Tms_alignedcoords_P2.maf.inter.ld

sed -i 's/+++++++++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld 
sed -i 's/+++++++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld 
sed -i 's/++++++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld

sed -i 's/+++++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/+++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/+++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/++++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld

sed -i 's/+++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/++++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/+++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/++/\t/g' Tms_alignedcoords_P2.maf.inter.ld
sed -i 's/+/\t/g' Tms_alignedcoords_P2.maf.inter.ld

# some '.' were left without tabs between them and the next column. To correct that:
sed -i 's/\.l/\.\tl/g' Tms_alignedcoords_P2.maf.inter.ld

# remove \t in the beginning of the line
sed -i 's/^\t//g' Tms_alignedcoords_P2.maf.inter.ld





# working file: Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/ /+/g' Tms_alignedcoords_P1.maf.inter.ld
# remove the last character
sed -i 's/.$//' Tms_alignedcoords_P1.maf.inter.ld

sed -i 's/+++++++++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/++++++++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld 
sed -i 's/+++++++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld 
sed -i 's/++++++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld

sed -i 's/+++++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/++++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/+++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/++++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/+++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/++++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld

sed -i 's/+++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/++++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/+++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/++/\t/g' Tms_alignedcoords_P1.maf.inter.ld
sed -i 's/+/\t/g' Tms_alignedcoords_P1.maf.inter.ld

# some '.' were left without tabs between them and the next column. To correct that:
sed -i 's/\.l/\.\tl/g' Tms_alignedcoords_P1.maf.inter.ld

# remove \t in the beginning of the line
sed -i 's/^\t//g' Tms_alignedcoords_P1.maf.inter.ld





# make het plots per individual
hetfb <- read.csv2(file = "Tms.fb.het", header=T, sep='\t')
het2 <- read.csv2(file = "Tms.trim2.het", header=T, sep='\t')
het3 <- read.csv2(file = "Tms.trim3.het", header=T, sep='\t')

colnames(hetfb) <- c("INDV", "O.HOM", "E.HOM", "N_SITES", "F", "indvno", "POP")
colnames(het2) <- c("INDV", "O.HOM", "E.HOM", "N_SITES", "F", "indvno", "POP")
colnames(het3) <- c("INDV", "O.HOM", "E.HOM", "N_SITES", "F", "indvno", "POP")
hetfb$ratio <- hetfb$O.HOM/hetfb$N_SITES
het2$ratio <- het2$O.HOM/het2$N_SITES
het3$ratio <- het3$O.HOM/het3$N_SITES


## plot by indv order
hetfb$INDV <- factor(hetfb$INDV, levels = hetfb$INDV)
het2$INDV <- factor(het2$INDV, levels = het2$INDV)
het3$INDV <- factor(het3$INDV, levels = het3$INDV)




##### plot allelic frequency
### all pop
## reading csv with repeated row names: row.names = NULL
Tms <- read.csv2(file = "Tms.coords.allclean.frq", header=T, sep='\t', row.names = NULL)
cols = c(2, 3, 4, 5, 6);    
Tms[,cols] = apply(Tms[,cols], 2, function(x) as.numeric(as.character(x)))
png("Tms_2pops_allelicfreqP.png", width = 480, height = 480)
ggplot(Tms, aes(x = CHROM, y = p, fill = CHROM)) +
  geom_boxplot()
dev.off()

### div pop
## eliminated all positions in each pop with freq = 1 (because pops were called together, several SNPs are fixed between pops)
Tms <- read.csv2(file = "Tms.pop1.coordsclean2.frq", header=T, sep='\t', row.names = NULL)
cols = c(2, 3, 4, 5, 6, 7);    
Tms[,cols] = apply(Tms[,cols], 2, function(x) as.numeric(as.character(x)))
png("Tms_divpop_pop1_allelicfreqP_withoutfreq1.png", width = 480, height = 480)
ggplot(Tms, aes(x = CHROM, y = max, fill = CHROM)) +
  geom_boxplot()
dev.off()



Tms2 <- read.csv2(file = "Tms.pop2.coordsclean2.frq", header=T, sep='\t', row.names = NULL)
cols = c(2, 3, 4, 5, 6, 7);    
Tms2[,cols] = apply(Tms2[,cols], 2, function(x) as.numeric(as.character(x)))
png("Tms_divpop_pop2_allelicfreqP_withoutfreq1.png", width = 480, height = 480)
ggplot(Tms2, aes(x = CHROM, y = max, fill = CHROM)) +
  geom_boxplot()
dev.off()



