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
	









