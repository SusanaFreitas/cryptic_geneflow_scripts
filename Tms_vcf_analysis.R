#### read and analyse a VCF file with R

library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)

#install.packages("StAMPP", dependencies = TRUE)
#install.packages("gplots", dependencies = TRUE)
#install.packages("ggmap", dependencies = TRUE)
#install.packages("vcfR", dependencies = TRUE)
#install.packages("adegenet", dependencies = TRUE)
#install.packages("adegraphics", dependencies = TRUE)
#install.packages("pegas", dependencies = TRUE)
#install.packages("StAMPP", dependencies = TRUE)
#install.packages("lattice", dependencies = TRUE)
#install.packages("gplots", dependencies = TRUE)
#install.packages("ape", dependencies = TRUE)


#### ADEGENET

### make a PCA plot

Tms <- read.genepop("populations.snps.gen", quiet = TRUE)
# Allele presence absencedata are extracted and NAs replaced usingtab:
X <- tab(Tms, NA.method="mean")
## make PCA
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(Tms))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]
## basic plot
pdf(file="PCA_Tms.pdf")
pdf(file="PCA_Tms_labels.pdf")
plot(pca1$li, col=myCol, cex=3, pch=myPch, xlim =c(-180, 170), ylim=c(-80, 60))
dev.off()

## use wordcloud for non-overlapping labels
library(wordcloud)
textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=0.7, new=FALSE)
dev.off()


### check eigenvalues and other PCA parameters.
## I followed this tutorial: http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization
library(factoextra)
pdf(file="Tms_pca_eig.pdf")
fviz_eig(pca1)
dev.off()


### calculate genetic distances
library(poppr)
## I have to use a genind obj
# Tmsnei <- nei.dist(Tms, warning = TRUE)
## calculate euclidean distance
D <- dist(tab(Tms))
## put them in a tree
library(ape)
tre <- nj(D)
par(xpd=TRUE)
pdf(file="Tms_dist.pdf")
pdf(file="Tms_dist_edgelab.pdf")
plot(tre, type="unrooted", edge.w=2, font =1)
dev.off()
edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))


## PCoA with the distances
pco1 <- dudi.pco(D, scannf=FALSE,nf=2)
s.label(pco1$li*1.1, clabel=0, pch=10)
textplot(pco1$li[,1], pco1$li[,2], words=rownames(pco1$li),cex=0.8, new=FALSE, xpd=TRUE)
title("Principal Coordinate Analysis\n-Tge-")


########################################################################################################################


## read vcf file
vcf3 <- read.vcfR("populations.snps.sed.vcf") #read in all data
vcf2 <- read.vcfR("corrected_output_allxuniq.vcf") #read in all data
head(vcf) #check the vcf object
vcf@fix[1:10,1:5] #check

## read annotation file
gff <- read.table("5_Tge_b3v06.max_arth_b2g_droso_b2g.gff", sep="\t", quote="")

## read sequence file (this is very heavy!!)
dna <- ape::read.dna("5_Tge_b3v08.fasta", format = "fasta")


## plot statistics summed over entire VCF
chrom <- create.chromR(name='Tge_RAD', vcf=vcf3)
plot(chrom) # plot the data

#quick check read depth distribution per individual
dp <- extract.gt(vcf3, element='DP', as.numeric=TRUE)

par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.4, cex.axis=0.5)

###plot read depth distribution per individual
pdf(file="Tms_readdepth.pdf")
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.4, cex.axis=0.5, ylim=c(0,50))
abline(h=8, col="red")
dev.off()


### plot FST pop1 - pop2

library(tidyverse)

# read in data
myData <- read_tsv("./populations.fst_pop1-pop2.tsv")

# clear up
colnames(myData) <- c("locus_ID", "Pop1_ID", "Pop2_ID", "chr", "bp", "col", "overall_pi", "fst", "fishersp", "odds_ratio", "ci_l", "ci_h", "LOD", "corr_amova_fst", "smooth_amova_fst", "smooth_amova_fst_pvalue",  "window_snp_count")


 [1] "locus_ID"                 "Pop1_ID"                  
 [3] "Pop2_ID"                   "chr"                       
 [5] "bp"                         "col"                    
 [7] "overall_pi"                 "fst"                 
 [9] "fishersp"                 "odds_ratio"                
[11] "ci_l"                     "ci_h"                   
[13] "LOD"                        "corr_amova_fst"       
[15] "smooth_amova_fst"         "smooth_amova_fst_pvalue"
[17] "window_snp_count"



# gather fst
fst <- myData %>% select(, AMOVA_fst) %>% gather(key = "stat", value = "fst")
fst$stat <- sub("_fst", "", fst$stat)

# make a plot
a <- ggplot(fst, aes(fst)) + geom_histogram(binwidth = 0.01, fill = "white", colour = "black")
a <- a + theme_light()  + xlab(expression(italic(F[ST])))
a <- a + facet_grid(stat~.)
a

# write it out
#dev.print(pdf, "./fst_plot.pdf", height = 6, width = 8)
ggsave("./fst_plot.pdf", a, width = 6, height = 8)




x <- vcfR2genlight(vcf3)

### give pop names
pop(x) <- as.factor(c((rep("pop1", 23)), (rep("pop2", 23))))


#### ADEGENET
### upload the same libraries from the beginning


## had to change the filename from .genpop to .gen
Tgeind <- read.genepop("populations.snps.gen")
## change the pop names in the genind file (only works with genind file)
pop(Tgeind) <- c((rep("pop1", 23)), (rep("pop2", 23)))
pop(Tgeind)
## convert gening to genpop
Tgepop <- genind2genpop(Tgeind)


# > nLoc(Tgepop)
# [1] 47079

#####################################
### calculate PCA with gening obj ###
#####################################

# Allele presence/absence data are extracted and NAs replaced using tab:
X <- tab(Tgeind, NA.method="mean")

pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(Tgeind))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]

barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
s.class(pca1$li, pop(Tgeind))






## tutorial I am doing:
http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf

#### other websites to follow
https://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html
http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2015/01/krumlov_2015_RAD_tutorial.pdf
http://simison.com/brian/Structure_notes.html
http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf


### ANALYSE LD RESULTS
LD <- read.delim(file = "test.geno.ld", header = TRUE)
ggplot(LD, aes(x = ((LD$POS2 + LD$POS1)/2), y = LD$R^2)) + geom_point()


test.geno.ld
CHR     POS1    POS2    N_INDV  R^2


### Manhattan plot
## taken from here: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
library(lattice)
manhattan.plot<-function(chr, pos, pvalue, 
	sig.level=NA, annotate=NULL, ann.default=list(),
	should.thin=T, thin.pos.places=2, thin.logp.places=2, 
	xlab="Scaffold", ylab="LD: R^2",
	col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {

	if (length(chr)==0) stop("chromosome vector is empty")
	if (length(pos)==0) stop("position vector is empty")
	if (length(pvalue)==0) stop("pvalue vector is empty")

	#make sure we have an ordered factor
	if(!is.ordered(chr)) {
		chr <- ordered(chr)
	} else {
		chr <- chr[,drop=T]
	}

	#make sure positions are in kbp
	if (any(pos>1e6)) pos<-pos/1e6;

	#calculate absolute genomic position
	#from relative chromosomal positions
	posmin <- tapply(pos,chr, min);
	posmax <- tapply(pos,chr, max);
	posshift <- head(c(0,cumsum(posmax)),-1);
	names(posshift) <- levels(chr)
	genpos <- pos + posshift[chr];
	getGenPos<-function(cchr, cpos) {
		p<-posshift[as.character(cchr)]+cpos
		return(p)
	}

	#parse annotations
	grp <- NULL
	ann.settings <- list()
	label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
		col=NULL, fontface=NULL, fontsize=NULL, show=F)
	parse.label<-function(rawval, groupname) {
		r<-list(text=groupname)
		if(is.logical(rawval)) {
			if(!rawval) {r$show <- F}
		} else if (is.character(rawval) || is.expression(rawval)) {
			if(nchar(rawval)>=1) {
				r$text <- rawval
			}
		} else if (is.list(rawval)) {
			r <- modifyList(r, rawval)
		}
		return(r)
	}

	if(!is.null(annotate)) {
		if (is.list(annotate)) {
			grp <- annotate[[1]]
		} else {
			grp <- annotate
		} 
		if (!is.factor(grp)) {
			grp <- factor(grp)
		}
	} else {
		grp <- factor(rep(1, times=length(pvalue)))
	}
  
	ann.settings<-vector("list", length(levels(grp)))
	ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

	if (length(ann.settings)>1) { 
		lcols<-trellis.par.get("superpose.symbol")$col 
		lfills<-trellis.par.get("superpose.symbol")$fill
		for(i in 2:length(levels(grp))) {
			ann.settings[[i]]<-list(pch=pch, 
				col=lcols[(i-2) %% length(lcols) +1 ], 
				fill=lfills[(i-2) %% length(lfills) +1 ], 
				cex=cex, label=label.default);
			ann.settings[[i]]$label$show <- T
		}
		names(ann.settings)<-levels(grp)
	}
	for(i in 1:length(ann.settings)) {
		if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
		ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
			parse.label(ann.settings[[i]]$label, levels(grp)[i]))
	}
	if(is.list(annotate) && length(annotate)>1) {
		user.cols <- 2:length(annotate)
		ann.cols <- c()
		if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
			ann.cols<-match(names(annotate)[-1], names(ann.settings))
		} else {
			ann.cols<-user.cols-1
		}
		for(i in seq_along(user.cols)) {
			if(!is.null(annotate[[user.cols[i]]]$label)) {
				annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
					levels(grp)[ann.cols[i]])
			}
			ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
				annotate[[user.cols[i]]])
		}
	}
 	rm(annotate)

	#reduce number of points plotted
	if(should.thin) {
		thinned <- unique(data.frame(
			logp=round(pvalue,thin.logp.places), 
			pos=round(genpos,thin.pos.places), 
			chr=chr,
			grp=grp)
		)
		logp <- thinned$logp
		genpos <- thinned$pos
		chr <- thinned$chr
		grp <- thinned$grp
		rm(thinned)
	} else {
		logp <- pvalue
	}
	rm(pos, pvalue)
	gc()

	#custom axis to print chromosome names
	axis.chr <- function(side,...) {
		if(side=="bottom") {
			panel.axis(side=side, outside=T,
				at=((posmax+posmin)/2+posshift),
				labels=levels(chr), 
				ticks=F, rot=0,
				check.overlap=F
			)
		} else if (side=="top" || side=="right") {
			panel.axis(side=side, draw.labels=F, ticks=F);
		}
		else {
			axis.default(side=side,...);
		}
	 }

	#make sure the y-lim covers the range (plus a bit more to look nice)
	prepanel.chr<-function(x,y,...) { 
		A<-list();
		maxy<-ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0)))+.5;
		A$ylim=c(0,maxy);
		A;
	}

	xyplot(logp~genpos, chr=chr, groups=grp,
		axis=axis.chr, ann.settings=ann.settings, 
		prepanel=prepanel.chr, scales=list(axs="i"),
		panel=function(x, y, ..., getgenpos) {
			if(!is.na(sig.level)) {
				#add significance line (if requested)
				panel.abline(h=sig.level, lty=2);
			}
			panel.superpose(x, y, ..., getgenpos=getgenpos);
			if(!is.null(panel.extra)) {
				panel.extra(x,y, getgenpos, ...)
			}
		},
		panel.groups = function(x,y,..., subscripts, group.number) {
			A<-list(...)
			#allow for different annotation settings
			gs <- ann.settings[[group.number]]
			A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
			A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
			A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
			A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
			A$x <- x
			A$y <- y
			do.call("panel.xyplot", A)
			#draw labels (if requested)
			if(gs$label$show) {
				gt<-gs$label
				names(gt)[which(names(gt)=="text")]<-"labels"
				gt$show<-NULL
				if(is.character(gt$x) | is.character(gt$y)) {
					peak = which.max(y)
					center = mean(range(x))
					if (is.character(gt$x)) {
						if(gt$x=="peak") {gt$x<-x[peak]}
						if(gt$x=="center") {gt$x<-center}
					}
					if (is.character(gt$y)) {
						if(gt$y=="peak") {gt$y<-y[peak]}
					}
				}
				if(is.list(gt$x)) {
					gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
				}
				do.call("panel.text", gt)
			}
		},
		xlab=xlab, ylab=ylab, 
		panel.extra=panel.extra, getgenpos=getGenPos, ...
	);
}

manhattan.plot(LD$CHR, LD$POS1, LD$R^2)
manhattan.plot(LD$CHR, LD$POS1, LD$R^2, col=c("orange","blue","purple"))

LDscaf <- read.delim(file = "scaf.test.geno.ld", header=TRUE)
manhattan.plot(LDscaf$scaf, LDscaf$POS1, LDscaf$R^2)

### because we are only interested in independent loci (and not dependent loci with R^2 = 1), we will remove all R^2 = 1
LD2 <- LDscaf[ which(LDscaf$R.2<1), ]
manhattan.plot(LD2$scaf, LD2$POS1, LD2$R.2)


## but we dont exactly want the x axis altered.
# Another option is to use this: https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html

LD2$new1 <- paste(LD2$CHR, "-", LD2$POS1, sep="") # new1 corresponds to the POS1
LD2$new2 <- paste(LD2$CHR, "-", LD2$POS2, sep="") # new2 corresponds to the POS2
head(LD2)

LD2$new1 <- as.character(LD2$new1)
ggplot(LD2, aes(R.2, y = new1)) + 
    geom_point(aes(y = y1, col = "y1")) + 
    geom_point(aes(y = y2, col = "y2"))

## alternate between 2 colours
LD2$col = ifelse(as.numeric(LD2$CHR) %% 2 == 0, 0, 1)
## from here: https://digibio.blogspot.com/2018/05/ggplot-and-alternate-shading-with.html
ggplot(LD2, aes(x = new1, y = R.2)) + 
	geom_point(aes(colour = factor(col))) +
	scale_colour_manual(values = c("azure4", "azure3"))


#scale_fill_manual(values = c(rep_len(c("azure4", "azure3"), length(unique(LD2$CHR)))))
#scale_colour_manual(values=rep(c("azure4", "azure3"), ceiling(length(LD2$CHR)/2))[1:length(LD2$CHR)])

## open allele freq file
## because row names are duplicates (same scaffold for several lines) we need to use the "row.names=NULL" option
frq <- read.delim(file = "scaf.test3.freq", header = TRUE, row.names=NULL)
frq$new2 <- paste(frq$CHROM, "-", frq$POS, sep="") # I will call it "new2" because is the equivalent at the POS2 in the LD2 df
ggplot(frq, aes(x = new2, y = freq_p, col = freq_p)) + 
	geom_point()

## The plots look very similar. Will try a correlation between R^2 and allelic frequency
# first I need to merge both dataframes by new
mrg <- merge(LD2,frq,by="new2")





### plot the flagstat values of read alignment before and after filtering
setwd('/home/susana/Dropbox/Timema_cryptic_geneflow/monikensis')
x <- read.table(file='flagstat.total', header=TRUE)


# plot the data
library(ggplot2)
ggplot(x, aes(x = files)) + 
 geom_line(aes(y = total, group = filter, col = filter), linetype = "dashed") +
 geom_line(aes(y = mapped, group = filter, col = filter)) +
 geom_point(aes(y = total, group = filter, col = filter)) +
 geom_point(aes(y = mapped, group = filter, col = filter), shape = 0) +
 theme(axis.text.x = element_text(angle = 90))
 
## plot the melt version of the dataset
library(reshape2)
xmelt <- melt(x, id.vars = c("files", "filter"), measure.vars = c("total", "mapped")

ggplot(xmelt, aes(x = files, y = value, group = filter)) + 
 geom_point(aes(col = filter, shape=variable)) +
 theme(axis.text.x = element_text(angle = 90))

## calculate proportion reads lost/mapped
x$p_read_total <- x$p_read_total
x$p_read_mapped
y <- c()
for (sample in unique(xmelt$files)) {
 line_total <- subset(xmelt, files==sample & filter=="map" & variable=='total')
 dig <- line_total[,4]
 print(dig)
 print(sample)
 line_map <- subset(xmelt, files==sample & filter=="trim3_map" & variable=='total')
 y$name <- sample
 y$p_total <- line_total - line_map
 }

subset(xmelt, filter=="map" & variable=='total')










########################################################################################################
################################# PLOT READ DEPTH AND OTHER VCF STATS ##################################
########################################################################################################


## plot DP for the first vcf produced by freebayes (before any filtering)
vcfpre <- read.vcfR('Tms.fb.vcf')
dp_pre <- extract.gt(vcfpre, element = "DP", as.numeric=TRUE)
dp_pre[1:4,1:6]
## set the plot parameters: margin size
par(mar=c(4,4,4,2))
boxplot(dp_pre, col=2:8, las=3)
pdf(file='Tms.fb.DP.pdf', height = 5, width = 10)
pdf(file='Tms.fb.DP-lowx.pdf', height = 5, width = 10)
boxplot(dp_pre, las=3, col=c("indianred3", "darkgoldenrod1", "cyan3", "hotpink3", "olivedrab3", "orangered3"),
         ylab="Read Depth (DP)", cex=0.6, cex.axis=0.8)
boxplot(dp_pre, las=3, col=c("indianred3", "darkgoldenrod1", "cyan3", "hotpink3", "olivedrab3", "orangered3"),
         ylab="Read Depth (DP)", cex=0.6, cex.axis=0.8, ylim=c(0,400))
## add horizontal lines
abline(h=10, col="red")
abline(h=400, col="tomato3")
abline(h=1000, col="tomato")
abline(h=2000, col="orange2")
## with labels
text(x=0.4, y=10, labels=10, adj=c(0.9,-0.1), col='red', cex=0.6)
text(x=0.4, y=400, labels=400, adj=c(0.9,-0.1), col='tomato3', cex=0.6)
text(x=0.4, y=1000, labels=1000, adj=c(0.9,-0.1), col='tomato', cex=0.6)
text(x=0.4, y=2000, labels=2000, adj=c(0.9,-0.1), col='orange2', cex=0.6)
dev.off()

## colours: colourful
"indianred3", "darkgoldenrod1", "cyan3", "hotpink3", "olivedrab3", "orangered3"
## colours: 2 shades of gray
"#C0C0C0", "#808080"


######## Repeat the same as before for the other vcf files - after filtering
vcf <- read.vcfR('Tms.fb.bial.vcf')
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
pdf(file='Tms.fb.bial.DP.pdf', height = 5, width = 10)
pdf(file='Tms.fb.bial.DP-lowx.pdf', height = 5, width = 10)
boxplot(dp, las=3, col=c("indianred3", "darkgoldenrod1", "cyan3", "hotpink3", "olivedrab3", "orangered3"),
         ylab="Read Depth (DP)", cex=0.6, cex.axis=0.8)
boxplot(dp, las=3, col=c("indianred3", "darkgoldenrod1", "cyan3", "hotpink3", "olivedrab3", "orangered3"),
         ylab="Read Depth (DP)", cex=0.6, cex.axis=0.8, ylim=c(0,400))
## add horizontal lines
abline(h=10, col="red")
abline(h=400, col="tomato3")
abline(h=1000, col="tomato")
abline(h=2000, col="orange2")
## with labels
text(x=0.4, y=10, labels=10, adj=c(0.9,-0.1), col='red', cex=0.6)
text(x=0.4, y=400, labels=400, adj=c(0.9,-0.1), col='tomato3', cex=0.6)
text(x=0.4, y=1000, labels=1000, adj=c(0.9,-0.1), col='tomato', cex=0.6)
text(x=0.4, y=2000, labels=2000, adj=c(0.9,-0.1), col='orange2', cex=0.6)
dev.off()








dp[1:4,1:6]

par(mar=c(12,4,4,2))
boxplot(dp, col=2:8, las=3)
title(ylab = "Depth (DP)")
par(mar=c(5,4,4,2))




































