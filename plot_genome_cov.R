### plot_genome_cov.R
### using command line args
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### aborts if no command line args provided
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}



library(ggplot2)


in_file_name = args[1]
#in_file_name = "/Users/dparker/Desktop/Tbi_F_CC87C_to_Tbi_v8_peUP_BWA_mapqfilt_30_coverage_genomecov.txt"


cov_plot <- function(in_file_n, min_cov, max_plot){
	samp_name = strsplit(in_file_n, "/")[[1]][length(strsplit(in_file_name, "/")[[1]])]
	dat1 <- read.table(in_file_name, head = T)
	dat1_filt2  = subset(dat1, dat1$depth >= min_cov)
	#dat1_filt2 = subset(dat1_filt, dat1_filt$depth <= max_cov)	
	d2 <- rep(dat1$depth, dat1$Nsites) 
	q99 <- quantile(d2, .99)	
	
	d2_filt <- rep(dat1_filt2$depth, dat1_filt2$Nsites) 	
	
	med_cov <- median(d2)
	
	med_cov_filt <- median(d2_filt)
	print(med_cov)	
	print(med_cov_filt)
	print(q99)
	
	p1 <- ggplot(dat1 , aes(depth,  Nsites )) + geom_col() + theme_bw() + xlim(c(-3,max_plot)) + 
		ggtitle(samp_name) + 
		geom_vline(xintercept = med_cov, linetype="dotted", color = "red", size=1.5) +
		geom_vline(xintercept = med_cov_filt, linetype="dashed", color = "blue", size=1.5) +    
		geom_vline(xintercept = med_cov * 2, linetype="dotted", color = "orange", size=1.5) +
		geom_vline(xintercept = med_cov_filt * 2, linetype="dashed", color = "green", size=1.5)  +
		geom_vline(xintercept = q99, color = "black", size=1.5)    
	
	p2 <- ggplot(dat1_filt2 , aes(depth,  Nsites )) + geom_col() + theme_bw() + xlim(c(-3,max_plot)) + 
		ggtitle(samp_name) + 
		geom_vline(xintercept = med_cov, linetype="dotted", color = "red", size=1.5) +
		geom_vline(xintercept = med_cov_filt, linetype="dashed", color = "blue", size=1.5) +    
		geom_vline(xintercept = med_cov * 2, linetype="dotted", color = "orange", size=1.5) +
		geom_vline(xintercept = med_cov_filt * 2, linetype="dashed", color = "green", size=1.5)  +
		geom_vline(xintercept = q99, color = "black", size=1.5)    
	

		
	out_df <- as.data.frame(rbind(c(samp_name, med_cov, med_cov_filt, med_cov * 2, med_cov_filt * 2,  q99)))
	colnames(out_df) <- c("sample", "med_cov", "med_filt_cov", "med_cov_x2", "med_filt_cov_x2", "q99")

	output = list("p1" = p1, "p2" = p2, "out_df" = out_df)
	return(output)		
		
			
}

output <- cov_plot(in_file_name, 1,65)

# output$p1
# output$out_df

write.csv(output$out_df, file = paste(in_file_name, "covest.csv", sep = ""), row.names = F)

png(filename = paste(in_file_name, "covest.png", sep = ""), width = 8, height = 6, units = "in", bg = "white", res = 300)
output$p1
dev.off()
getwd() ## where has my plot gone....

png(filename = paste(in_file_name, "covest_2.png", sep = ""), width = 8, height = 6, units = "in", bg = "white", res = 300)
output$p2
dev.off()
getwd() ## where has my plot gone....







