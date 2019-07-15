#!/usr/bin/env Rscript

#to be used with putput of 4067_countOH_bam_v3.pl

require(optparse)
require(reshape2)
require(ggplot2)

option_list = list(
	make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="directory with input files", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

indir=opt$indir

plots="plots"
plotsdir=file.path(indir,plots)
dir.create(plotsdir)



all_file_stats=list()



filenames <- list.files(indir, pattern="*.summary", full.names=TRUE)


for (i in filenames){
	name=basename(i)
	print(name)

	dat=read.table(i,sep="\t",header=TRUE)
	rownames(dat)=dat[,1]
	#dat=dat[,-1]

	all_file_stats <- c(all_file_stats,list(dat))

}

all_file_stats=as.data.frame(all_file_stats)

all_file_stats=all_file_stats[, -grep("sample_id", colnames(all_file_stats))]


df_plots=all_file_stats[c("number_of_clusters","cluster_size_perc","perc_unique_tags_in_clusters","perc_reads_in_clusters","avg_perc_1T_in_clusters"),]

df_plots=as.data.frame(t(df_plots))
df_plots$sample=rownames(df_plots)


p1=ggplot(df_plots, aes(x=sample,y=number_of_clusters) ) + 
	theme_bw() +
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
	geom_col(color="gray40",fill="gray70") + 
	labs(title="Number of piRNA clusters detected by proTRAC \n", 
  			y="Number of clusters", x="") 


p2=ggplot(df_plots, aes(x=sample,y=cluster_size_perc) ) + 
	theme_bw() +
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
	geom_col(color="gray40",fill="gray70") + 
	labs(title="Percentage of genome comprising piRNA clusters detected by proTRAC \n", 
  			y="Percentage of genome", x="") 

p3=ggplot(df_plots, aes(x=sample,y=perc_unique_tags_in_clusters) ) + 
	theme_bw() +
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
	geom_col(color="gray40",fill="gray70") + 
	labs(title="Percentage of unique sequence tags\nmapped to piRNA clusters detected by proTRAC \n", 
  			y="Percentage of unique tags", x="") 
p4=ggplot(df_plots, aes(x=sample,y=perc_reads_in_clusters) ) + 
	theme_bw() +
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
	geom_col(color="gray40",fill="gray70") + 
	labs(title="Percentage of reads\nmapped to piRNA clusters detected by proTRAC \n", 
  			y="Percentage of reads", x="") 

p5=ggplot(df_plots, aes(x=sample,y=avg_perc_1T_in_clusters) ) + 
	theme_bw() +
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
	geom_col(color="gray40",fill="gray70") + 
	labs(title="Average percentage of 1T - reads\nin reads mapped to piRNA clusters detected by proTRAC \n", 
  			y="Average percentage of reads", x="") 

pdf(paste(plotsdir,"proTRAC_stats_plots.pdf",sep="/"))
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()



