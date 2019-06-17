#!/usr/bin/env Rscript

#to be used with putput of 4067_countOH_bam_v3.pl

require(optparse)
require(reshape2)
require(ggplot2)

option_list = list(
	make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="input directory with bam files", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

indir=opt$indir

plots="plots"
plotsdir=file.path(indir,plots)
dir.create(plotsdir)


filenames <- list.files(indir, pattern="*.tab", full.names=TRUE)


for (i in filenames){
	#get name
	name=basename(i)
	print(name)
	plotfile=paste(name,".plot.pdf",sep="")
	plotpth=file.path(plotsdir,plotfile)

	nt_plot=sub("\\S+.trimmed.","",name, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
	nt_plot=sub("_nt_distr.tab","",nt_plot, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)


	xlabel=paste("Read length",sep=" ")
	title=paste("Base composition at", nt_plot,"nt in file \n",name,sep=" ")

	dat=read.table(i,sep="\t",header=T)
	dat.plot=as.matrix(rbind(dat$A,dat$G,dat$C,dat$T))
	rownames(dat.plot)=c("A","G","C","U")
	colnames(dat.plot)=dat$read_length
	ylimit=as.numeric(max(dat$total))+0.1*(max(dat$total))

	pdf(plotpth)
	barplot(dat.plot, ylab="number of reads", xlab=xlabel, ylim=c(0,ylimit), legend.text=TRUE, col=c("cadetblue2","gold","chartreuse2","mediumpurple3"), main=title)
	dev.off()
} 

