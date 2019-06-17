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



filenames.bam <- list.files(indir, pattern="*.bam.txt", full.names=TRUE)


for (i in filenames.bam){
	name=basename(i)
	print(name)
	plotfile=paste(name,"OHcomposition.pdf",sep=".")
	plotpth=file.path(plotsdir,plotfile)

	oh_type=sub("_ntcomposition.+.bam.txt","",name, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)

	sample=sub("OH\\d_ntcomposition.","",name, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
	sample=sub(".bam.txt","",sample, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)


	dat=read.table(i,sep="\t",header=T)
	rownames(dat)=dat[,1]

	dat.plot=dat[,c(1,2,8:11)]
	colnames(dat.plot)=c("Read_Length","alignments","A","T","G","C")



	dat.plot.m=melt(dat.plot,id.vars=c("Read_Length","alignments"))
	colnames(dat.plot.m)=c("Read_Length","alignments","Base","Fraction")


	df_labs=as.data.frame(cbind(dat.plot$Read_Length,dat.plot$alignments))
	colnames(df_labs)=c("Read_Length","alignments")


	p1=ggplot(dat.plot.m, aes(x=Read_Length, y=Fraction)) + geom_col(aes(fill=Base))  + ylim(0, 1.3) + 
		scale_fill_manual(values=c("chartreuse4","orange","red","slategrey","blue"),name="Base") + theme_bw() +
  		geom_text(data = df_labs, aes(label = alignments, y=1.04) , vjust = 0.5, hjust= -0.08, angle=90, size=rel(3.5)) +
  		labs(title=paste("Composition of nontemplated read overhangs in file \n", name, sep=" "), 
  			caption="Numbers indicate number of alignments with 5' overhangs for read length",
  			y="Fraction of bases in 5' overhangs",
  			x="Read length") 

  	pdf(plotpth)
  	print(p1)
 	dev.off()

}

