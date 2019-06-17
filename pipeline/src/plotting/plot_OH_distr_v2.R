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



all_file_stats_3=list()
all_file_stats_5=list()


filenames.bam <- list.files(indir, pattern="*.bam.txt", full.names=TRUE)



for (i in filenames.bam){
	name=basename(i)
	print(name)
	plotfile=paste(name,"plot.pdf",sep=".")
	plotpth=file.path(plotsdir,plotfile)

	oh_type=sub("_length.+.bam.txt","",name, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)

	sample=sub("OH\\d_length.","",name, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
	sample=sub(".bam.txt","",sample, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)


	dat=read.table(i,sep="\t",header=T)
	rownames(dat)=dat[,1]

	dat$f_1nt=dat$OH_1nt/dat$alignments
	dat$f_2nt=dat$OH_2nt/dat$alignments
	dat$f_3nt=dat$OH_3nt/dat$alignments
	dat$f_4nt=dat$OH_4nt/dat$alignments
	dat$f_5nt=dat$OH_5nt/dat$alignments
	dat$f_6nt_and_longer=dat$OH_6nt_and_longer/dat$alignments
	dat$f_allOH=dat$allOH/dat$alignments


	if(oh_type == "OH3"){
		all_file_stats_3[[sample]]=list()
		all_file_stats_3[[sample]][["totaln"]]=sum(dat$alignments)
		all_file_stats_3[[sample]][["oh_cnt"]]=sum(dat$allOH)
		all_file_stats_3[[sample]][["oh_f"]]=sum(dat$allOH)/sum(dat$alignments)
	}

	if(oh_type == "OH5"){
		all_file_stats_5[[sample]]=list()
		all_file_stats_5[[sample]][["totaln"]]=sum(dat$alignments)
		all_file_stats_5[[sample]][["oh_cnt"]]=sum(dat$allOH)
		all_file_stats_5[[sample]][["oh_f"]]=sum(dat$allOH)/sum(dat$alignments)
	}
	

	dat=dat[,-c(3:9)]
	
	df_labs=as.data.frame(cbind(dat$Read_Length,dat$alignments,dat$f_allOH))
	colnames(df_labs)=c("Read_Length","alignments","frac_allOH")
	
	dat=dat[,-9]

	
	dat.plot.m=melt(dat,id.vars=c("Read_Length","alignments"))
	colnames(dat.plot.m)=c("Read_Length","alignments","Overhang_Length","Fraction")
	dat.plot.m$Overhang_Length=sub("f_", "", dat.plot.m$Overhang_Length, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)

	max(dat.plot.m$Fraction)

	ylim_plot=max(dat.plot.m$Fraction)+0.1

	p2=ggplot(dat.plot.m, aes(x=Read_Length, y=Fraction)) + geom_col(aes(fill=Overhang_Length)) + ylim(0,ylim_plot) +
		scale_fill_manual(values=c("chartreuse4","orange","red","slategrey","blue","skyblue"),name="Overhang length")+ theme_bw()+
  		geom_text(data = df_labs, aes(label = alignments, y = frac_allOH), vjust = 0.5, hjust= -0.08, angle=90, size=rel(3.5)) + #Here the text is added #https://stackoverflow.com/questions/50026999/how-to-add-a-single-label-to-a-stacked-bar-chart-in-ggplot2
  		labs(title=paste("Distribution of nontemplated read overhangs in file \n", name, sep=" "), 
  			caption="Numbers indicate total number of alignments for read length",
  			y="Fraction of all alignments",
  			x="Read length") 

  	pdf(plotpth)
  	print(p2)
 	dev.off()
	
}




all_file_stats5=data.frame(matrix(unlist(all_file_stats_5), nrow=length(all_file_stats_5), byrow=T),stringsAsFactors=FALSE)
rownames(all_file_stats5)=names(all_file_stats_5)
colnames(all_file_stats5)=c("totaln","oh_cnt_5","oh_f_5")
all_file_stats5$sample=rownames(all_file_stats5)

all_file_stats3=data.frame(matrix(unlist(all_file_stats_3), nrow=length(all_file_stats_3), byrow=T),stringsAsFactors=FALSE)
rownames(all_file_stats3)=names(all_file_stats_3)
colnames(all_file_stats3)=c("totaln","oh_cnt_3","oh_f_3")
all_file_stats3$sample=rownames(all_file_stats3)


all_file_stats=merge(all_file_stats5,all_file_stats3,by="sample",sort=FALSE)

summary=all_file_stats[,-c(3,5,6)]
colnames(summary)=c("sample","alignments","OH5","OH3")

summary.m=melt(summary,id.vars=c("sample","alignments"))

p3=ggplot(summary.m, aes(x=sample, y=value)) + geom_col(aes(fill=variable), position=position_dodge()) + 
		scale_fill_manual(values=c("chartreuse4","orange"),name="Overhang type")+ theme_bw()+ coord_flip()+ 
		labs(title=paste("Overhang type distribution in all files"), y="Fraction of all alignments", y="File")


plotfile=paste("all_files_OH","plot.pdf",sep=".")
plotpth_all=file.path(plotsdir,plotfile)

pdf(plotpth_all)
print(p3)
dev.off()

