#!/usr/bin/env Rscript

#to be used with putput of 4067_countOH_bam_v3.pl

require(optparse)
require(reshape2)
require(ggplot2)

option_list = list(
	make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="input file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

infile=opt$infile

plots="plots"
plotsdir=file.path(dirname(normalizePath(infile)),plots)
dir.create(plotsdir)


dat=read.table(infile,sep="\t",header=TRUE)

newcolnames=colnames(dat)
newcolnames=sub(".+mapped.bowtie_v.","",newcolnames, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
newcolnames=sub(".bwt.+.bam","",newcolnames, ignore.case = FALSE, perl = FALSE,fixed = FALSE, useBytes = FALSE)
colnames(dat)=newcolnames
rownames(dat)=dat$Geneid
dat=dat[,-1]

dat.f=sweep(dat,2,colSums(dat),"/") 



df_plots=as.data.frame(t(dat.f))


df_plots$sample=rownames(df_plots)
df_to_plot=melt(df_plots,id.vars="sample")
colnames(df_to_plot)=c("sample","biotype","value")



p1=ggplot(df_to_plot, aes(x=sample, y=value)) + geom_col(aes(fill=biotype)) + 
	theme_bw()+
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
	coord_flip()+ 
	scale_fill_manual(values=c("maroon3","chartreuse4","orange","chocolate","slategrey","blue","skyblue","limegreen","red"),name="Biotype")+ 
  	labs(title="Biotype representation", 
  			y="Fraction of all alignments\n(normalised to number of alignments for each read)",
  			x="") 

plotfile=paste("biotype_representation","plot.pdf",sep=".")
plotpth=file.path(plotsdir,plotfile)

pdf(plotpth)
print(p1)
dev.off()
