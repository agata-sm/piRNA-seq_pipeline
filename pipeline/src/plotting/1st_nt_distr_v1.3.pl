#!c:/perl/bin/perl.exe


# script to provide a metric for assesment of sRNA seq data
# 1st nt distribution / read length (& number / fraction of reads table)

# 21 viii 2015
# 7 iii 2018

# 17 vi 2019
# use fastq.gz directly
# output all tables to the same directory
# takes infile rather than indir (compatible with a snakemake workflow style)

# usage: perl 1st_nt_distr.pl --indir /path/2/indir --outdir /path/2/outdir

# this has to be run as a slurm job for fastq.gz: the files are copied to SNIC_TMP, ungzipped locally and processed

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Find;
use File::Copy;
use File::Basename;
	
#commandline parsing for parameters
GetOptions(
	'infile=s'		=>	\(my $path2infile),		#path/to/input/fastq_files
	'outdir=s'		=>	\(my $path2out_dir),		#path/to/output/files
) or die "Error in command line arguments";



if (! -d $path2out_dir) {
	mkdir($path2out_dir) || die "Cannot create output directory $path2out_dir $!\n";
}



#######################################################
##### find fastq files
#######################################################

#test whether fastq or fastq.gz	

#######################################################
# regex to find files

my $gz;
if ($path2infile=~m/\S+\.fastq$/){
	$gz="no";
}
elsif ($path2infile=~m/\S+\.fastq\.gz$/){
	$gz="yes";
}


#######################################################
##### misc
#######################################################
#SNIC temp dir env variable
my $snic_tmp='$SNIC_TMP';



#######################################################
##### main
#######################################################



#file names and paths
my($infile_name, $dirs, $suffix) = fileparse($path2infile);
print "processing file $infile_name\n";

$infile_name=~m/(^\S+)\.fastq/;
my $library_name=$1;
print "$library_name\n";

my $out_file_1nt="$library_name\.1_nt_distr.tab";
my $out_file_10nt="$library_name\.10_nt_distr.tab";

my $path2out_file_1="$path2out_dir\/$out_file_1nt";
my $path2out_file_10="$path2out_dir\/$out_file_10nt";


if($gz eq qw /no/){
	open(FH,"<$path2infile") or die "cannot open input file $path2infile : $!";  
}elsif($gz eq qw /yes/){
	open(FH,"gzip -dc < $path2infile |") or die "cannot open input file $path2infile : $!"; 
}

#1st nt info saved in the hash of arrays		
my %nt1_len; #key is read length, value is array of 1st nts    
my %nt10_len; 	#for piRNAs

#read in fastq file in chunks 4 lines each (1 record)
while ((my @lines = map $_ = <FH>, 1..4) [0]) {
	chomp($lines[0],$lines[1],$lines[2],$lines[3]);
	my $sequence=$lines[1];
	my $seq_length=length($sequence);
	if ($seq_length>0){
		my @seq=split //, $sequence;
		my $first_letter=$seq[0];
		my $tenth_letter=$seq[9];

		#print "$sequence\t $first_letter\n"; # it's working; N throws a problem (in Illumina BaseSpace data adapter only reads are marked polyN)

		push @{$nt1_len{$seq_length}}, $first_letter;
		push @{$nt10_len{$seq_length}}, $tenth_letter;
		@lines=undef;
}	}
close (FH);

#process the letters in %nt_len; sort hash keys 
my $header="read_length\ttotal\tA\tG\tC\tT\tfreqA\tfreqG\tfreqC\tfreqT";

open(OUTFILE1,">$path2out_file_1") || die "Couldn't open file $path2out_file_1, $!";
print OUTFILE1 "$header\n";

open(OUTFILE10,">$path2out_file_10") || die "Couldn't open file $path2out_file_10, $!";
print OUTFILE10 "$header\n";

my $total_infile=0;
foreach my $readlength (sort {$a <=> $b} keys %nt1_len){
	my $A_count=grep $_ eq qw 'A', @{$nt1_len{$readlength}};
	my $G_count=grep $_ eq qw 'G', @{$nt1_len{$readlength}};
	my $C_count=grep $_ eq qw 'C', @{$nt1_len{$readlength}};
	my $T_count=grep $_ eq qw 'T', @{$nt1_len{$readlength}};
	my $total=scalar @{$nt1_len{$readlength}};

	my $A_freq=$A_count/$total;
	my $G_freq=$G_count/$total;
	my $C_freq=$C_count/$total;
	my $T_freq=$T_count/$total;

	$total_infile=$total_infile+$total;
	my $line="$readlength\t$total\t$A_count\t$G_count\t$C_count\t$T_count\t$A_freq\t$G_freq\t$C_freq\t$T_freq";
	print OUTFILE1 "$line\n";
}

foreach my $readlength (sort {$a <=> $b} keys %nt10_len){
	my $A_count=grep $_ eq qw 'A', @{$nt10_len{$readlength}};
	my $G_count=grep $_ eq qw 'G', @{$nt10_len{$readlength}};
	my $C_count=grep $_ eq qw 'C', @{$nt10_len{$readlength}};
	my $T_count=grep $_ eq qw 'T', @{$nt10_len{$readlength}};
	my $total=scalar @{$nt10_len{$readlength}};

	my $A_freq=$A_count/$total;
	my $G_freq=$G_count/$total;
	my $C_freq=$C_count/$total;
	my $T_freq=$T_count/$total;

	my $line="$readlength\t$total\t$A_count\t$G_count\t$C_count\t$T_count\t$A_freq\t$G_freq\t$C_freq\t$T_freq";
	print OUTFILE10 "$line\n";
}

close (OUTFILE10);
close (OUTFILE1);
print "total number of reads in file $path2infile : $total_infile\n";



exit;
