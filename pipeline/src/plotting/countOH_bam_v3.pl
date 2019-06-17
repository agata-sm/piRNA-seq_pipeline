#!c:/perl/bin/perl.exe

# script to count the distribution of 5´ and 3´ overhangs (non-templated nucleotides) in bam files from STAR
# requires samtools installed and in $PATH

# v.2. more memory-efficent compared to v.1 > still, it crashes on 3 cores on rackham; as written is not robust to larger data sets
#
# v.3.handles each bam file separately, and the summary across all files is then produced in R



# output (for each 5' and 3' OH separately):

# 1. OH<n>_length.<file.sam>.txt - distribution of OH length vs read length in each bam file


use warnings;
use strict;
use diagnostics;
use Getopt::Long;
#use File::Find;
use File::Basename;
#use File::Path qw(make_path);
#use File::Copy;



#commandline parsing for parameters
GetOptions(
	'infile=s'		=>	\(my $path2in_file),		#path/to/input/file.bam
	'outdir=s'		=>	\(my $path2out_dir),		#path/to/outdir
) or die "Error in command line arguments";



my($infile_name, $dirs, $suffix) = fileparse($path2in_file);
my $datestring = localtime();

print "$datestring: processing $infile_name\n";

my $outfile5="$path2out_dir\/OH5_length.$infile_name.txt";
my $outfile3="$path2out_dir\/OH3_length.$infile_name.txt";


my %OH3_length; #overhang length (arrays) for each alignment length in the file
my %OH5_length;

my %alignments_length; #to count all accepted alignments for each read length


open BAM,"samtools view $path2in_file |" or die "cannot open input file $path2in_file : $!";

while (<BAM>){
	next if(/^(\@)/);
	chomp $_;
	if($_=~m/^(\S+)\t(\d+)\t(\w+)\t(\d+)/){
		my @line=split(/\t/,$_);
		unless($line[5]=~m/N/){ #non-split alignmemts allowed
			my $cigar=$line[5];
			my @numbers=$cigar=~/(\d+)/g; 
			my $readlength=eval join '+', @numbers;
			#my $flag=$line[1];

			if ($readlength<51){#cap the readlength 

				if( exists($alignments_length{$readlength}) ){
					++$alignments_length{$readlength};
				}else{
					$alignments_length{$readlength}=1;
				}

				if($cigar=~m/^(\d+)S(\d+)M/){
					push @{ $OH5_length{$readlength}}, $1;
				}

				if($cigar=~m/(\d+)M(\d+)S$/){
					push @{ $OH3_length{$readlength}}, $2;
				}
			}
		}
	}
}
close(BAM);

	
#resolve the file-specific hashes

open(OH5,">",$outfile5) or die "Cannot open output file $outfile5:  $!\n";
print OH5 "Read_Length\talignments\tOH_1nt\tOH_2nt\tOH_3nt\tOH_4nt\tOH_5nt\tOH_6nt_and_longer\tallOH\n";


foreach my $readlength (sort {$a <=> $b} keys %OH5_length){
	my $totOH_5=scalar @{$OH5_length{$readlength}}; #all alignments with 5 OH for given read length
	my $all_aln=$alignments_length{$readlength}; #all alignments for given read length

	my $OH_1=grep $_ eq qw '1', @{$OH5_length{$readlength}};
	my $OH_2=grep $_ eq qw '2', @{$OH5_length{$readlength}};
	my $OH_3=grep $_ eq qw '3', @{$OH5_length{$readlength}};
	my $OH_4=grep $_ eq qw '4', @{$OH5_length{$readlength}};
	my $OH_5=grep $_ eq qw '5', @{$OH5_length{$readlength}};
	
	my $OH_6m=$totOH_5-$OH_1-$OH_2-$OH_3-$OH_4-$OH_5;

	print OH5 "$readlength\t$all_aln\t$OH_1\t$OH_2\t$OH_3\t$OH_4\t$OH_5\t$OH_6m\t$totOH_5\n";
}
close(OH5);
undef %OH5_length;


open(OH3,">",$outfile3) or die "Cannot open output file $outfile3:  $!\n";
print OH3 "Read_Length\talignments\tOH_1nt\tOH_2nt\tOH_3nt\tOH_4nt\tOH_5nt\tOH_6nt_and_longer\tallOH\n";

foreach my $readlength (sort {$a <=> $b} keys %OH3_length){		
	my $totOH_3=scalar @{$OH3_length{$readlength}};
	my $all_aln=$alignments_length{$readlength};

	my $OH_1=grep $_ eq qw '1', @{$OH3_length{$readlength}};
	my $OH_2=grep $_ eq qw '2', @{$OH3_length{$readlength}};
	my $OH_3=grep $_ eq qw '3', @{$OH3_length{$readlength}};
	my $OH_4=grep $_ eq qw '4', @{$OH3_length{$readlength}};
	my $OH_5=grep $_ eq qw '5', @{$OH3_length{$readlength}};
		
	my $OH_6m=$totOH_3-$OH_1-$OH_2-$OH_3-$OH_4-$OH_5;

	print OH3 "$readlength\t$all_aln\t$OH_1\t$OH_2\t$OH_3\t$OH_4\t$OH_5\t$OH_6m\t$totOH_3\n";
}
close(OH3);

my $datestring2 = localtime();
print "$datestring2: summarising OH distribution done\n";

exit;
