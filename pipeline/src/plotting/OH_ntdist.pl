#!c:/perl/bin/perl.exe

# script to summarise nt composition of 5´ and 3´ overhangs (non-templated nucleotides) in bam files from STAR
# requires samtools installed and in $PATH

# 14vi2019

# output (for each 5' only):

# 1. OH<n>_ntcomposition.<file.sam>.txt - distribution of OH length vs read length in each bam file


use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use List::Util qw(sum);
use File::Basename;



#commandline parsing for parameters
GetOptions(
	'infile=s'		=>	\(my $path2in_file),		#path/to/input/file.bam
	'outdir=s'		=>	\(my $path2out_dir),		#path/to/outdir
) or die "Error in command line arguments";


if (! -d $path2out_dir) {
	mkdir($path2out_dir) || die "Cannot create output directory $path2out_dir $!\n";
}


my($infile_name, $dirs, $suffix) = fileparse($path2in_file);
my $datestring = localtime();

print "$datestring: processing $infile_name\n";

my $outfile5="$path2out_dir\/OH5_ntcomposition.$infile_name.txt";


my %OH5_nts_readlen; #hash to include info on nt composition for each read length
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

				#only the reads with an 5' overhang
				if($cigar=~m/^(\d+)S(\d+)M/){
					my $softclip=$1;

					#get the sequence of the OH
					my $readseq=$line[9];
					my $overhang=substr($readseq, 0, $softclip);
					#print "$readseq\t$cigar\t$overhang\n";

					my %counts;
					foreach my $char ( split //, $overhang ) { $counts{$char}++ }
						
					my $countA;
					if(exists ($counts{"A"})){
						$countA = sum( @counts{ qw(A) } );
					}else{$countA=0};

					my $countT;
					if(exists ($counts{"T"})){
						$countT = sum( @counts{ qw(T) } );
					}else{$countT=0};


					my $countG;
					if(exists ($counts{"G"})){
						$countG = sum( @counts{ qw(G) } );
					}else{$countG=0};


					my $countC;
					if(exists ($counts{"C"})){
						$countC = sum( @counts{ qw(C) } );
					}else{$countC=0};
					
					#print "$readseq\t$cigar\t$overhang\t$countA $countC $countT $countG\n";

					undef %counts;

					#collect the info on overhangs in a hash
					if(exists($OH5_nts_readlen{$readlength})){

						my ($OHlenALL,$countAall,$countTall,$countGall,$countCall)=split/:/,$OH5_nts_readlen{$readlength};

						$OHlenALL=$OHlenALL+$softclip;
						$countAall=$countAall+$countA;
						$countTall=$countTall+$countT;
						$countGall=$countGall+$countG;
						$countCall=$countCall+$countC;

						my $oh_stats_aln="$OHlenALL:$countAall:$countTall:$countGall:$countCall";

						$OH5_nts_readlen{$readlength}=$oh_stats_aln;
					}else{
						my $oh_stats_aln="$softclip:$countA:$countT:$countG:$countC";
						$OH5_nts_readlen{$readlength}=$oh_stats_aln;
					}
		
					# add the info on the number of accepted alignments i.e. alns with 5OH overhangs
					if( exists($alignments_length{$readlength}) ){
						++$alignments_length{$readlength};
					}else{
						$alignments_length{$readlength}=1;
					}



				}
			}
		}
	}
}
close(BAM);



open(OH5,">",$outfile5) or die "Cannot open output file $outfile5:  $!\n";
print OH5 "Read_Length\talignments_5oh\tbases_OH\tA\tT\tG\tC\tfrA\tfrT\tfrG\tfrC\n";

foreach my $readlength (sort {$a <=> $b} keys %OH5_nts_readlen){

	my ($OHlenALL,$countAall,$countTall,$countGall,$countCall)=split/:/,$OH5_nts_readlen{$readlength};
	
	my $frA=$countAall/$OHlenALL;
	my $frT=$countTall/$OHlenALL;
	my $frG=$countGall/$OHlenALL;
	my $frC=$countCall/$OHlenALL;
	
	my $all_aln=$alignments_length{$readlength}; #all alignments for given read length

	print OH5 "$readlength\t$all_aln\t$OHlenALL\t$countAall\t$countTall\t$countGall\t$countCall\t$frA\t$frT\t$frG\t$frC\n";

}

exit;
