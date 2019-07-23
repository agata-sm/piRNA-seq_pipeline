#!c:/perl/bin/perl.exe

# script to parse output of proTRAC results.table into 5 field bed files


# 19vi2019 - 25vi2019
# based on older scripts for the same purpose

# outputs
# 5-filed bed file for merging
# summary file with table info from results.table from each cluster
# summary file with library-wide stats



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use List::Util qw(sum);
use File::Basename;
#use Statistics::Basic qw(:all);


#commandline parsing for parameters
GetOptions(
	'infile=s'		=>	\(my $path2in_file),		#path/to/input/file.bam
	'outdir=s'		=>	\(my $path2out_dir),		#path/to/outdir
) or die "Error in command line arguments";


if (! -d $path2out_dir) {
	mkdir($path2out_dir) || die "Cannot create output directory $path2out_dir $!\n";
}

# my $path2out_dir_bed="$path2out_dir\/bed";
# if (! -d $path2out_dir_bed) {
# 	mkdir($path2out_dir_bed) || die "Cannot create output directory $path2out_dir_bed $!\n";
# }

# my $path2out_dir_tab="$path2out_dir\/tab";
# if (! -d $path2out_dir_tab) {
# 	mkdir($path2out_dir_tab) || die "Cannot create output directory $path2out_dir_tab $!\n";
# }

# my $path2out_dir_summary="$path2out_dir\/summary";
# if (! -d $path2out_dir_summary) {
# 	mkdir($path2out_dir_summary) || die "Cannot create output directory $path2out_dir_summary $!\n";
# }

my($infile_name, $dirs, $suffix) = fileparse($path2in_file);
my $datestring = localtime();

print "$datestring: processing $infile_name\n";

$infile_name=~m/results.table.(\S+)$/;

my $smpl_id=$1;
my $outfile_bed="$path2out_dir\/$smpl_id.proTRAC.bed";
my $outfile_tab="$path2out_dir\/$smpl_id.proTRAC.tab";
my $outfile_sum="$path2out_dir\/$smpl_id.proTRAC.summary";




open(INPUT,"<","$path2in_file") or die "Cannot open input file $path2in_file $!";


open(OUTPUT_BED,">", "$outfile_bed") or die "Cannot open input file $outfile_bed $!";
open(OUTPUT_TAB,">", "$outfile_tab") or die "Cannot open input file $outfile_tab $!";
open(OUTPUT_SUM,">", "$outfile_sum") or die "Cannot open input file $outfile_sum $!";

my $header_tab="cluster_id\tcontig\tstart\tend\tsize_bp\thits\thits_normalised\thits_norm_kb\tnorm_hits_1T\tnorm_hits_10A\tnorm_hits_24_35nt\tnorm_hits_main_strand\tdirectionality\tbinding_sites";

print OUTPUT_TAB "$header_tab\n";


my @percnt_1T;
my $clusters_in_library;

while (<INPUT>){
	chomp $_;


	if ($_ =~m/^Cluster.+/){
		my @line=split/\t/,$_;

		$line[0]=~m/Cluster (\d+)$/;
		my $cluster_id="cluster_$1";
		
		$line[1]=~m/Location: (.+)$/;
		my $contig=$1;

		$line[2]=~m/Coordinates: (\d+)-(\d+)/;
		my ($start,$end)=($1,$2);

		my $start_0=$start-1;
		my $line_bed="$contig\t$start_0\t$end\t$cluster_id\t0";
		print OUTPUT_BED "$line_bed\n";

		$line[3]=~m/(\d+)$/;
		my $length=$1;

		my $TFsites;
		my $masked;
		if ($_=~m/Binding sites/){
			if ($line[12] =~m/Binding sites:/){
				$TFsites=$line[12];
				$TFsites=~s/Binding sites: //;
			}else{$TFsites="na";}
		}else{$TFsites="na";}
		
			
		s/.+:\s// for @line;

		my $line_tab1="$cluster_id\t$contig\t$start\t$end\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]";
		my $line_tab2="$line[11]\t$TFsites";
		print OUTPUT_TAB "$line_tab1\t$line_tab2\n";

		my $perc_1T=$line[7];
		$perc_1T=~s/%//;
		#print "$line[7]\t$perc_1T\n";
		push @percnt_1T, $perc_1T;
	}


	#library-wide summary
	if($_ =~m/Total size of (\d+) predicted piRNA clusters: (\d+) bp \((\d+\.?\d*)%\)/ ){
		my ($cluster_number,$total_cluster_size,$cluster_size_perc)=($1,$2,$3);
		$clusters_in_library=$cluster_number;
		print OUTPUT_SUM "sample_id\t$smpl_id\nnumber_of_clusters\t$cluster_number\ntotal_cluster_size_bp\t$total_cluster_size\ncluster_size_perc\t$cluster_size_perc\n";
	}


	if($_ =~m/Non identical sequences that can be assigned to clusters: (\d+) \((\d+\.?\d*)%\)/ ){
		my ($total_uniq_reads,$uniq_reads_perc)=($1,$2);
		print OUTPUT_SUM "unique_tags_in_clusters\t$total_uniq_reads\nperc_unique_tags_in_clusters\t$uniq_reads_perc\n";
	}

	if($_ =~m/Sequence reads that can be assigned to clusters: (\d+) \((\d+\.?\d*)%\)/ ){
		my ($total_reads,$reads_perc)=($1,$2);
		print OUTPUT_SUM "reads_in_clusters\t$total_reads\nperc_reads_in_clusters\t$reads_perc\n";
	}

}
close(INPUT);


my $sum = sum @percnt_1T;
my $avg_perc_1T;

if ($clusters_in_library>0){
	 $avg_perc_1T= $sum /scalar(@percnt_1T);
}else{
	$avg_perc_1T=0;
}


print OUTPUT_SUM "avg_perc_1T_in_clusters\t$avg_perc_1T\n";

exit;


