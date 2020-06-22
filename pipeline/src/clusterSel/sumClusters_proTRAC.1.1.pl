#!c:/perl/bin/perl.exe


# script to aggregate cluster detection metrics from proTRAC into "db" file
# input are files  processed by piRNA pipelne *.proTRAC.tab and merged_clusters.bed


##############
### workflow:

# merged clusters from bedops (or similar)
# 	give the unique id > "master list"

# check if given cluster in individual sample file corresponds to a merged cluster from the "master list"
# 	add to master cluster info

##############
### output:

# master_cluster_ID - master_cluster_contig -  master_cluster_start - master_cluster_end - samples_cluster_detected - cluster_info - hits_norm_kb - norm_hits_1T - norm_hits_10A norm_hits_24_35nt - cluster_directionality
# cluster_info for individual clusters
# cluster_id_proTRAC::start::end


##########
### comment:

# due to the fact that sometimes the "master" (merged) cluster encompassed several narrower clusters in individual samples, such sample would then be present multiple times in the output;
# in such cases the information is preserved in a way such that data for each narrow cluster is included
# however, downstream teh information is flattened: each sample ID is present only once and the metrics are averaged
# this problem is best solved by providing "master cluster" bed file that was generated taking fraction of overlap into account rather than merging all clusters in all samples together

##########
### in this version (vs v 0.2):
# multiple --indir paths supported



use warnings;
use strict;
use diagnostics;
use File::Basename;
use Getopt::Long;
use File::Find;


my $script_name="sumClusters_proTRAC.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--indir: /path/to/directory/processed_proTRAC_clusters/*.proTRAC.tab\n";
	print "--merged: /path/to/merged_clusters.bed (bed file; produced by bedops merge or bedtools intersect)\n";	
	print "--outfile: /path/to/outfile (recommended to create new each time new data set is added)\n";
}




else{
	my $parameters=join(' ', @ARGV);
	
	my @indirs;
	GetOptions(
		'indir=s' => \(my $input_dir),
		'merged=s'		=>	\(my $merged_clusters_bed),
		'outfile=s'	=>	\(my $outfile)
	) or die "Error in command line arguments";
#		'indir=s' => \@indirs,


	foreach my $indir_ (@indirs){
		print "TO PROCESS $indir_\n";
	}

	##################
	## parse merged clusters

	my %merged_clusters;
	my $cluster_id=1;

	open(INPUT_MERGED,"<","$merged_clusters_bed") or die "Cannot open input file $merged_clusters_bed $!";
	while(<INPUT_MERGED>){
		chomp $_; 

		my @cluster_coords=split/\t/,$_;
		my $cluster_dta="$cluster_coords[0]::$cluster_coords[1]::$cluster_coords[2]";
		$merged_clusters{$cluster_id}=$cluster_dta;
		++$cluster_id;
	}
	close(INPUT_MERGED);



	##################
	## parse sample cluster data

	#in all the hashes (HoA) below the key is the master cluster ID
	my %clusters_smpls;
	my %clusters_hits_kb;
	my %clusters_nhits_1T;
	my %clusters_nhits_10A;
	my %clusters_nhits_24_35;
	my %clusters_direct;
	my %clusters_info;

#	foreach my $input_dir (@indirs){
		print "process sample data in\n";
		print "$input_dir\n";

		opendir DIR, $input_dir or die "Cannot open directory $input_dir $!";

		my @files=readdir(DIR); #all files & subdirs, local (no paths)
		closedir DIR;

		my @tab_files=grep(/\.tab$/, @files);
		foreach my $filetab (@tab_files){
			$filetab="$input_dir/$filetab";
		}

		#my $number=scalar(@tab_files);
		#print "$number\n@tab_files\n";
		
		foreach my $infile_tab (@tab_files){
			print "processing $infile_tab\n";

			my($infile_name, $dirs, $suffix) = fileparse($infile_tab);
			print "processing file $infile_name\n";

			$infile_name=~m/(\S+)\.proTRAC\.tab/;
			my $sample_id=$1;

			open (INFILE, "<",$infile_tab) or die "Cannot open input file $infile_tab $!";
			while(<INFILE>){
				chomp $_;
				unless ($_=~m/cluster_id/){

					my @line=split/\t/,$_; #cluster_id	contig	start	end	size_bp	hits	hits_normalised	hits_norm_kb	norm_hits_1T	norm_hits_10A	norm_hits_24_35nt	norm_hits_main_strand	directionality	binding_sites

					my $cluster_id_proTRAC=$line[0];
					my $contig=$line[1];
					my $start=$line[2];
					my $end=$line[3];

					my $directionality_long=$line[12];
					$directionality_long=~m/(^\S+)/;# to remove trailing info bi:plus-minus (split between 20214126 and 20214137)
					my $directionality=$1;

					my $nhits_1T =$line[8];
					$nhits_1T=~s/%//;
					$nhits_1T=$nhits_1T/100;

					my $nhits_10A=$line[9];
					$nhits_10A=~s/%//;
					$nhits_10A=$nhits_10A/100;

					my $nhits_24_35=$line[10];
					$nhits_24_35=~s/%//;
					$nhits_24_35=$nhits_24_35/100;

					while (my($m_clusterID,$m_cluster_coords) = each (%merged_clusters) ){

						my($m_cluster_chr,$m_cluster_start,$m_cluster_end)=split/::/,$m_cluster_coords;
							if($contig eq $m_cluster_chr){
								if ( ($start >= $m_cluster_start) && ($end <= $m_cluster_end) ) {

									my $start_gff=$start+1;
									my $cluster_info="$cluster_id_proTRAC::$start_gff::$end";

									push @{$clusters_smpls{$m_clusterID}},$sample_id;
									push @{$clusters_hits_kb{$m_clusterID}},$line[7];
									push @{$clusters_nhits_1T{$m_clusterID}},$nhits_1T;
									push @{$clusters_nhits_10A{$m_clusterID}},$nhits_10A;
									push @{$clusters_nhits_24_35{$m_clusterID}},$nhits_24_35;
									push @{$clusters_direct{$m_clusterID}},$directionality;
									push @{$clusters_info{$m_clusterID}},$cluster_info;
							}
						}
						
					}
				}

			}#INFILE
			close(INFILE);
		}	

#	}#indir


	##### process clusters info collected from all sample.tab files	

	open(OUTFILE,">","$outfile") or die "Cannot open output file $outfile $!";
	my $header="#clusterID\tchr\tstart\tend\tsamples\tclusters\tnorm_hits_kb\tnorm_hits_1T\tnorm_hits_10A\tnorm_hits_24_35\tdirectionality";
	print OUTFILE "$header\n";

	foreach my $master_clusterID (sort { $a <=> $b } keys %merged_clusters){

		my($m_cluster_chr,$m_cluster_start,$m_cluster_end)=split/::/,$merged_clusters{$master_clusterID};

		my @samples=@{$clusters_smpls{$master_clusterID}};

		my $smpl_txt=join ',',@samples;

		my @nhits_kb=@{$clusters_hits_kb{$master_clusterID}};
		my $nhitskb_txt=join ',',@nhits_kb;

		my @nhits_1T=@{$clusters_nhits_1T{$master_clusterID}};
		my $nhits1T_txt=join ',',@nhits_1T;

		my @nhits_10A=@{$clusters_nhits_10A{$master_clusterID}};
		my $nhits10A_txt=join ',',@nhits_10A;

		my @nhits_24_35=@{$clusters_nhits_24_35{$master_clusterID}};
		my $nhits24_35_txt=join ',',@nhits_24_35;

		my @direct=@{$clusters_direct{$master_clusterID}};
		my $direct_txt=join ',',@direct;

		my @clusters_info=@{$clusters_info{$master_clusterID}};
		my $clusters_info_txt=join ',',@clusters_info;

		my $line="$master_clusterID\t$m_cluster_chr\t$m_cluster_start\t$m_cluster_end\t$smpl_txt\t$clusters_info_txt\t$nhitskb_txt\t$nhits1T_txt\t$nhits10A_txt\t$nhits24_35_txt\t$direct_txt";

		print OUTFILE "$line\n";

	}
	close(OUTFILE);
}
exit;




