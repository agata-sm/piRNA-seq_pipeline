#!c:/perl/bin/perl.exe


# script to produce gff file of cluster coordinates from aggregate cluster detection metrics from proTRAC into "db" file
# input are: (i) output of sumClusters_proTRAC.pl 
# #clusterID  chr start end samples norm_hits_kb  norm_hits_1T  norm_hits_10A norm_hits_24_35 directionality
# and (ii) metadata 
# metadata in `/Users/agata.smialowska/NBISproj/4067_sRNAseq_piRNAs/4067_analysis/pipeline_tsts_rackham/processed/sample_data.txt`
# header in metadata must start with a `#` and contain fields: sample condition replicate; other fields optional
#sample  condition tissue  replicate

# comment on command line parameters
# n is a minimal number of samples (which satisfy this keyword requirement)
# --clusters keyword:value:n for sample annotation
# the keywords and factor levels MUST be identical to column names and factor levels in the sample metadata file
# --clusters cluster:value:n for cluster metrics
# the values must be identical to the following (header of the input file)
# norm_hits_kb	norm_hits_1T	norm_hits_10A	norm_hits_24_35
# directionality is not supported at the moment

# the average of the cluster metrics is calculated for samples satisfying selection criteria (if any)

# OBS! the resulting gff needs to be resorted ? (does not seem so)

# --ox option: --ox factor_in_samples.txt:value for example: --ox condition:ox
# using it will only pull clusters with `ox` samples assigned to them (any number which satisfies the criteria dictated by --samples)




use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(sum);
use List::MoreUtils('indexes','uniq','first_index');


my $script_name="getClusters_gff.0.5.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/infile output by sumClusters_proTRAC.pl\n";
	print "--meta: /path/to/metadata.txt tab-delimted file with sample info\n";	
	print "--outfile: /path/to/outfile.gff\n";
	print "--clusters value:n where value is the name of the metric of interest and n is its average value to satisfy the requirement\n";
	print "--samples keyword:value:n where keyword is column name in metadata.txt, value is its desired value and n is minimal number of samples to satisfy the requirement\n";
	print "--ox factor_name:value where factor_name is column name in metadata.txt that lists status for ox, value is the value of this factor which corresponds to ox samples\n";
	print "--ox setting this option will extract clusters only if they are present in oxidised samples, regardless of the samples filter set in --samples\n";

}


else{
	my $parameters=join(' ', @ARGV);
	
	GetOptions(
		'infile=s' => \(my $infile),
		'meta=s'		=>	\(my $metadata),
		'outfile=s'	=>	\(my $outfile),
		'clusters=s' => \(my $clusters),
		'samples=s' => \(my $samples),
		'ox=s' => \(my $ox)

	) or die "Error in command line arguments";


	my %sample_metadata;

	open(INFILE_META,"<",$metadata) or die "Cannot open input file $metadata : $!";
	
	my @header;
	my $idx_rep;
	my $idx_tis;
	my $idx_cond;
	my %factors;

	while(<INFILE_META>){
		chomp $_;
		my @columns=split/\t/,$_;
		if ($_=~m/^#/)	{
			$columns[0]=~s/#//;
			@header=@columns;

			foreach my $factor (@header){
				my $idx_factor= first_index { lc $_ eq lc $factor } @header;
				$factors{$idx_factor}=$factor;
			}

		}else{
			my $idx_name= first_index { lc $_ eq 'sample' } @header;
			my $sample_name=$columns[$idx_name];

			for my $i (0 .. $#columns) {
      			my $value=$columns[$i];
      			my $name=$factors{$i};

      			$sample_metadata{$sample_name}{$name}=$value;

      			print "$sample_name\t$name\t$value\n";
			}
			
		}
	}
	close(INFILE_META);
			
	#print Dumper(%sample_metadata);


	sub mean {
   		return sum(@_)/@_;
	}


	my $ox_factor_name;
	my $ox_value;
	if (defined $ox){
		($ox_factor_name,$ox_value)=split/:/,$ox;
	}


	print "\n\nprocessing clusters\n\n";

	open(INFILE_CLUST, "<",$infile) or die "Cannot open infile $infile: $!";

	my @column;
	my @cluster_data;

	while(<INFILE_CLUST>){
		chomp $_;
		
		if ($_=~m/^#/){

			@column = split /\t/;
			$column[0]=~s/#//;
		}

		#read in the file into array of hashes
		unless ($_=~m/^#/){

			my @row_line=split /\t/;			

			#select samples
			my @sample_ids=split/,/,$row_line[4];
			print "@sample_ids\n";


			
			#arrays to be subset
			my @clusters_info=split/,/,$row_line[5];
			my @directionalities=split/,/,$row_line[10];


			my @sample_filter_values;
			my @accepted_sample_ids;
			my @ox_status;
			
			if (defined $samples){
				my ($smpl_factor,$smpl_factor_value,$smpl_number)=split/:/,$samples;
				print "$smpl_factor,$smpl_factor_value,$smpl_number\n";

				foreach my $sample_id (@sample_ids){
					print "$sample_id\n";
					my %sample_data=%{$sample_metadata{$sample_id}};
		
					if (exists $sample_data{$smpl_factor}){

						my $sample_factor_value= $sample_data{$smpl_factor};

						if($sample_factor_value eq $smpl_factor_value){

							push @accepted_sample_ids, $sample_id; #better than storing the factor value
						}
						if(defined $ox){
							my $ox_status_sample=$sample_data{$ox_factor_name};
							push @ox_status,$ox_status_sample;
						}			
					}else{
						die "Unknown keyword $smpl_factor\n";
					}
				}
			
				### line selection
				my @uniq_accepted_sample_ids=uniq(@accepted_sample_ids);
				my $selected_sample_number=scalar(@uniq_accepted_sample_ids);

				if ($selected_sample_number >= $smpl_number){

			 		#### sample filtering
					#get the indices of the samples satisfying the $smpl_factor_value
					my @sample_indices;

					for my $accepted_sample (@uniq_accepted_sample_ids){
						my @sample_indices_smpl = indexes { $_ eq $accepted_sample } @sample_ids;
						@sample_indices=(@sample_indices,@sample_indices_smpl);			
					}

					#calculate averages
					for my $metric_idx (6 .. $#row_line-1) {

						my $metric=$row_line[$metric_idx];

						my @metrics=split/,/,$metric; 
						#subset this array if only some samples are selected based on sample annotation
						@metrics=@metrics[@sample_indices];
						my $avg_metric=mean(@metrics);
						$row_line[$metric_idx]=$avg_metric;
					}



					#subset also directionality and samples and clusters

					my $smpl_subset=join ',',@sample_ids[@sample_indices];
					$row_line[4]=$smpl_subset;

					my $clusters_subset=join ',',@clusters_info[@sample_indices];
					$row_line[5]=$clusters_subset;

					my $direc_subset=join ',',@directionalities[@sample_indices];
					$row_line[10]=$direc_subset;

					#push modfied values to the %row hash
					my %row;
					@row{@column} = @row_line;

					#check if there are ox samples amongst ALL samples
					if(defined $ox){
						my @ox_samples=grep{ $_ eq $ox_value } @ox_status;
						my $number_of_ox_samples=scalar(@ox_samples);
						if($number_of_ox_samples>0){
							push @cluster_data, \%row;
						}
					}else{
						push @cluster_data, \%row;	
					}
				}


			}else{ #option samples not used, i.e. all samples are accepted
				#place entire line in the %

				# check ox status
				if (defined $ox){
					foreach my $sample_id (@sample_ids){
						my %sample_data=%{$sample_metadata{$sample_id}};
						my $ox_status_sample=$sample_data{$ox_factor_name};
						push @ox_status,$ox_status_sample;
					}
				}

				for my $metric_idx (6 .. $#row_line-1) {

					my $metric=$row_line[$metric_idx];

					my @metrics=split/,/,$metric; 
					my $avg_metric=mean(@metrics);
					$row_line[$metric_idx]=$avg_metric;
				}


				my %row;
				@row{@column} = @row_line;

				#check if there are ox samples
				if(defined $ox){
					my @ox_samples=grep{ $_ eq $ox_value } @ox_status;
					my $number_of_ox_samples=scalar(@ox_samples);
					if($number_of_ox_samples>0){
						push @cluster_data, \%row;
					}
				}else{
					push @cluster_data, \%row;	
				}
			}
		}
	}
	close(INFILE_CLUST);
	

	#print Dumper(@cluster_data);

	########### filter based on cluster properties > in selected samples only (see above)!

		my @filtered;

		if (defined $clusters){
			my ($filter_clust,$filter_clust_value)=split/:/,$clusters;

			@filtered = grep { $_->{$filter_clust} >= $filter_clust_value } @cluster_data;

		}else{#clusters option not used, i.e. all selected samples are output
			@filtered=@cluster_data;
		}


	open(OUTFILE,">",$outfile) or die "Cannot open output file $outfile : $!";


	for my $hash_ref (@filtered){

	 	my $chr=$$hash_ref{'chr'};
	 	my $start=$$hash_ref{'start'};
	 	my $end=$$hash_ref{'end'};

	 	my $clust_annot="clusterID: cluster_$$hash_ref{'clusterID'} ;samples: $$hash_ref{'samples'}; clusters:$$hash_ref{'clusters'}; filter: $clusters; avg_norm_hits_kb: $$hash_ref{'norm_hits_kb'}; avg_norm_hits_1T: $$hash_ref{'norm_hits_1T'}; avg_norm_hits_10A: $$hash_ref{'norm_hits_10A'}; avg_norm_hits_24_35: $$hash_ref{'norm_hits_24_35'}; directionality: $$hash_ref{'directionality'}";

	 	my $line_gff="$chr\tselected_cluster\tproTRAC\t$start\t$end\t.\t.\t.\t$clust_annot";
		
	 	if (defined $clusters){
	 		my ($filter_clust,$filter_clust_value)=split/:/,$clusters;
	 		for my $metric (keys %$hash_ref){
	 			my $metric_value=$$hash_ref{$metric};
	 			if ( ($metric eq $filter_clust) && ($metric_value >=$filter_clust_value) ){
					print OUTFILE "$line_gff\n";
				}
			}
	 	}else{ #no cluster filter, all lines are printed
	 		print OUTFILE "$line_gff\n";
	 	}
	}

	close(OUTFILE);
}
exit;
