#!/usr/bin/env perl

# 2022-04-12

# 1. Function
# Get totoal discordant reads, total clusters, (discordant reads)% in each cluster, (discordant reads)% in fusion cluster, (discordant reads)% support ratio within --flanking_length bp
# Select the minimum sum distance pair as minimum blat distance
# Change chr23 to chrX, and chr24 to chrY
# Drop cluster id

# 2. Input 
# column 1: chromosome of the left segment
# column 2: left segment site
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment site
# column 6: strand of the right segment
# column 7: cluster id
# column 8: split read count (star)
# column 9: read pair count (star)
# column 10: split read count (processed by tophat)
# column 11: potential split read count (processed by tophat)
# column 12: read pair count (processed by tophat)
# column 13: split read count (blat tophat split reads and tophat potential split reads)
# column 14: split read count (blat tophat split reads)
# column 15: minimum distance of read pair to left breakpoint 
# column 16: minimum distance of read pair to right breakpoint 
# column 17: identity
# column 18: minimum blat distace of left breakpoint (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 19: minimum blat distace of right breakpoint (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 20: minimum blat distace of left breakpoint (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 21: minimum blat distace of right breakpoint (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 22: is canonical split site
# column 23: split reads (star)
# column 24: read pairs (star)
# column 25: split reads (processed by tophat)
# column 26: potential split reads (processed by tophat)
# column 27: read pairs (processed by tophat)
# column 28: split reads (blat tophat split read and tophat potential split read)
# column 29: distance of read pair to left breakpoint
# column 30: distance of read pair to right breakpoint
# column 31: distace of split read to left breakpoint when blating to artifact reference (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 32: distace of split read to right breakpoint when blating to artifact reference (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 33: distace of split read to left breakpoint when blating to artifact reference (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 34: distace of split read to right breakpoint when blating to artifact reference (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 35: alignment of left segment
# column 36: alignment of right segment

# 3. Output
# column 1: chromosome of the left segment
# column 2: left segment site
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment site
# column 6: strand of the right segment
# column 7: split read count (star)
# column 8: read pair count (star)
# column 9: split read count (processed by tophat)
# column 10: potential split read count (processed by tophat)
# column 11: read pair count (processed by tophat)
# column 12: split read count (blat tophat split reads and tophat potential split reads)
# column 13: split read count (blat tophat split reads)
# column 14: minimum distance of read pair to left breakpoint 
# column 15: minimum distance of read pair to right breakpoint 
# column 16: identity
# column 17: minimum blat distace of left breakpoint (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 18: minimum blat distace of right breakpoint (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 19: minimum blat distace of left breakpoint (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 20: minimum blat distace of right breakpoint (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 21: is canonical split site
# column 22: total clusters in cluster (left)
# column 23: total clusters in cluster (right)
# column 24: total clusters in cluster (merge)
# column 25: (discordant read)% support fusion
# column 26: sd of (discordant read)% in each of cluster (left)
# column 27: sd of (discordant read)% in each of cluster (right)
# column 28: split reads (star)
# column 29: read pairs (star)
# column 30: split reads (processed by tophat)
# column 31: potential split reads (processed by tophat)
# column 32: read pairs (processed by tophat)
# column 33: split reads (blat tophat split read and tophat potential split read)
# column 34: distance of read pair to left breakpoint
# column 35: distance of read pair to right breakpoint
# column 36: distace of split read to left breakpoint when blating to artifact reference (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 37: distace of split read to right breakpoint when blating to artifact reference (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 38: distace of split read to left breakpoint when blating to artifact reference (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 39: distace of split read to right breakpoint when blating to artifact reference (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 40: alignment of left segment
# column 41: alignment of right segment
# column 42: overlapped cluster ids (left)
# column 43: overlapped cluster ids (right)
# column 44: total reads in each cluster (left)
# column 45: total reads in each cluster (right)


use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $flanking_length=100;

GetOptions(
    'f|flanking_length=i'  => \$flanking_length,
    'h|help'    => sub{usage()}
)||usage();


system("awk -v flanking_length=$flanking_length 'BEGIN{FS=OFS=\"\t\"} NR>1{print \$1,\$2-flanking_length,\$2+flanking_length,\$1,\$2,\$3; print \$4,\$5-flanking_length,\$5+flanking_length,\$4,\$5,\$6;}' $ARGV[0] | awk 'BEGIN{FS=OFS=\"\t\"} {if(\$2<0) \$2=0; print \$0;}' | sort -k1,1 -k2,2n | uniq >temp_fusion_breakpoint.bed");
system("bedtools intersect -wa -wb -a temp_fusion_breakpoint.bed -b temp_cluster.bed | cut -f4-7,9- >temp_overlapped_discordant_reads_in_fusion_breakpoint.tsv");

my %overlapped_discordant_reads;
open IN, "temp_overlapped_discordant_reads_in_fusion_breakpoint.tsv" or die "Can't open temp_overlapped_discordant_reads_in_fusion_breakpoint.tsv:$!";
while(<IN>){
    chomp;
    my ($chr_fusion, $pos_fusion, $strand_fusion, $chr_discordant_read, $pos_discordant_read, $strand_discordant_read, $junc_type_discordant_read,
        $discordant_read_name, $cluster_id_discordant_read)=split "\t",$_;
    my $temp_key=join ":", ($chr_fusion, $pos_fusion, $strand_fusion);
    $overlapped_discordant_reads{$temp_key}{"cluster_reads"}{$cluster_id_discordant_read}{$discordant_read_name}=0;
    next if $strand_fusion ne $strand_discordant_read;
    if($junc_type_discordant_read==1){ # split read, may support fusion, need confirm with another end 
        $overlapped_discordant_reads{$temp_key}{"split_reads"}{$discordant_read_name}=0 if $pos_fusion==$pos_discordant_read;
    }else{ # read pair, may support fusion, need confirm with another end
        $overlapped_discordant_reads{$temp_key}{"read_pairs"}{$discordant_read_name}=0 if (($strand_fusion eq '+') && ($pos_discordant_read<=$pos_fusion)) || (($strand_fusion eq '-') && ($pos_discordant_read>=$pos_fusion));
    }
}
close(IN);

say join "\t", ("Chr_left", "Pos_left", "Strand_left", "Chr_right", "Pos_right", "Strand_right",
    "Split_read_count_(star)", "Read_pair_count_(star)", "Split_read_count_(tophat)", "Potential_split_read_count_(tophat)", "Read_pair_count_(tophat)", "Split_read_count_(blat_tophat_split_and_tophat_potential_split_reads)", "Split_read_count_(blat_tophat_split_reads)",
    "Minimum_read_pair_distance_to_left", "Minimum_read_pair_distance_to_right",
    "Identity", "Minimum_blat_distance_to_left_(tophat_split_and_potential_split_reads)", "Minimum_blat_distance_to_right_(tophat_split_and_potential_split_reads)", "Minimum_blat_distance_to_left_(tophat_split_reads)", "Minimum_blat_distance_to_right_(tophat_split_reads)",
    "Is_canonical_motif",
    "Total_clusters_(left)", "Total_clusters_(right)", "Total_clusters_(merge)",
    "(discordant_reads)%_support_fusion",
    "SD_(discordant_reads)%_in_clusters_(left)", "SD_(discordant_reads)%_in_clusters_(right)",
    "Split_reads_(star)", "Read_pairs_(star)", "Split_reads_(tophat)", "Potential_split_reads_(tophat)", "Read_pairs_(tophat)", "Split_reads_(blat_tophat_split_and_tophat_potential_split_reads)",
    "Read_pair_distance_to_left", "Read_pair_distance_to_right",
    "Blat_distance_to_left_(tophat_split_and_tophat_potential_split_reads)", "Blat_distance_to_right_(tophat_split_and_tophat_potential_split_reads)", "Blat_distance_to_left_(tophat_split_reads)", "Blat_distance_to_right_(tophat_split_eads)",
    "Identity_align_left", "Identity_align_right",
    "Cluster_ids_left", "Cluster_ids_right",
    "Total_reads_in_each_cluster_left", "Total_reads_in_each_cluster_right");

open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $cluster_id,
        $split_read_count_star, $read_pair_count_star, $split_read_count_tophat, $potential_split_read_count_tophat, $read_pair_count_tophat, $split_read_count_blat, $split_read_count_blat_base_on_tophat_split_read,
        $min_read_pair_distance_left, $min_read_pair_distance_right, $identity_output,
	$minimum_blat_distance_left, $minimum_blat_distance_right, $minimum_blat_distance_base_on_tophat_split_read_left, $minimum_blat_distance_base_on_tophat_split_read_right,
        $is_split_site_contain_canonical_motif,
	$split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat, $split_reads_blat,
        $read_pair_distance_to_left, $read_pair_distance_to_right, 
	$blat_distance_left, $blat_distance_right, $blat_distance_base_on_tophat_split_read_left, $blat_distance_base_on_tophat_split_read_right,
	$identity_seq_left, $identity_seq_right)=split "\t", $_;
    my ($total_clusters_left, $percentage_discordant_read_in_each_cluster_left, $all_discordant_reads_left, $clusters_left, $cluster_ids_left, $total_reads_in_each_cluster_left)=&get_cluster_and_discordant_statistics($chr_left, $pos_left, $strand_left, $cluster_id);
    my ($total_clusters_right, $percentage_discordant_read_in_each_cluster_right, $all_discordant_reads_right, $clusters_right, $cluster_ids_right, $total_reads_in_each_cluster_right)=&get_cluster_and_discordant_statistics($chr_right, $pos_right, $strand_right, $cluster_id);
    
    my $key_left=join ":", ($chr_left, $pos_left, $strand_left);
    my $key_right=join ":", ($chr_right, $pos_right, $strand_right);
    my %processed_reads;
    my ($total_discordant_reads, $fusion_supported_reads)=(0) x 2;
    foreach my $read (split ",", $all_discordant_reads_left){
        next if exists $processed_reads{$read};
        $processed_reads{$read}=0;
        $total_discordant_reads++;
        $fusion_supported_reads++ if (exists $overlapped_discordant_reads{$key_left}{"split_reads"}{$read} && exists $overlapped_discordant_reads{$key_right}{"split_reads"}{$read}) ||
            (exists $overlapped_discordant_reads{$key_left}{"read_pairs"}{$read} && exists $overlapped_discordant_reads{$key_right}{"read_pairs"}{$read});
    }
    foreach my $read (split ",", $all_discordant_reads_right){
        next if exists $processed_reads{$read};
        $processed_reads{$read}=0;
        $total_discordant_reads++;
    }
    my $percentage_support_fusion=sprintf("%.2f", $fusion_supported_reads/$total_discordant_reads);
    
    my %prcessed_cluster_ids;
    my $total_clusters=0;
    foreach my $cluster_id(split ",", $clusters_left){
        next if exists $prcessed_cluster_ids{$cluster_id};
        $prcessed_cluster_ids{$cluster_id}=0;
        $total_clusters++;
    }
    foreach my $cluster_id(split ",", $clusters_right){
        next if exists $prcessed_cluster_ids{$cluster_id};
        $prcessed_cluster_ids{$cluster_id}=0;
        $total_clusters++;
    }
    
    # sd
    my @Temp_array=split ",", $percentage_discordant_read_in_each_cluster_left;
    my $sd_percentage_discorant_read_in_each_cluster_left="NA";
    if(@Temp_array>1){
        $sd_percentage_discorant_read_in_each_cluster_left=&get_sd(\@Temp_array);
    }
    @Temp_array=split ",", $percentage_discordant_read_in_each_cluster_right;
    my $sd_percentage_discorant_read_in_each_cluster_right="NA";
    if(@Temp_array>1){
        $sd_percentage_discorant_read_in_each_cluster_right=&get_sd(\@Temp_array);
    }
    
    $chr_left=~s/chr23/chrX/;
    $chr_left=~s/chr24/chrY/;
    $chr_right=~s/chr23/chrX/;
    $chr_right=~s/chr24/chrY/;
    say join "\t", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right,
        $split_read_count_star, $read_pair_count_star, $split_read_count_tophat, $potential_split_read_count_tophat, $read_pair_count_tophat, $split_read_count_blat, $split_read_count_blat_base_on_tophat_split_read,
        $min_read_pair_distance_left, $min_read_pair_distance_right, $identity_output,
	$minimum_blat_distance_left, $minimum_blat_distance_right, $minimum_blat_distance_base_on_tophat_split_read_left, $minimum_blat_distance_base_on_tophat_split_read_right,
        $is_split_site_contain_canonical_motif,
	$total_clusters_left, $total_clusters_right, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right,
        $split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat, $split_reads_blat,
        $read_pair_distance_to_left, $read_pair_distance_to_right, 
	$blat_distance_left, $blat_distance_right, $blat_distance_base_on_tophat_split_read_left, $blat_distance_base_on_tophat_split_read_right,
	$identity_seq_left, $identity_seq_right,
        $cluster_ids_left, $cluster_ids_right,
        $total_reads_in_each_cluster_left, $total_reads_in_each_cluster_right);
}
close(IN);

system("rm -f temp_fusion_breakpoint.bed temp_overlapped_discordant_reads_in_fusion_breakpoint.tsv");

sub get_cluster_and_discordant_statistics{
    my ($chr, $pos, $strand, $cluster_id)=@_;
    my $temp_key=join ":", ($chr, $pos, $strand);
    my @Clusters=sort keys %{$overlapped_discordant_reads{$temp_key}{"cluster_reads"}};
    my $clusters=join ",", @Clusters;
    my %temp_hash;
    foreach my $cluster_id (@Clusters){
        foreach my $read (keys %{$overlapped_discordant_reads{$temp_key}{"cluster_reads"}{$cluster_id}}){
            $temp_hash{$read}=0;
        }
    }
    
    my $total_clusters=@Clusters;
    my $all_discordant_reads=join ",",(keys %temp_hash);
    my $total_discordant_reads=split ",", $all_discordant_reads;
    my $percentage_discordant_read_in_each_cluster="NA";
    my $cluster_ids="NA";
    my $total_reads_in_each_cluster="NA";
    foreach my $temp_cluster_id (@Clusters){
        my $total_cluster_reads=keys %{$overlapped_discordant_reads{$temp_key}{"cluster_reads"}{$temp_cluster_id}};
        my $temp_percentage=sprintf("%.2f", $total_cluster_reads/$total_discordant_reads);
        $percentage_discordant_read_in_each_cluster.=",".$temp_percentage;
        $cluster_ids.=",".$temp_cluster_id;
        $total_reads_in_each_cluster.=",".$total_cluster_reads;
    }
    
    $percentage_discordant_read_in_each_cluster=~s/NA,//;
    $cluster_ids=~s/NA,//;
    $total_reads_in_each_cluster=~s/NA,//;
    return ($total_clusters, $percentage_discordant_read_in_each_cluster, $all_discordant_reads, $clusters, $cluster_ids, $total_reads_in_each_cluster);
}

sub get_sd{
    my($input_array) = @_;
    my @numbers=@$input_array;
    #Prevent division by 0 error in case you get junk data
    return undef unless(scalar(@numbers));

    # Step 1, find the mean of the numbers
    my $total1 = 0;
    foreach my $num (@numbers) {
        $total1 += $num;
    }
    my $mean1 = $total1 / (scalar @numbers);

    # Step 2, find the mean of the squares of the differences
    # between each number and the mean
    my $total2 = 0;
    foreach my $num (@numbers) {
        $total2 += ($mean1-$num)**2;
    }
    my $mean2 = $total2 / (scalar @numbers);

    # Step 3, standard deviation is the square root of the
    # above mean
    my $std_dev = sprintf("%.2f", sqrt($mean2));
    
    return $std_dev;
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
This script was used to get cluster statistics in the breakpoint
Usage: perl $scriptName processed_with_blat.tsv >output
Options:

    -f --flanking_length    calculate cluster statitics within this distance [default: $flanking_length]
    -h --help     	    print this help information
HELP
    exit(-1);
}