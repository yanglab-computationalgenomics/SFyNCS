#!/usr/bin/env perl

# 2022-03-15

# 1. Function
# Get preliminary fusion statistics, preliminary fusion must meet split_read_tophat_and_blat+read_pair_tophat>=3

# 2. Input
# 2.1. fusion with cluster statistics
# column 1: chromosome of the left segment
# column 2: left segment site
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment site
# column 6: strand of the right segment
# column 7: split read count (processed by tophat and blat)
# column 8: read pair count (processed by tophat)
# column 9: minimum distance of read pair to left breakpoint 
# column 10: minimum distance of read pair to right breakpoint 
# column 11: identity
# column 12: minimum blat distace of left breakpoint (use mean if both read 1 and read 2 are split read)
# column 13: minimum blat distace of right breakpoint (use mean if both read 1 and read 2 are split read)
# column 14: total clusters in cluster (left)
# column 15: total clusters in cluster (right)
# column 16: total clusters in cluster (merge)
# column 17: (discordant read)% support fusion
# column 18: sd of (discordant read)% in each of cluster (left)
# column 19: sd of (discordant read)% in each of cluster (right)
# column 20: split reads (processed by tophat and blat)
# column 21: read pairs (processed by tophat)
# column 22: overlapped cluster ids (left)
# column 23: overlapped cluster ids (right)
# column 24: total reads in each cluster (left)
# column 25: total reads in each cluster (right)

# 2.2. preliminary_candidates.tsv
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: junction type: -1=read pair (between the mates), 1=split reads
# column 8: cluster id
# column 9: split read count
# column 10: read pair count
# column 11: split reads
# column 12: read pairs

# 2.3. processed by tophat
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: cluster id
# column 8: split read count (processed by tophat)
# column 9: read pair count (processed by tophat)
# column 10: minimum distance of read pair to left breakpoint 
# column 11: minimum distance of read pair to right breakpoint 
# column 12: split reads (processed by tophat)
# column 13: read pairs (processed by tophat)
# column 14: distance of read pair to left breakpoint
# column 15: distance of read pair to right breakpoint

# 2.4. processed by blat
# column 1: chromosome of the left segment
# column 2: left segment site
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment site
# column 6: strand of the right segment
# column 7: cluster id
# column 8: split read count (processed by tophat and blat)
# column 9: read pair count (processed by tophat)
# column 10: minimum distance of read pair to left breakpoint 
# column 11: minimum distance of read pair to right breakpoint 
# column 12: identity
# column 13: minimum blat distace of left breakpoint (use mean if both read 1 and read 2 are split read)
# column 14: minimum blat distace of right breakpoint (use mean if both read 1 and read 2 are split read)
# column 15: split reads (processed by tophat and blat)
# column 16: read pairs (processed by tophat)
# column 17: distace of split read to left breakpoint when blating to artifact reference (use mean if both read 1 and read 2 are split read)
# column 18: distace of split read to right breakpoint when blating to artifact reference  (use mean if both read 1 and read 2 are split read)


# 3. Output
# column 1: chromosome of the left segment
# column 2: left segment site
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment site
# column 6: strand of the right segment
# column 7: cluster id
# column 8: split read count (processed by star)
# column 9: read pair count (processed by star)
# column 10: split read count (processed by tophat)
# column 11: read pair count (processed by tophat)
# column 12: split read count (processed by tophat and blat)
# column 13: minimum distance of read pair to left breakpoint 
# column 14: minimum distance of read pair to right breakpoint 
# column 15: identity
# column 16: minimum blat distace of left breakpoint (use mean if both read 1 and read 2 are split read)
# column 17: minimum blat distace of right breakpoint (use mean if both read 1 and read 2 are split read)
# column 18: total clusters in cluster (left)
# column 19: total clusters in cluster (right)
# column 20: total clusters in cluster (merge)
# column 21: (discordant read)% support fusion
# column 22: sd of (discordant read)% in each of cluster (left)
# column 23: sd of (discordant read)% in each of cluster (right)
# column 24: split reads (processed by star)
# column 25: read pairs (processed by star)
# column 26: split reads (processed by tophat)
# column 27: read pairs (processed by tophat)
# column 28: split reads (processed by tophat and blat)
# column 29: overlapped cluster ids (left)
# column 30: overlapped cluster ids (right)
# column 31: total reads in each cluster (left)
# column 32: total reads in each cluster (right)
# column 33: distance of read pair to left breakpoint
# column 34: distance of read pair to right breakpoint
# column 35: distace of split read to left breakpoint when blating to artifact reference (use mean if both read 1 and read 2 are split read)
# column 36: distace of split read to right breakpoint when blating to artifact reference  (use mean if both read 1 and read 2 are split read)


use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

GetOptions(
    'h|help'    => sub{usage()}
)||usage();


my %select_fusions;
open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right,
        $split_read_count_blat, $read_pair_count_tophat,
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $identity_output, $minimum_blat_distance_left, $minimum_blat_distance_right,
        $total_clusters_left, $total_clusters_right, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right,
        $split_reads_blat, $read_pairs_tophat,
        $cluster_ids_left, $cluster_ids_right,
        $total_reads_in_each_cluster_left, $total_reads_in_each_cluster_right)=split "\t", $_;
    $chr_left=~s/chrX/chr23/;
    $chr_left=~s/chrY/chr24/;
    $chr_right=~s/chrX/chr23/;
    $chr_right=~s/chrY/chr24/;
    my $temp_key=join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    $select_fusions{$temp_key}{temp}=0;
}
close(IN);

open IN, $ARGV[1] or die "Can't open $ARGV[1]:$!"; # preliminary_candidates.tsv
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $cluster_id,
        $split_read_count, $read_pair_count, $split_reads, $read_pairs)=split "\t", $_;
    my $temp_key=join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    next if ! exists $select_fusions{$temp_key};
    $select_fusions{$temp_key}{cluster_id}=$cluster_id;
    $select_fusions{$temp_key}{split_read_count_star}=$split_read_count;
    $select_fusions{$temp_key}{read_pair_count_star}=$read_pair_count;
    $select_fusions{$temp_key}{split_reads_star}=$split_reads;
    $select_fusions{$temp_key}{read_pairs_star}=$read_pairs;
}
close(IN);

open IN, $ARGV[2] or die "Can't open $ARGV[2]:$!"; # processed_with_tophat.tsv
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $cluster_id,
        $split_read_count_tophat, $read_pair_count_tophat,
        $min_distance_left, $min_distance_right,
        $split_reads_tophat, $read_pairs_tophat,
        $read_distance_to_left, $read_distance_to_right)=split "\t", $_;
    my $temp_key=join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    next if ! exists $select_fusions{$temp_key};
    $select_fusions{$temp_key}{split_read_count_tophat}=$split_read_count_tophat;
    $select_fusions{$temp_key}{split_reads_tophat}=$split_reads_tophat;
    $select_fusions{$temp_key}{read_pair_distance_to_left}=$read_distance_to_left;
    $select_fusions{$temp_key}{read_pair_distance_to_right}=$read_distance_to_right;
}
close(IN);

open IN, $ARGV[3] or die "Can't open $ARGV[3]:$!"; # processed_with_blat.tsv
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $cluster_id,
        $split_read_count_blat, $read_pair_count_tophat, 
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $identity, $minimum_blat_distance_left, $minimum_blat_distance_right,
        $split_reads_blat, $read_pairs_tophat,
        $blat_distance_left, $blat_distance_right)=split "\t", $_;
    my $temp_key=join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    next if ! exists $select_fusions{$temp_key};
    $select_fusions{$temp_key}{blat_distance_to_left}=$blat_distance_left;
    $select_fusions{$temp_key}{blat_distance_to_right}=$blat_distance_right;
}
close(IN);


say join "\t", ("Chr_left", "Pos_left", "Strand_left", "Chr_right", "Pos_right", "Strand_right", "Cluster_id",
     "Split_read_count_(star)", "Read_pair_count_(star)", "Split_read_count_(tophat)", "Read_pair_count_(tophat)", "Split_read_count_(tophat_and_blat)",
     "Minimum_read_distance_to_left", "Minimum_read_distance_to_right",
     "Identity", "Minimum_blat_distance_to_left", "Minimum_blat_distance_to_right",
     "Total_clusters_(left)", "Total_clusters_(right)", "Total_clusters_(merge)",
     "(discordant_reads)%_support_fusion",
     "SD_(discordant_reads)%_in_clusters_(left)", "SD_(discordant_reads)%_in_clusters_(right)",
     "Split_reads_(star)", "Read_pairs_(star)", "Split_reads_(tophat)", "Read_pairs_(tophat)", "Split_reads_(tophat_and_blat)",
     "Cluster_ids_left", "Cluster_ids_right",
     "Total_reads_in_each_cluster_left", "Total_reads_in_each_cluster_right",
     "Read_pair_distance_to_left", "Read_pair_distance_to_right",
     "Blat_distance_to_left", "Blat_distance_to_right");

open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right,
        $split_read_count_blat, $read_pair_count_tophat,
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $identity_output, $minimum_blat_distance_left, $minimum_blat_distance_right,
        $total_clusters_left, $total_clusters_right, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right,
        $split_reads_blat, $read_pairs_tophat,
        $cluster_ids_left, $cluster_ids_right,
        $total_reads_in_each_cluster_left, $total_reads_in_each_cluster_right)=split "\t", $_;
    my ($temp_chr_left, $temp_chr_right)=($chr_left, $chr_right);
    $temp_chr_left=~s/chrX/chr23/;
    $temp_chr_left=~s/chrY/chr24/;
    $temp_chr_right=~s/chrX/chr23/;
    $temp_chr_right=~s/chrY/chr24/;
    
    my $temp_key=join ",", ($temp_chr_left, $pos_left, $strand_left, $temp_chr_right, $pos_right, $strand_right);
    say join "\t",($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $select_fusions{$temp_key}{cluster_id},
        $select_fusions{$temp_key}{split_read_count_star}, $select_fusions{$temp_key}{read_pair_count_star},
        $select_fusions{$temp_key}{split_read_count_tophat}, $read_pair_count_tophat,
        $split_read_count_blat,
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $identity_output, $minimum_blat_distance_left, $minimum_blat_distance_right,
        $total_clusters_left, $total_clusters_right, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right,
        $select_fusions{$temp_key}{split_reads_star}, $select_fusions{$temp_key}{read_pairs_star},
        $select_fusions{$temp_key}{split_reads_tophat},$read_pairs_tophat,
        $split_reads_blat,
        $cluster_ids_left, $cluster_ids_right,
        $total_reads_in_each_cluster_left, $total_reads_in_each_cluster_right,
        $select_fusions{$temp_key}{read_pair_distance_to_left}, $select_fusions{$temp_key}{read_pair_distance_to_right},
        $select_fusions{$temp_key}{blat_distance_to_left}, $select_fusions{$temp_key}{blat_distance_to_right});
}
close(IN);

sub usage{
    my $scriptName = basename $0;
print <<HELP;
This script was used to get all fusion statistics
Usage: perl $scriptName fusion_statistics.tsv preliminary_candidates.tsv processed_with_tophat.tsv processed_with_blat.tsv >output
Options:

    -h --help     	    print this help information
HELP
    exit(-1);
}