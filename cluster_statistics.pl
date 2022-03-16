#!/usr/bin/env perl

# 2022-03-15

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
# column 8: split read count (processed by tophat and blat)
# column 9: read pair count (processed by tophat)
# column 10: minimum distance of read pair to left breakpoint 
# column 11: minimum distance of read pair to right breakpoint 
# column 12: identity
# column 13: minimum blat distace of left breakpoint (use mean if both read 1 and read 2 are split read)
# column 14: minimum blat distace of right breakpoint (use mean if both read 1 and read 2 are split read)
# column 15: split reads (processed by tophat and blat)
# column 16: read pairs (processed by tophat)

# 3. Output
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
     "Split_read_count_(Tophat_and_Blat)", "Read_pair_count_(Tophat)",
     "Minimum_read_distance_to_left", "Minimum_read_distance_to_right",
     "Identity", "Minimum_blat_distance_to_left", "Minimum_blat_distance_to_right",
     "Total_clusters_(left)", "Total_clusters_(right)", "Total_clusters_(merge)",
     "(discordant_reads)%_support_fusion",
     "SD_(discordant_reads)%_in_clusters_(left)", "SD_(discordant_reads)%_in_clusters_(right)",
     "Split_reads_(tophat_and_blat)", "Read_pairs_(tophat)");

open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $cluster_id,
        $split_read_count_blat, $read_pair_count_tophat, 
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $identity_output, $minimum_blat_distance_left, $minimum_blat_distance_right,
        $split_reads_blat, $read_pairs_tophat)=split "\t", $_;
    my ($total_clusters_left, $percentage_discordant_read_in_each_cluster_left, $all_discordant_reads_left, $clusters_left)=&get_cluster_and_discordant_statistics($chr_left, $pos_left, $strand_left, $cluster_id);
    my ($total_clusters_right, $percentage_discordant_read_in_each_cluster_right, $all_discordant_reads_right, $clusters_right)=&get_cluster_and_discordant_statistics($chr_right, $pos_right, $strand_right, $cluster_id);
    
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
        $split_read_count_blat, $read_pair_count_tophat,
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $identity_output, $minimum_blat_distance_left, $minimum_blat_distance_right,
        $total_clusters_left, $total_clusters_right, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right,
        $split_reads_blat, $read_pairs_tophat);
}
close(IN);

system("rm temp_fusion_breakpoint.bed temp_overlapped_discordant_reads_in_fusion_breakpoint.tsv");

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
    foreach my $temp_cluster_id (@Clusters){
        my $total_cluster_reads=keys %{$overlapped_discordant_reads{$temp_key}{"cluster_reads"}{$temp_cluster_id}};
        my $temp_percentage=sprintf("%.2f", $total_cluster_reads/$total_discordant_reads);
        $percentage_discordant_read_in_each_cluster.=",".$temp_percentage;
    }
    
    $percentage_discordant_read_in_each_cluster=~s/NA,//;
    return ($total_clusters, $percentage_discordant_read_in_each_cluster, $all_discordant_reads, $clusters);
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