#!/usr/bin/env perl

# 2023-03-08

# 1. Function
# Get totoal discordant reads, total clusters, (discordant reads)% in each cluster, (discordant reads)% in fusion cluster, (discordant reads)% support ratio within --flanking_length bp
# Drop cluster id

# to do in next version:
# Change chr23 to chrX, and chr24 to chrY

# 2. Input 
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 4: chromosome of breakpoint 2
# column 5: breakpoint 2 position (1-based)
# column 6: breakpoint 2 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 7: junction type: -1=encompassing junction (between the mates), 1=split reads
# column 8: cluster id
# column 9: split read count (e.g., 1)
# column 10: read pair count (e.g., 4)
# column 11: minimal distance of read pair to breakpoint 1 
# column 12: minimal distance of read pair to breakpoint 2
# column 13: split reads (e.g., read_23)
# column 14: read pairs (e.g., read_38,read_62,read_70,read_8)
# column 15: distance of read pair to breakpoint 1 
# column 16: distance of read pair to breakpoint 2

# 3. Output
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 4: chromosome of breakpoint 2
# column 5: breakpoint 2 position (1-based)
# column 6: breakpoint 2 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 7: split read count (e.g., 1)
# column 8: read pair count (e.g., 4)
# column 9: minimal distance of read pair to breakpoint 1 
# column 10: minimal distance of read pair to breakpoint 2
# column 11: standard deviation for candidate fusion clusters around breakpoint 1 (e.g., NA, NA for only one cluster)
# column 12: standard deviation for candidate fusion clusters around breakpoint 2 (e.g., NA, NA for only one cluster)
# column 13: split reads (e.g., read_23)
# column 14: read pairs (e.g., read_38,read_62,read_70,read_8)
# column 15: distance of read pair to breakpoint 1 
# column 16: distance of read pair to breakpoint 2
# column 17: discordant read cluster IDs around breakpoint 1 (e.g., 1)
# column 18: discordant read cluster IDs around breakpoint 2 (e.g., 1)
# column 19: supporting reads count in each cluster around breakpoint 1 (e.g., 4)
# column 20: supporting reads count in each cluster around breakpoint 2 (e.g., 3)


use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $flanking_length=100;
my $sd_cutoff=0.1;


GetOptions(
    'f|flanking_length=i'   => \$flanking_length,
    's|sd_cutoff=f'         => \$sd_cutoff,
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

say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2",
    "Split_read_count", "Read_pair_count", 
    "Minimum_read_pair_distance_to_breakpoint_1", "Minimum_read_pair_distance_to_breakpoint_2",
    "SD_(discordant_reads)%_in_breakpoint_1", "SD_(discordant_reads)%_in_breakpoint_2",
    "Split_reads", "Read_pairs",
    "Read_pair_distance_to_breakpoint_1", "Read_pair_distance_to_breakpoint_2",
    "Cluster_ids_breakpoint_1", "Cluster_ids_breakpoint_2",
    "Read_count_in_each_cluster_breakpoint_1", "Read_count_in_each_cluster_breakpoint_2");

open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $cluster_id,
        $split_read_count, $read_pair_count,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
        $split_reads, $read_pairs,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2)=split "\t", $_;
    my ($percentage_discordant_read_in_each_cluster_breakpoint_1, $cluster_ids_breakpoint_1, $total_reads_in_each_cluster_breakpoint_1)=&get_cluster_and_discordant_statistics($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cluster_id);
    my ($percentage_discordant_read_in_each_cluster_breakpoint_2, $cluster_ids_breakpoint_2, $total_reads_in_each_cluster_breakpoint_2)=&get_cluster_and_discordant_statistics($chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cluster_id);
    
    
    # sd
    my @Temp_array=split ",", $percentage_discordant_read_in_each_cluster_breakpoint_1;
    my $sd_percentage_discorant_read_in_each_cluster_breakpoint_1="NA";
    if(@Temp_array>1){
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1=&get_sd(\@Temp_array);
    }
    @Temp_array=split ",", $percentage_discordant_read_in_each_cluster_breakpoint_2;
    my $sd_percentage_discorant_read_in_each_cluster_breakpoint_2="NA";
    if(@Temp_array>1){
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_2=&get_sd(\@Temp_array);
    }
    
    next if ($sd_percentage_discorant_read_in_each_cluster_breakpoint_1 ne "NA" && $sd_percentage_discorant_read_in_each_cluster_breakpoint_1 < $sd_cutoff) ||
        ($sd_percentage_discorant_read_in_each_cluster_breakpoint_2 ne "NA" && $sd_percentage_discorant_read_in_each_cluster_breakpoint_2 < $sd_cutoff);
    
    #$chr_breakpoint_1=~s/chr23/chrX/;
    #$chr_breakpoint_1=~s/chr24/chrY/;
    #$chr_breakpoint_2=~s/chr23/chrX/;
    #$chr_breakpoint_2=~s/chr24/chrY/;
    say join "\t", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2,
        $split_read_count, $read_pair_count,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
        $split_reads, $read_pairs,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2,
        $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
        $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2);
}
close(IN);

system("rm -f temp_fusion_breakpoint.bed temp_overlapped_discordant_reads_in_fusion_breakpoint.tsv");

sub get_cluster_and_discordant_statistics{
    my ($chr, $pos, $strand, $cluster_id)=@_;
    my $temp_key=join ":", ($chr, $pos, $strand);
    my @Clusters=sort keys %{$overlapped_discordant_reads{$temp_key}{"cluster_reads"}};
    my %temp_hash;
    foreach my $cluster_id (@Clusters){
        foreach my $read (keys %{$overlapped_discordant_reads{$temp_key}{"cluster_reads"}{$cluster_id}}){
            $temp_hash{$read}=0;
        }
    }
    
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
    return ($percentage_discordant_read_in_each_cluster, $cluster_ids, $total_reads_in_each_cluster);
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

    -f --flanking_length    INT     Breakpoint flanking size to calculate standard deviation [default: $flanking_length]
    -s --sd_cutoff          FLOAT   Standard deviation cutoff to filter fusions  [default: $sd_cutoff]
    -h --help     	                 Print this help information
HELP
    exit(-1);
}