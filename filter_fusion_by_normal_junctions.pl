#!/usr/bin/env perl

# 2023-03-08

# 1. function
# fusion will be filtered if there are normal junctions within defined window size

# 2. input
# 2.1. discordant read count in normal sample
#column 1: chromosome of breakpoint 1
#column 2: breakpoint 1 position (1-based)
#column 3: breakpoint 1 strand (+ means left of the site will be used, while - means right of site will be used)
#column 4: chromosome of breakpoint 2
#column 5: breakpoint 2 position (1-based)
#column 6: breakpoint 2 strand (+ means left of the site will be used, while - means right of site will be used)
#column 7: read count

# 2.2. filtered fusions
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
# column 13: presence of canonical splice site
# column 14: whether fusion transcript locates in the same gene (e.g., N)
# column 15: fusion type
# column 16: gene name of breakpoint 1
# column 17: gene type of breakpoint 1
# column 18: gene strand of breakpoint 1
# column 19: breakpoint location of breakpoint 1
# column 20: breakpoint region type of breakpoint 1
# column 21: Frame_extra_base of breakpoint 1
# column 22: gene name of breakpoint 2
# column 23: gene type of breakpoint 2
# column 24: gene strand of breakpoint 2
# column 25: breakpoint location of breakpoint 2
# column 26: breakpoint region type of breakpoint 2
# column 27: Frame_extra_base of breakpoint 2
# column 28: Fusion frame (in-frame, out-frame, or unknown)
# column 29: split reads (e.g., read_23)
# column 30: read pairs (e.g., read_38,read_62,read_70,read_8)
# column 31: distance of read pair to breakpoint 1 
# column 32: distance of read pair to breakpoint 2
# column 33: discordant read cluster IDs around breakpoint 1 (e.g., 1)
# column 34: discordant read cluster IDs around breakpoint 2 (e.g., 1)
# column 35: supporting reads count in each cluster around breakpoint 1 (e.g., 4)
# column 36: supporting reads count in each cluster around breakpoint 2 (e.g., 3)


# 3. output
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
# column 13: presence of canonical splice site
# column 14: whether fusion transcript locates in the same gene (e.g., N)
# column 15: fusion type
# column 16: gene name of breakpoint 1
# column 17: gene type of breakpoint 1
# column 18: gene strand of breakpoint 1
# column 19: breakpoint location of breakpoint 1
# column 20: breakpoint region type of breakpoint 1
# column 21: Frame_extra_base of breakpoint 1
# column 22: gene name of breakpoint 2
# column 23: gene type of breakpoint 2
# column 24: gene strand of breakpoint 2
# column 25: breakpoint location of breakpoint 2
# column 26: breakpoint region type of breakpoint 2
# column 27: Frame_extra_base of breakpoint 2
# column 28: Fusion frame (in-frame, out-frame, or unknown)
# column 29: split reads (e.g., read_23)
# column 30: read pairs (e.g., read_38,read_62,read_70,read_8)
# column 31: distance of read pair to breakpoint 1 
# column 32: distance of read pair to breakpoint 2
# column 33: discordant read cluster IDs around breakpoint 1 (e.g., 1)
# column 34: discordant read cluster IDs around breakpoint 2 (e.g., 1)
# column 35: supporting reads count in each cluster around breakpoint 1 (e.g., 4)
# column 36: supporting reads count in each cluster around breakpoint 2 (e.g., 3)



use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $help;
my $adjacent_distance=10000;
my $threshould_normal_read_count=2;
GetOptions(
    'a|adjacent_distance=i'=>	              \$adjacent_distance,
    't|threshould_normal_read_count=i'=>    \$threshould_normal_read_count,
    'h|help'=>        \$help
)||usage(); 
usage () if defined $help;


my %read_count;
if($ARGV[0]=~/gz$/){
    open IN, "zcat $ARGV[0] | " or die "Can't open $ARGV[0]:$!";
}else{
    open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
}
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $read_count)=split "\t", $_;
    $read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}{$pos_left}{$pos_right}=$read_count;
}
close(IN);

say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2", 
    "Split_read_count", "Read_pair_count", 
    "Minimum_read_pair_distance_to_breakpoint_1", "Minimum_read_pair_distance_to_breakpoint_2",
    "SD_(discordant_reads)%_in_breakpoint_1", "SD_(discordant_reads)%_in_breakpoint_2",
    "Have_canonical_motif",
    "Whether_in_the_same_gene",
    "Fusion_type",
    "Gene_name_breakpoint_1", "Gene_type_breakpoint_1", "Gene_strand_breakpoint_1", "Breakpoint_location_breakpoint_1", "Breakpoint_region_type_breakpoint_1", "Frame_extra_base_breakpoint_1", 
				"Gene_name_breakpoint_2", "Gene_type_breakpoint_2", "Gene_strand_breakpoint_2", "Breakpoint_location_breakpoint_2", "Breakpoint_region_type_breakpoint_2", "Frame_extra_base_breakpoint_2", 
    "Fusion_frame",
    "Split_reads", "Read_pairs",
    "Read_pair_distance_to_breakpoint_1", "Read_pair_distance_to_breakpoint_2",
    "Cluster_ids_breakpoint_1", "Cluster_ids_breakpoint_2",
    "Read_count_in_each_cluster_breakpoint_1", "Read_count_in_each_cluster_breakpoint_2");

if($ARGV[1]=~/gz$/){
    open IN, "zcat $ARGV[1] | " or die "Can't open $ARGV[1]:$!";
}else{
    open IN, $ARGV[1] or die "Can't open $ARGV[1]:$!";
}
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, $input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2,
        $split_read_count, $read_pair_count,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
        $is_split_site_contain_canonical_motif,
        $is_in_the_same_gene,
        $fusion_type,
        $gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1,
        $gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2,
        $fusion_frame,
        $split_reads, $read_pairs,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2,
        $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
        $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2)=split "\t", $_;
    
    # format chr_left and chr_right as they may switch after gene annotation
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right)=($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, $input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2);
    if($input_chr_breakpoint_1 eq $input_chr_breakpoint_2){
        ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right)=($input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2, $input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1) if $input_pos_breakpoint_2 < $input_pos_breakpoint_1;
    }else{
        ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right)=($input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2, $input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1) if substr($input_chr_breakpoint_2, 3) < substr($input_chr_breakpoint_1, 3);
    }
    
    
    my $nearby_discordant_read_count_in_normal=0;
    my ($debug_pos_left, $debug_pos_right)=("NA") x 2; # debug
    foreach my $pos_left_normal (sort {$a <=> $b} keys %{$read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}}){
        my $distance_left=abs($pos_left_normal-$pos_left);
        next if $distance_left > $adjacent_distance && $pos_left_normal < $pos_left;
        last if $distance_left > $adjacent_distance && $pos_left_normal > $pos_left;
        foreach my $pos_right_normal (sort {$a <=> $b} keys %{$read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}{$pos_left_normal}}){
            my $distance_right=abs($pos_right_normal-$pos_right);
            next if $distance_right > $adjacent_distance && $pos_right_normal < $pos_right;
            last if $distance_right> $adjacent_distance && $pos_right_normal > $pos_right;
            $nearby_discordant_read_count_in_normal+=$read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}{$pos_left_normal}{$pos_right_normal};
            last if $nearby_discordant_read_count_in_normal >= $threshould_normal_read_count;
            # debug
            #$debug_pos_left.=",".$pos_left_normal;
            #$debug_pos_right.=",".$pos_right_normal;
            #if($nearby_discordant_read_count_in_normal >= $threshould_normal_read_count){
            #    $debug_pos_left=~s/NA,//;
            #    $debug_pos_right=~s/NA,//;
            #    say STDERR join "\t", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right,
            #                            $debug_pos_left, $debug_pos_right);
            #    goto OUTPUT;
            #}
        }
        last if $nearby_discordant_read_count_in_normal >= $threshould_normal_read_count;
    }

    next if $nearby_discordant_read_count_in_normal >= $threshould_normal_read_count;
    
    say join "\t", ($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, $input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2,
        $split_read_count, $read_pair_count,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
        $is_split_site_contain_canonical_motif,
        $is_in_the_same_gene,
        $fusion_type,
        $gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1,
        $gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2,
        $fusion_frame,
        $split_reads, $read_pairs,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2,
        $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
        $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2);
}
close(IN);

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to filter fusions by normal discordant read
Usage: perl $scriptName normal_discordant_read_count.tsv filtered_fusions.tsv >output 

    -a	--adjacent_distance	             INT     Window size to search for breakpoints in normal samples [default: $adjacent_distance]
    -t	--threshould_normal_read_count	  INT     Filter fusions if numbers of reads (discordant pairs or split reads) in normal samples are equal to or more than the specified value [default: $threshould_normal_read_count]
    -h	--help			                                Print this help information screen

HELP
    exit(-1);
}
