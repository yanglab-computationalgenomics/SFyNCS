#!/usr/bin/env perl

# 2022-11-04

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
# column 3: breakpoint 1 strand 
# column 4: chromosome of breakpoint 2
# column 5: breakpoint 2 position (1-based)
# column 6: breakpoint 2 strand
# column 7: split read count reported by STAR (e.g., 1)
# column 8: read pair count reported by STAR (e.g., 4)
# column 9: split read count supported by TopHat (e.g., 1)
# column 10: potential split read count not supported by TopHat (e.g., 0)
# column 11: read pair count reported by TopHat
# column 12: split read count supported by Blat (e.g., 1)
# column 13: split read count supported by Blat (considering split read supported by TopHat only, e.g., 1)
# column 14: minimal distance between read pair and breakpoint 1 after aligning by TopHat (e.g., 13)
# column 15: minimal distance between read pair and breakpoint 2 after aligning by TopHat (e.g., 22)
# column 16: sequence identity (e.g., 0.52)
# column 17: minimal distance between split read and breakpoint 1 after aligning by Blat (e.g., 0, use mean if both read 1 and read 2 are split read)
# column 18: minimal distance between split read and breakpoint 2 after aligning by Blat (e.g., 0, use mean if both read 1 and read 2 are split read)
# column 19: minimal distance between split read and breakpoint 1 after aligning by Blat (considering split read supported by TopHat only, e.g., 0, use mean if both read 1 and read 2 are split read)
# column 20: minimal distance between split read and breakpoint 2 after aligning by Blat (considering split read supported by TopHat only, e.g., 0, use mean if both read 1 and read 2 are split read)
# column 21: presence of canonical splice site
# column 22: cluster count around breakpoint 1 (e.g., 1)
# column 23: cluster count around breakpoint 2 (e.g., 1)
# column 24: cluster count around breakpoint 1 and breakpoint 2 (e.g., 1)
# column 25: percentage of split reads and read pairs supported fusion transcript’s cluster around breakpoint 1 and breakpoint 2 (e.g., 0.75)
# column 26: standard deviation for candidate fusion clusters around breakpoint 1 (e.g., NA, NA for only one cluster)
# column 27: standard deviation for candidate fusion clusters around breakpoint 2 (e.g., NA, NA for only one cluster)
# column 28: whether fusion transcript locates in the same gene (e.g., N)
# column 29: fusion type
# column 30: gene name of breakpoint 1
# column 31: gene type of breakpoint 1
# column 32: gene strand of breakpoint 1
# column 33: breakpoint location of breakpoint 1
# column 34: breakpoint region type of breakpoint 1
# column 35: Frame_extra_base of breakpoint 2
# column 36: gene name of breakpoint 2
# column 37: gene type of breakpoint 2
# column 38: gene strand of breakpoint 2
# column 39: breakpoint location of breakpoint 2
# column 40: breakpoint region type of breakpoint 2
# column 41: Frame_extra_base of breakpoint 2
# column 42: Fusion frame (in-frame, out-frame, or unknown)
# column 43: split reads reported by STAR (e.g., read_23)
# column 44: read pairs reported by STAR (e.g., read_38,read_62,read_70,read_8)
# column 45: split reads supported by TopHat (e.g., read_23)
# column 46: potential split reads not supported by TopHat (e.g., NA)
# column 47: read pairs supported by TopHat (e.g., read_38,read_70,read_8)
# column 48: split reads supported by Blat(blat tophat split read and tophat potential split read)
# column 49: distance between each read pair and breakpoint 1 after aligning by TopHat (e.g., 740,23,13)
# column 50: distance between each read pair and breakpoint 2 after aligning by TopHat (e.g., 292,2592,22)
# column 51: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 52: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 53: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 54: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 55: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 1 (e.g., --AGTGGGCCAGGTAG-GGCTGG)
# column 56: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 2 (e.g., CCACT--GCCAGG-AGAACCTCA)
# column 57: discordant read cluster IDs around breakpoint 1 (e.g., 1)
# column 58: discordant read cluster IDs around breakpoint 2 (e.g., 1)
# column 59: supporting reads count in each cluster around breakpoint 1 (e.g., 4)
# column 60: supporting reads count in each cluster around breakpoint 2 (e.g., 3)


# 3. output
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand 
# column 4: chromosome of breakpoint 2
# column 5: breakpoint 2 position (1-based)
# column 6: breakpoint 2 strand
# column 7: split read count reported by STAR (e.g., 1)
# column 8: read pair count reported by STAR (e.g., 4)
# column 9: split read count supported by TopHat (e.g., 1)
# column 10: potential split read count not supported by TopHat (e.g., 0)
# column 11: read pair count reported by TopHat
# column 12: split read count supported by Blat (e.g., 1)
# column 13: split read count supported by Blat (considering split read supported by TopHat only, e.g., 1)
# column 14: minimal distance between read pair and breakpoint 1 after aligning by TopHat (e.g., 13)
# column 15: minimal distance between read pair and breakpoint 2 after aligning by TopHat (e.g., 22)
# column 16: sequence identity (e.g., 0.52)
# column 17: minimal distance between split read and breakpoint 1 after aligning by Blat (e.g., 0, use mean if both read 1 and read 2 are split read)
# column 18: minimal distance between split read and breakpoint 2 after aligning by Blat (e.g., 0, use mean if both read 1 and read 2 are split read)
# column 19: minimal distance between split read and breakpoint 1 after aligning by Blat (considering split read supported by TopHat only, e.g., 0, use mean if both read 1 and read 2 are split read)
# column 20: minimal distance between split read and breakpoint 2 after aligning by Blat (considering split read supported by TopHat only, e.g., 0, use mean if both read 1 and read 2 are split read)
# column 21: presence of canonical splice site
# column 22: cluster count around breakpoint 1 (e.g., 1)
# column 23: cluster count around breakpoint 2 (e.g., 1)
# column 24: cluster count around breakpoint 1 and breakpoint 2 (e.g., 1)
# column 25: percentage of split reads and read pairs supported fusion transcript’s cluster around breakpoint 1 and breakpoint 2 (e.g., 0.75)
# column 26: standard deviation for candidate fusion clusters around breakpoint 1 (e.g., NA, NA for only one cluster)
# column 27: standard deviation for candidate fusion clusters around breakpoint 2 (e.g., NA, NA for only one cluster)
# column 28: whether fusion transcript locates in the same gene (e.g., N)
# column 29: fusion type
# column 30: gene name of breakpoint 1
# column 31: gene type of breakpoint 1
# column 32: gene strand of breakpoint 1
# column 33: breakpoint location of breakpoint 1
# column 34: breakpoint region type of breakpoint 1
# column 35: Frame_extra_base of breakpoint 2
# column 36: gene name of breakpoint 2
# column 37: gene type of breakpoint 2
# column 38: gene strand of breakpoint 2
# column 39: breakpoint location of breakpoint 2
# column 40: breakpoint region type of breakpoint 2
# column 41: Frame_extra_base of breakpoint 2
# column 42: Fusion frame (in-frame, out-frame, or unknown)
# column 43: split reads reported by STAR (e.g., read_23)
# column 44: read pairs reported by STAR (e.g., read_38,read_62,read_70,read_8)
# column 45: split reads supported by TopHat (e.g., read_23)
# column 46: potential split reads not supported by TopHat (e.g., NA)
# column 47: read pairs supported by TopHat (e.g., read_38,read_70,read_8)
# column 48: split reads supported by Blat(blat tophat split read and tophat potential split read)
# column 49: distance between each read pair and breakpoint 1 after aligning by TopHat (e.g., 740,23,13)
# column 50: distance between each read pair and breakpoint 2 after aligning by TopHat (e.g., 292,2592,22)
# column 51: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 52: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 53: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 54: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 55: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 1 (e.g., --AGTGGGCCAGGTAG-GGCTGG)
# column 56: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 2 (e.g., CCACT--GCCAGG-AGAACCTCA)
# column 57: discordant read cluster IDs around breakpoint 1 (e.g., 1)
# column 58: discordant read cluster IDs around breakpoint 2 (e.g., 1)
# column 59: supporting reads count in each cluster around breakpoint 1 (e.g., 4)
# column 60: supporting reads count in each cluster around breakpoint 2 (e.g., 3)



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
    "Split_read_count_(star)", "Read_pair_count_(star)", "Split_read_count_(tophat)", "Potential_split_read_count_(tophat)", "Read_pair_count_(tophat)",
    "Split_read_count_(blat_tophat_split_and_tophat_potential_split_reads)", "Split_read_count_(blat_tophat_split_reads)",
    "Minimum_read_pair_distance_to_breakpoint_1", "Minimum_read_pair_distance_to_breakpoint_2",
    "Sequence_identity", "Minimum_blat_distance_to_breakpoint_1_(tophat_split_and_potential_split_reads)", "Minimum_blat_distance_to_breakpoint_2_(tophat_split_and_potential_split_reads)",
    "Minimum_blat_distance_to_breakpoint_1_(tophat_split_reads)", "Minimum_blat_distance_to_breakpoint_2_(tophat_split_reads)",
    "Have_canonical_motif",
    "Cluster_count_breakpoint_1", "Cluster_count_breakpoint_2", "Cluster_count_breakpoint_1_and_2",
    "(discordant_reads)%_support_fusion",
    "SD_(discordant_reads)%_in_breakpoint_1", "SD_(discordant_reads)%_in_breakpoint_2",
    "Whether_in_the_same_gene",
    "Fusion_type",
    "Gene_name_breakpoint_1", "Gene_type_breakpoint_1", "Gene_strand_breakpoint_1", "Breakpoint_location_breakpoint_1", "Breakpoint_region_type_breakpoint_1", "Frame_extra_base_breakpoint_1", 
				"Gene_name_breakpoint_2", "Gene_type_breakpoint_2", "Gene_strand_breakpoint_2", "Breakpoint_location_breakpoint_2", "Breakpoint_region_type_breakpoint_2", "Frame_extra_base_breakpoint_2", 
    "Fusion_frame",
    "Split_reads_(star)", "Read_pairs_(star)", "Split_reads_(tophat)", "Potential_split_reads_(tophat)", "Read_pairs_(tophat)",
    "Split_reads_(blat_tophat_split_and_tophat_potential_split_reads)",
    "Read_pair_distance_to_breakpoint_1", "Read_pair_distance_to_breakpoint_2",
    "Blat_distance_to_breakpoint_1_(tophat_split_and_tophat_potential_split_reads)", "Blat_distance_to_breakpoint_2_(tophat_split_and_tophat_potential_split_reads)",
    "Blat_distance_to_breakpoint_1_(tophat_split_reads)", "Blat_distance_to_breakpoint_2_(tophat_split_eads)",
    "Sequence_identity_alignment_breakpoint_1", "Sequence_identity_alignment_breakpoint_2",
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
        $split_read_count_star, $read_pair_count_star, $split_read_count_tophat, $potential_split_read_count_tophat, $read_pair_count_tophat, $split_read_count_blat, $split_read_count_blat_base_on_tophat_split_read,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2, $identity_output,
        $minimum_blat_distance_breakpoint_1, $minimum_blat_distance_breakpoint_2, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_1, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_2,
        $is_split_site_contain_canonical_motif,
        $total_clusters_breakpoint_1, $total_clusters_breakpoint_2, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
        $is_in_the_same_gene,
        $fusion_type,
        $gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1,
        $gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2,
        $fusion_frame,
        $split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat, $split_reads_blat,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2, 
        $blat_distance_breakpoint_1, $blat_distance_breakpoint_2, $blat_distance_base_on_tophat_split_read_breakpoint_1, $blat_distance_base_on_tophat_split_read_breakpoint_2,
        $identity_seq_breakpoint_1, $identity_seq_breakpoint_2,
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
        $split_read_count_star, $read_pair_count_star, $split_read_count_tophat, $potential_split_read_count_tophat, $read_pair_count_tophat, $split_read_count_blat, $split_read_count_blat_base_on_tophat_split_read,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2, $identity_output,
        $minimum_blat_distance_breakpoint_1, $minimum_blat_distance_breakpoint_2, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_1, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_2,
        $is_split_site_contain_canonical_motif,
        $total_clusters_breakpoint_1, $total_clusters_breakpoint_2, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
        $is_in_the_same_gene,
        $fusion_type,
        $gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1,
        $gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2,
        $fusion_frame,
        $split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat, $split_reads_blat,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2, 
        $blat_distance_breakpoint_1, $blat_distance_breakpoint_2, $blat_distance_base_on_tophat_split_read_breakpoint_1, $blat_distance_base_on_tophat_split_read_breakpoint_2,
        $identity_seq_breakpoint_1, $identity_seq_breakpoint_2,
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
