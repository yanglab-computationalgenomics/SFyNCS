#!/usr/bin/env perl

# 2022-03-17

# 1. Function
# Annotate fusions
# Filter fusion locate in the same gene

# 2. Input
# 2.1. Annotation file (must have header)
# column 1: chromosome
# column 2: start (1-base)
# column 3: end (1-base)
# column 4: strand
# column 5: symbol
# column 6: gene type

# 2.2. filtered fusions
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
# column 20: fusion annotations
# column 21: split reads (processed by tophat and blat)
# column 22: read pairs (processed by tophat)



use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;


GetOptions(
    'h|help'	=> sub{usage()}
)||usage();

system("awk 'BEGIN{FS=OFS=\"\t\"} NR>1{print \$1,\$2-1,\$3,\$4,\$5,\$6;}' $ARGV[0] | sort -k1,1 -k2,2n | uniq >temp_annotation.bed");
system("awk 'BEGIN{FS=OFS=\"\t\"} NR>1{print \$1,\$2-1,\$2; print \$4,\$5-1,\$5;}' $ARGV[1] | sort -k1,1 -k2,2n | uniq >temp_fusion_breakpoint.bed");
system("bedtools intersect -wa -wb -a temp_fusion_breakpoint.bed -b temp_annotation.bed | cut -f1,3- >temp_overlapped_gene.tsv");

open IN, "temp_overlapped_gene.tsv" or die "Can't open temp_overlapped_gene.tsv:$!";
my %breakpoint2gene;
while(<IN>){
    chomp;
    # symbol may contain "@", "." and "-"
    # gene type may contain "_" and "."
    my ($chr_fusion, $pos_fusion, $chr_gene, $start_gene, $end_gene, $strand_gene, $symbol, $gene_type)=split "\t", $_;
    my $temp_pos=join ";", ($chr_fusion, $pos_fusion);
    my $temp_gene_info=join ";", ($symbol, $gene_type, $strand_gene);
    $breakpoint2gene{$temp_pos}{$temp_gene_info}=0;
}
close(IN);

say join "\t", ("Chr_left", "Pos_left", "Strand_left", "Chr_right", "Pos_right", "Strand_right",
     "Split_read_count_(Tophat_and_Blat)", "Read_pair_count_(Tophat)",
     "Minimum_read_distance_to_left", "Minimum_read_distance_to_right",
     "Identity", "Minimum_blat_distance_to_left", "Minimum_blat_distance_to_right",
     "Total_clusters_(left)", "Total_clusters_(right)", "Total_clusters_(merge)",
     "(discordant_reads)%_support_fusion",
     "SD_(discordant_reads)%_in_clusters_(left)", "SD_(discordant_reads)%_in_clusters_(right)",
     "Fusion_annotations",
     "Split_reads_(tophat_and_blat)", "Read_pairs_(tophat)");
open IN, $ARGV[1] or die "Can't open $ARGV[1]:$!";
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
        $split_reads_blat, $read_pairs_tophat)=split "\t", $_;
    
    my $temp_key=join ";", ($chr_left, $pos_left);
    my @genes_left_info=sort keys %{$breakpoint2gene{$temp_key}};
    push @genes_left_info, "NA;NA;NA"; # add NA;NA;NA to have at least one for loop
    $temp_key=join ";", ($chr_right, $pos_right);
    my @genes_right_info=sort keys %{$breakpoint2gene{$temp_key}};
    push @genes_right_info, "NA;NA;NA"; # add NA;NA;NA to have at least one for loop
    
    
    my (@valid_fusions, @fusion_to_unannotated);
    my $is_in_the_same_gene=0; 
    foreach my $left_gene_info (@genes_left_info){
	$left_gene_info=~/(.*);(.*);(.*)/;
	my ($symbol_left, $gene_type_left, $strand_gene_left)=($1, $2, $3);
	foreach my $right_gene_info (@genes_right_info){
	    $right_gene_info=~/(.*);(.*);(.*)/;
	    my ($symbol_right, $gene_type_right, $strand_gene_right)=($1, $2, $3);
	    if($strand_left eq '+' && $strand_right eq '-' && $symbol_left eq $symbol_right && $symbol_left ne "NA"){
		$is_in_the_same_gene=1;
		goto OUTPUT;
	    }
	    
	    # 1. if $strand_left ne $strand_right, $strand_gene_left eq $strand_gene_right will be a valid fusion
	    # 2. if $strand_left eq $strand_right, $strand_gene_left ne $strand_gene_right will be a valid fusion
	    # 3. if $strand_left ne $strand_gene_left, switch left and right symbol
	    # symbol may contain "@", "." and "-"
	    # gene type may contain "_" and "."
	    next if $symbol_left eq "NA" && $symbol_right eq "NA";
	    if($symbol_left ne "NA" && $symbol_right ne "NA"){
		if($strand_left ne $strand_right){
		    if($strand_gene_left eq $strand_gene_right){ # valid fusion
			if($strand_left ne $strand_gene_left){ # swith left and right symbol
			    push @valid_fusions, $symbol_right."_"."(".$gene_type_right.")"."--".$symbol_left."_"."(".$gene_type_left.")";
			}else{
			    push @valid_fusions, $symbol_left."_"."(".$gene_type_left.")"."--".$symbol_right."_"."(".$gene_type_right.")";
			}
		    }else{
			push @fusion_to_unannotated, $symbol_left."_"."(".$gene_type_left.")"."--"."unannotated"."_"."("."unknown".")";
			push @fusion_to_unannotated, $symbol_right."_"."(".$gene_type_right.")"."--"."unannotated"."_"."("."unknown".")";
		    }
		}else{
		    if($strand_gene_left ne $strand_gene_right){ # valid fusion
			if($strand_left ne $strand_gene_left){ # swith left and right symbol
			    push @valid_fusions, $symbol_right."_"."(".$gene_type_right.")"."--".$symbol_left."_"."(".$gene_type_left.")";
			}else{
			    push @valid_fusions, $symbol_left."_"."(".$gene_type_left.")"."--".$symbol_right."_"."(".$gene_type_right.")";
			}
		    }else{
			push @fusion_to_unannotated, $symbol_left."_"."(".$gene_type_left.")"."--"."unannotated"."_"."("."unknown".")";
			push @fusion_to_unannotated, $symbol_right."_"."(".$gene_type_right.")"."--"."unannotated"."_"."("."unknown".")";
		    }
		}
	    }else{
		push @fusion_to_unannotated, $symbol_left."_"."(".$gene_type_left.")"."--"."unannotated"."_"."("."unknown".")" if $symbol_left ne "NA";
		push @fusion_to_unannotated, $symbol_right."_"."(".$gene_type_right.")"."--"."unannotated"."_"."("."unknown".")" if $symbol_right ne "NA";
	    }
	}
    }
    
    my $fusion_annotation="NA";
    my $validate_fusion_count=@valid_fusions;
    my $fusion_to_unannotated_count=@fusion_to_unannotated;
    if($validate_fusion_count==0 && $fusion_to_unannotated_count==0){
	$fusion_annotation="unannotated_(unknown)--unannotated_(unknown)";
    }else{
	if($validate_fusion_count>0){
	    foreach my $temp_value (@valid_fusions){
		$fusion_annotation.=";".$temp_value;
	    }
	}else{
	    foreach my $temp_value (@fusion_to_unannotated){
		$fusion_annotation.=";".$temp_value;
	    }
	}
	$fusion_annotation=~s/NA;//;
    }
    
    OUTPUT: say join "\t",($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right,
        $split_read_count_blat, $read_pair_count_tophat,
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $identity_output, $minimum_blat_distance_left, $minimum_blat_distance_right,
        $total_clusters_left, $total_clusters_right, $total_clusters,
        $percentage_support_fusion,
        $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right,
	$fusion_annotation,
        $split_reads_blat, $read_pairs_tophat) if $is_in_the_same_gene==0;
}

system('rm temp_annotation.bed temp_fusion_breakpoint.bed temp_overlapped_gene.tsv');



sub usage{
    my $scriptName = basename $0;
print <<HELP;
This script was used to annotate fusions
Usage: perl $scriptName annotation_file.tsv fusions.tsv >output
Options:

    -h --help     		print this help information
HELP
    exit(-1);
}


