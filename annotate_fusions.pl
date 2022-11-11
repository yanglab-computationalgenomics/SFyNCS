#!/usr/bin/env perl

# 2022-11-04

# 1. Function
# Annotate fusions
# select annotation: 1) both annotate > only one breakpoint annotate 2) splice_site > exon > intron 3) protein-coding > non-coding
# if more than 1 frame, in frame > out of frame > unknown frame
# Filter fusion locate in the same gene
# Swith breakpoint so fusion is in sense-sense orientation
# Genes were classify to two categories: protein_coding_gene and non_protein_coding_gene

# 2. Input
# 2.1. Annotation file (must have header, Gene Predictions Extended or gpe format)
# column 1: name of gene (usually transcript_id from GTF)
# column 2: chromosome
# column 3: strand
# column 4: transcription start position (0-base)
# column 5: transcription end position (1-base)
# column 6: coding region start (0-base)
# column 7: coding region end (1-base)
# column 8: number of exons
# column 9: exon start positions (0-base)
# column 10: exon end positions (1-base)
# column 11: score
# column 12: alternate name (e.g. gene_id from GTF)
# column 13: status of CDS start annotation (none, unknown, incomplete, or complete)
# column 14: status of CDS end annotation (none, unknown, incomplete, or complete)
# column 15: exon frame offsets {0,1,2}

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
# column 28: split reads reported by STAR (e.g., read_23)
# column 29: read pairs reported by STAR (e.g., read_38,read_62,read_70,read_8)
# column 30: split reads supported by TopHat (e.g., read_23)
# column 31: potential split reads not supported by TopHat (e.g., NA)
# column 32: read pairs supported by TopHat (e.g., read_38,read_70,read_8)
# column 33: split reads supported by Blat(blat tophat split read and tophat potential split read)
# column 34: distance between each read pair and breakpoint 1 after aligning by TopHat (e.g., 740,23,13)
# column 35: distance between each read pair and breakpoint 2 after aligning by TopHat (e.g., 292,2592,22)
# column 36: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 37: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 38: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 39: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 40: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 1 (e.g., --AGTGGGCCAGGTAG-GGCTGG)
# column 41: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 2 (e.g., CCACT--GCCAGG-AGAACCTCA)
# column 42: discordant read cluster IDs around breakpoint 1 (e.g., 1)
# column 43: discordant read cluster IDs around breakpoint 2 (e.g., 1)
# column 44: supporting reads count in each cluster around breakpoint 1 (e.g., 4)
# column 45: supporting reads count in each cluster around breakpoint 2 (e.g., 3)

# 3. Output
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
# column 35: Frame_extra_base of breakpoint 1
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

my $fasta;
my $filter_in_the_same_gene="Y";

GetOptions(
    'f|filter_in_the_same_gene=s'   =>  \$filter_in_the_same_gene,
    'a|fasta=s'                     => 	\$fasta,
    'h|help'	                       =>  sub{usage()}
)||usage();

# delete chr*_ by $2!~/_/
system("awk 'BEGIN{FS=OFS=\"\t\"} NR>1 && \$2!~/_/{print \$2,\$4,\$5,\$1,\$3;}' $ARGV[0] | sort -k1,1 -k2,2n | uniq >temp_annotation.bed");
system("awk 'BEGIN{FS=OFS=\"\t\"} NR>1{print \$1,\$2-1,\$2; print \$4,\$5-1,\$5;}' $ARGV[1] | sort -k1,1 -k2,2n | uniq >temp_fusion_breakpoint.bed");
system("bedtools intersect -wa -wb -a temp_fusion_breakpoint.bed -b temp_annotation.bed | cut -f1,3- >temp_overlapped_gene.tsv");

open IN, "temp_overlapped_gene.tsv" or die "Can't open temp_overlapped_gene.tsv:$!";
my %breakpoint2gene;
my %select_transcripts;
my %hash_gene_info;
while(<IN>){
    chomp;
    # symbol may contain "@", "." and "-"
    # gene type may contain "_" and "."
    my ($chr_fusion, $pos_fusion, $chr_gene, $start_gene, $end_gene, $transcript_id, $strand_gene)=split "\t", $_;
    my $temp_pos=join ";", ($chr_fusion, $pos_fusion);
    # one transcript_id may have two locus (NM_000071 in chr21:43053190-3075835,- and chr21:6444868-6467513,-)
    my $temp_gene_info=join ";", ($chr_gene, $start_gene, $end_gene, $transcript_id, $strand_gene);
    $breakpoint2gene{$temp_pos}{$temp_gene_info}=0;
    $select_transcripts{$temp_gene_info}=0;
}
close(IN);

open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($transcript_id, $chr, $strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames)=split "\t", $_;
    my $gene_type=($cds_start==$cds_end) ? "non_coding_gene" : "protein_coding_gene";
				if((! exists $hash_gene_info{"gene_type"}{$name_2}) || ($hash_gene_info{"gene_type"}{$name_2} ne "protein_coding_gene")){
								$hash_gene_info{"gene_type"}{$name_2}=$gene_type;
								$hash_gene_info{"strand"}{$name_2}=$strand;
				}
    
    my $temp_gene_info=join ";", ($chr, $transcript_start, $transcript_end, $transcript_id, $strand);
    next if ! exists $select_transcripts{$temp_gene_info};
    $select_transcripts{$temp_gene_info}=$_;
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
open IN, $ARGV[1] or die "Can't open $ARGV[1]:$!";
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
        $split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat, $split_reads_blat,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2, 
        $blat_distance_breakpoint_1, $blat_distance_breakpoint_2, $blat_distance_base_on_tophat_split_read_breakpoint_1, $blat_distance_base_on_tophat_split_read_breakpoint_2,
        $identity_seq_breakpoint_1, $identity_seq_breakpoint_2,
        $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
        $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2)=split "\t", $_;
    
    # is contain canonical split motif
    my $flanking_search_sequence_breakpoint_1=&getFlankingSequence($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, 5, 5, "upstream");
				#$flanking_search_sequence_breakpoint_1=~s/-//g;
    #$flanking_search_sequence_breakpoint_1=substr($flanking_search_sequence_breakpoint_1, 5, 10);
    my $is_breakpoint_1_contain_donor_motif=($flanking_search_sequence_breakpoint_1=~/GT/) ? "Y" : "N";
    $flanking_search_sequence_breakpoint_1=~tr/[ACTG]/[TGAC]/;
    $flanking_search_sequence_breakpoint_1=reverse($flanking_search_sequence_breakpoint_1);
    my $is_breakpoint_1_contain_acceptor_motif=($flanking_search_sequence_breakpoint_1=~/[CTA]AG/) ? "Y" : "N";
    
    my $flanking_search_sequence_breakpoint_2=&getFlankingSequence($input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2, 5, 5, "downstream");
				#$flanking_search_sequence_breakpoint_2=~s/-//g;
    #$flanking_search_sequence_breakpoint_2=substr($flanking_search_sequence_breakpoint_2, 5, 10);
    my $is_breakpoint_2_contain_acceptor_motif=($flanking_search_sequence_breakpoint_2=~/[CTA]AG/) ? "Y" : "N";
    $flanking_search_sequence_breakpoint_2=~tr/[ACTG]/[TGAC]/;
    $flanking_search_sequence_breakpoint_2=reverse($flanking_search_sequence_breakpoint_2);
    my $is_breakpoint_2_contain_donor_motif=($flanking_search_sequence_breakpoint_2=~/GT/) ? "Y" : "N";
    
    
    my $fusion_type_must_be_non_coding_as_unable_to_select=0;
				#my ($is_breakpoint_1_contain_donor_motif, $is_breakpoint_2_contain_acceptor_motif, $is_breakpoint_2_contain_donor_motif, $is_breakpoint_1_contain_acceptor_motif)=("Y") x 4;
				my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
				my ($gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1,
        $gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2);
				if(($is_breakpoint_1_contain_donor_motif eq "Y" && $is_breakpoint_2_contain_acceptor_motif eq "Y") &&
								($is_breakpoint_2_contain_donor_motif eq "Y" && $is_breakpoint_1_contain_acceptor_motif eq "Y")){
								my ($gene_name_breakpoint_1_1, $gene_type_breakpoint_1_1, $gene_strand_breakpoint_1_1, $fusion_pos_location_breakpoint_1_1, $fusion_pos_region_type_breakpoint_1_1, $fusion_extra_base_breakpoint_1_1)=&get_gene_name_and_type($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, "donor");
								my ($gene_name_breakpoint_2_1, $gene_type_breakpoint_2_1, $gene_strand_breakpoint_2_1, $fusion_pos_location_breakpoint_2_1, $fusion_pos_region_type_breakpoint_2_1, $fusion_extra_base_breakpoint_2_1)=&get_gene_name_and_type($input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2, "acceptor");
								my ($gene_name_breakpoint_1_2, $gene_type_breakpoint_1_2, $gene_strand_breakpoint_1_2, $fusion_pos_location_breakpoint_1_2, $fusion_pos_region_type_breakpoint_1_2, $fusion_extra_base_breakpoint_1_2)=&get_gene_name_and_type($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, "acceptor");
								my ($gene_name_breakpoint_2_2, $gene_type_breakpoint_2_2, $gene_strand_breakpoint_2_2, $fusion_pos_location_breakpoint_2_2, $fusion_pos_region_type_breakpoint_2_2, $fusion_extra_base_breakpoint_2_2)=&get_gene_name_and_type($input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2, "donor");
								my ($select_index, $able_to_select);
								
								if($gene_name_breakpoint_1_1 eq "Unannotated" && $gene_name_breakpoint_2_1 eq "Unannotated" && $gene_name_breakpoint_1_2 eq "Unannotated" && $gene_name_breakpoint_2_2 eq "Unannotated"){ # all unannotated
												$select_index=1;
												$able_to_select="N";
								}elsif($gene_name_breakpoint_1_1 ne "Unannotated" && $gene_name_breakpoint_2_1 ne "Unannotated" && $gene_name_breakpoint_1_2 ne "Unannotated" && $gene_name_breakpoint_2_2 ne "Unannotated"){ # all annotated
												# select splice site
												my ($have_splice_site_in_1, $have_splice_site_in_2)=(0) x 2;
												$have_splice_site_in_1=1 if $fusion_pos_location_breakpoint_1_1=~/splice_site/ || $fusion_pos_location_breakpoint_2_1=~/splice_site/;
												$have_splice_site_in_2=1 if $fusion_pos_location_breakpoint_1_2=~/splice_site/ || $fusion_pos_location_breakpoint_2_2=~/splice_site/;
												if($have_splice_site_in_1==1 && $have_splice_site_in_2==1){
																my ($both_have_splice_site_in_1, $both_have_splice_site_in_2)=(0) x 2;
																$both_have_splice_site_in_1=1 if $fusion_pos_location_breakpoint_1_1=~/splice_site/ && $fusion_pos_location_breakpoint_2_1=~/splice_site/;
																$both_have_splice_site_in_2=1 if $fusion_pos_location_breakpoint_1_2=~/splice_site/ && $fusion_pos_location_breakpoint_2_2=~/splice_site/;
																if($both_have_splice_site_in_1==1 && $both_have_splice_site_in_2==1){
																				my ($have_protein_in_1, $have_protein_in_2)=(0) x 2;
																				$have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ || $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																				$have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ || $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																				if($have_protein_in_1==1 && $have_protein_in_2==1){
																								my ($both_have_protein_in_1, $both_have_protein_in_2)=(0) x 2;
																								$both_have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ && $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																								$both_have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ && $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																								if($both_have_protein_in_1==1 && $both_have_protein_in_2==1){
																												$select_index=1;
																												$able_to_select="N";
																								}elsif($both_have_protein_in_1==1){
																												$select_index=1;
																												$able_to_select="Y";
																								}elsif($both_have_protein_in_2==1){
																												$select_index=2;
																												$able_to_select="Y";
																								}else{
																												$select_index=1;
																												$able_to_select="N";
																								}
																				}elsif($have_protein_in_1==1){
																								$select_index=1;
																								$able_to_select="Y";
																				}elsif($have_protein_in_2==1){
																								$select_index=2;
																								$able_to_select="Y";
																				}else{
																								$select_index=1;
																								$able_to_select="N";
																				}
																}elsif($both_have_splice_site_in_1==1){
																				$select_index=1;
																				$able_to_select="Y";
																}elsif($both_have_splice_site_in_2==1){
																				$select_index=2;
																				$able_to_select="Y";
																}else{
																				my ($have_protein_in_1, $have_protein_in_2)=(0) x 2;
																				$have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ || $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																				$have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ || $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																				if($have_protein_in_1==1 && $have_protein_in_2==1){
																								my ($both_have_protein_in_1, $both_have_protein_in_2)=(0) x 2;
																								$both_have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ && $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																								$both_have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ && $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																								if($both_have_protein_in_1==1 && $both_have_protein_in_2==1){
																												$select_index=1;
																												$able_to_select="N";
																								}elsif($both_have_protein_in_1==1){
																												$select_index=1;
																												$able_to_select="Y";
																								}elsif($both_have_protein_in_2==1){
																												$select_index=2;
																												$able_to_select="Y";
																								}else{
																												$select_index=1;
																												$able_to_select="N";
																								}
																				}elsif($have_protein_in_1==1){
																								$select_index=1;
																								$able_to_select="Y";
																				}elsif($have_protein_in_2==1){
																								$select_index=2;
																								$able_to_select="Y";
																				}else{
																								$select_index=1;
																								$able_to_select="N";
																				}
																}
												}elsif($have_splice_site_in_1==1){
																$select_index=1;
																$able_to_select="Y";
												}elsif($have_splice_site_in_2==1){
																$select_index=2;
																$able_to_select="Y";
												}else{
																# select base on exon
																my ($have_exon_in_1, $have_exon_in_2)=(0) x 2;
																$have_exon_in_1=1 if $fusion_pos_location_breakpoint_1_1=~/exon/ || $fusion_pos_location_breakpoint_2_1=~/exon/;
																$have_exon_in_2=1 if $fusion_pos_location_breakpoint_1_2=~/exon/ || $fusion_pos_location_breakpoint_2_2=~/exon/;
																if($have_exon_in_1==1 && $have_exon_in_2==1){
																				my ($both_have_exon_in_1, $both_have_exon_in_2)=(0) x 2;
																				$both_have_exon_in_1=1 if $fusion_pos_location_breakpoint_1_1=~/exon/ && $fusion_pos_location_breakpoint_2_1=~/exon/;
																				$both_have_exon_in_2=1 if $fusion_pos_location_breakpoint_1_2=~/exon/ && $fusion_pos_location_breakpoint_2_2=~/exon/;
																				if($both_have_exon_in_1==1 && $both_have_exon_in_2==1){
																								my ($have_protein_in_1, $have_protein_in_2)=(0) x 2;
																								$have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ || $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																								$have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ || $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																								if($have_protein_in_1==1 && $have_protein_in_2==1){
																												my ($both_have_protein_in_1, $both_have_protein_in_2)=(0) x 2;
																												$both_have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ && $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																												$both_have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ && $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																												if($both_have_protein_in_1==1 && $both_have_protein_in_2==1){
																																$select_index=1;
																																$able_to_select="N";
																												}elsif($both_have_protein_in_1==1){
																																$select_index=1;
																																$able_to_select="Y";
																												}elsif($both_have_protein_in_2==1){
																																$select_index=2;
																																$able_to_select="Y";
																												}else{
																																$select_index=1;
																																$able_to_select="N";
																												}
																								}elsif($have_protein_in_1==1){
																												$select_index=1;
																												$able_to_select="Y";
																								}elsif($have_protein_in_2==1){
																												$select_index=2;
																												$able_to_select="Y";
																								}else{
																												$select_index=1;
																												$able_to_select="N";
																								}
																				}elsif($both_have_exon_in_1==1){
																								$select_index=1;
																								$able_to_select="Y";
																				}elsif($both_have_exon_in_2==1){
																								$select_index=2;
																								$able_to_select="Y";
																				}else{
																								$select_index=1;
																								$able_to_select="N";
																				}
																}elsif($have_exon_in_1==1){
																				$select_index=1;
																				$able_to_select="Y";
																}elsif($have_exon_in_2==1){
																				$select_index=2;
																				$able_to_select="Y";
																}else{
																				$select_index=1;
																				$able_to_select="N";
																}
												}
								}else{
												my ($both_annotated_in_1, $both_annotated_in_2)=(0) x 2;
												$both_annotated_in_1=1 if $gene_name_breakpoint_1_1 ne "Unannotated" && $gene_name_breakpoint_2_1 ne "Unannotated";
												$both_annotated_in_2=1 if $gene_name_breakpoint_1_2 ne "Unannotated" && $gene_name_breakpoint_2_2 ne "Unannotated";
												if($both_annotated_in_1==1){ # it can't be $both_annotated_in_1==1 and $both_annotated_in_2==1
																$select_index=1;
																$able_to_select="Y";
												}elsif($both_annotated_in_2==1){
																$select_index=2;
																$able_to_select="Y";
												}else{ # at least one annotated
																my ($both_unannotated_in_1, $both_unannotated_in_2)=(0) x 2;
																$both_unannotated_in_1=1 if $gene_name_breakpoint_1_1 eq "Unannotated" && $gene_name_breakpoint_2_1 eq "Unannotated";
																$both_unannotated_in_2=1 if $gene_name_breakpoint_1_2 eq "Unannotated" && $gene_name_breakpoint_2_2 eq "Unannotated";
																if($both_unannotated_in_1==1){
																				$select_index=2;
																				$able_to_select="Y";
																}elsif($both_unannotated_in_2==1){
																				$select_index=1;
																				$able_to_select="Y";
																}else{ # both with one breakpoint annotated and another not annotated
																				# select splice site
																				my ($have_splice_site_in_1, $have_splice_site_in_2)=(0) x 2;
																				$have_splice_site_in_1=1 if $fusion_pos_location_breakpoint_1_1=~/splice_site/ || $fusion_pos_location_breakpoint_2_1=~/splice_site/;
																				$have_splice_site_in_2=1 if $fusion_pos_location_breakpoint_1_2=~/splice_site/ || $fusion_pos_location_breakpoint_2_2=~/splice_site/;
																				if($have_splice_site_in_1==1 && $have_splice_site_in_2==1){
																								my ($have_protein_in_1, $have_protein_in_2)=(0) x 2;
																								$have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ || $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																								$have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ || $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																								if($have_protein_in_1==1 && $have_protein_in_2==1){
																												$select_index=1;
																												$able_to_select="N";
																								}elsif($have_protein_in_1==1){
																												$select_index=1;
																												$able_to_select="Y";
																								}elsif($have_protein_in_2==1){
																												$select_index=2;
																												$able_to_select="Y";
																								}else{
																												$select_index=1;
																												$able_to_select="N";
																								}
																				}elsif($have_splice_site_in_1==1){
																								$select_index=1;
																								$able_to_select="Y";
																				}elsif($have_splice_site_in_2==1){
																								$select_index=2;
																								$able_to_select="Y";
																				}else{
																								# select base on exon
																								my ($have_exon_in_1, $have_exon_in_2)=(0) x 2;
																								$have_exon_in_1=1 if $fusion_pos_location_breakpoint_1_1=~/exon/ || $fusion_pos_location_breakpoint_2_1=~/exon/;
																								$have_exon_in_2=1 if $fusion_pos_location_breakpoint_1_2=~/exon/ || $fusion_pos_location_breakpoint_2_2=~/exon/;
																								if($have_exon_in_1==1 && $have_exon_in_2==1){
																												my ($have_protein_in_1, $have_protein_in_2)=(0) x 2;
																												$have_protein_in_1=1 if $gene_type_breakpoint_1_1=~/protein_coding_gene/ || $gene_type_breakpoint_2_1=~/protein_coding_gene/;
																												$have_protein_in_2=1 if $gene_type_breakpoint_1_2=~/protein_coding_gene/ || $gene_type_breakpoint_2_2=~/protein_coding_gene/;
																												if($have_protein_in_1==1 && $have_protein_in_2==1){
																																$select_index=1;
																																$able_to_select="N";
																												}elsif($have_protein_in_1==1){
																																$select_index=1;
																																$able_to_select="Y";
																												}elsif($have_protein_in_2==1){
																																$select_index=2;
																																$able_to_select="Y";
																												}else{
																																$select_index=1;
																																$able_to_select="N";
																												}
																								}elsif($have_exon_in_1==1){
																												$select_index=1;
																												$able_to_select="Y";
																								}elsif($have_exon_in_2==1){
																												$select_index=2;
																												$able_to_select="Y";
																								}else{
																												$select_index=1;
																												$able_to_select="N";
																								}
																				}
																}
												}
								}
								
								if($select_index==1){
												($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2)=($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, $input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2);
												($gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1)=($gene_name_breakpoint_1_1, $gene_type_breakpoint_1_1, $gene_strand_breakpoint_1_1, $fusion_pos_location_breakpoint_1_1, $fusion_pos_region_type_breakpoint_1_1, $fusion_extra_base_breakpoint_1_1);
												($gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2)=($gene_name_breakpoint_2_1, $gene_type_breakpoint_2_1, $gene_strand_breakpoint_2_1, $fusion_pos_location_breakpoint_2_1, $fusion_pos_region_type_breakpoint_2_1, $fusion_extra_base_breakpoint_2_1);
												
												if($able_to_select ne "Y"){
																if($gene_name_breakpoint_1_2 ne "Unannotated"){
																				if($gene_name_breakpoint_1 ne "Unannotated"){
																								$gene_name_breakpoint_1.="/".$gene_name_breakpoint_1_2;
																								$gene_type_breakpoint_1.="/".$gene_type_breakpoint_1_2;
																								$gene_strand_breakpoint_1.="/".$gene_strand_breakpoint_1_2;
																								$fusion_pos_location_breakpoint_1.="/".$fusion_pos_location_breakpoint_1_2;
                        $fusion_pos_region_type_breakpoint_1.="/".$fusion_pos_region_type_breakpoint_1_2;
																				}else{
																								$gene_name_breakpoint_1=$gene_name_breakpoint_1_2;
																								$gene_type_breakpoint_1=$gene_type_breakpoint_1_2;
																								$gene_strand_breakpoint_1=$gene_strand_breakpoint_1_2;
																								$fusion_pos_location_breakpoint_1=$fusion_pos_location_breakpoint_1_2;
                        $fusion_pos_region_type_breakpoint_1=$fusion_pos_region_type_breakpoint_1_2;
																				}
																				
																}
																if($gene_name_breakpoint_2_2 ne "Unannotated"){
																				if($gene_name_breakpoint_2 ne "Unannotated"){
																								$gene_name_breakpoint_2.="/".$gene_name_breakpoint_2_2;
																								$gene_type_breakpoint_2.="/".$gene_type_breakpoint_2_2;
																								$gene_strand_breakpoint_2.="/".$gene_strand_breakpoint_2_2;
																								$fusion_pos_location_breakpoint_2.="/".$fusion_pos_location_breakpoint_2_2;
                        $fusion_pos_region_type_breakpoint_2.="/".$fusion_pos_region_type_breakpoint_2_2;
																				}else{
																								$gene_name_breakpoint_2=$gene_name_breakpoint_2_2;
																								$gene_type_breakpoint_2=$gene_type_breakpoint_2_2;
																								$gene_strand_breakpoint_2=$gene_strand_breakpoint_2_2;
																								$fusion_pos_location_breakpoint_2=$fusion_pos_location_breakpoint_2_2;
                        $fusion_pos_region_type_breakpoint_2=$fusion_pos_region_type_breakpoint_2_2;
																				}
																}
                $fusion_extra_base_breakpoint_1="N";
                $fusion_extra_base_breakpoint_2="N";
												}
								}else{
												($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2)=($input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2, $input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1);
            ($min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
                $minimum_blat_distance_breakpoint_1, $minimum_blat_distance_breakpoint_2,
                $minimum_blat_distance_base_on_tophat_split_read_breakpoint_1, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_2,
                $total_clusters_breakpoint_1, $total_clusters_breakpoint_2,
                $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
                $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2,
                $blat_distance_breakpoint_1, $blat_distance_breakpoint_2,
                $blat_distance_base_on_tophat_split_read_breakpoint_1, $blat_distance_base_on_tophat_split_read_breakpoint_2,
                $identity_seq_breakpoint_1, $identity_seq_breakpoint_2,
                $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
                $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2)=($min_read_pair_distance_breakpoint_2, $min_read_pair_distance_breakpoint_1,
                    $minimum_blat_distance_breakpoint_2, $minimum_blat_distance_breakpoint_1, 
                    $minimum_blat_distance_base_on_tophat_split_read_breakpoint_2, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_1, 
                    $total_clusters_breakpoint_2, $total_clusters_breakpoint_1, 
                    $sd_percentage_discorant_read_in_each_cluster_breakpoint_2, $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, 
                    $read_pair_distance_to_breakpoint_2, $read_pair_distance_to_breakpoint_1, 
                    $blat_distance_breakpoint_2, $blat_distance_breakpoint_1, 
                    $blat_distance_base_on_tophat_split_read_breakpoint_2, $blat_distance_base_on_tophat_split_read_breakpoint_1, 
                    $identity_seq_breakpoint_2, $identity_seq_breakpoint_1, 
                    $cluster_ids_breakpoint_2, $cluster_ids_breakpoint_1, 
                    $total_reads_in_each_cluster_breakpoint_2, $total_reads_in_each_cluster_breakpoint_1);
            ($gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1)=($gene_name_breakpoint_2_2, $gene_type_breakpoint_2_2, $gene_strand_breakpoint_2_2, $fusion_pos_location_breakpoint_2_2, $fusion_pos_region_type_breakpoint_2_2, $fusion_extra_base_breakpoint_2_2);
												($gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2)=($gene_name_breakpoint_1_2, $gene_type_breakpoint_1_2, $gene_strand_breakpoint_1_2, $fusion_pos_location_breakpoint_1_2, $fusion_pos_region_type_breakpoint_1_2, $fusion_extra_base_breakpoint_1_2);
								}
        # output to STDERR to debug
								#if($able_to_select eq "Y"){
								#				say STDERR "able_to_select:";
								#				say STDERR $_;
								#}else{
        #    $fusion_type_must_be_non_coding_as_unable_to_select=1;
								#				say STDERR "unable_to_select:";
								#				say STDERR $_;
								#}
				}elsif($is_breakpoint_1_contain_donor_motif eq "Y" && $is_breakpoint_2_contain_acceptor_motif eq "Y"){
								($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2)=($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, $input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2);
								($gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1)=&get_gene_name_and_type($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, "donor");
								($gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2)=&get_gene_name_and_type($chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, "acceptor");
				}elsif($is_breakpoint_2_contain_donor_motif eq "Y" && $is_breakpoint_1_contain_acceptor_motif eq "Y"){
								($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2)=($input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2, $input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1);
								($min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
                $minimum_blat_distance_breakpoint_1, $minimum_blat_distance_breakpoint_2,
                $minimum_blat_distance_base_on_tophat_split_read_breakpoint_1, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_2,
                $total_clusters_breakpoint_1, $total_clusters_breakpoint_2,
                $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
                $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2,
                $blat_distance_breakpoint_1, $blat_distance_breakpoint_2,
                $blat_distance_base_on_tophat_split_read_breakpoint_1, $blat_distance_base_on_tophat_split_read_breakpoint_2,
                $identity_seq_breakpoint_1, $identity_seq_breakpoint_2,
                $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
                $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2)=($min_read_pair_distance_breakpoint_2, $min_read_pair_distance_breakpoint_1,
                    $minimum_blat_distance_breakpoint_2, $minimum_blat_distance_breakpoint_1, 
                    $minimum_blat_distance_base_on_tophat_split_read_breakpoint_2, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_1, 
                    $total_clusters_breakpoint_2, $total_clusters_breakpoint_1, 
                    $sd_percentage_discorant_read_in_each_cluster_breakpoint_2, $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, 
                    $read_pair_distance_to_breakpoint_2, $read_pair_distance_to_breakpoint_1, 
                    $blat_distance_breakpoint_2, $blat_distance_breakpoint_1, 
                    $blat_distance_base_on_tophat_split_read_breakpoint_2, $blat_distance_base_on_tophat_split_read_breakpoint_1, 
                    $identity_seq_breakpoint_2, $identity_seq_breakpoint_1, 
                    $cluster_ids_breakpoint_2, $cluster_ids_breakpoint_1, 
                    $total_reads_in_each_cluster_breakpoint_2, $total_reads_in_each_cluster_breakpoint_1);
        ($gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1)=&get_gene_name_and_type($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, "donor");
								($gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2)=&get_gene_name_and_type($chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, "acceptor");
				}else{
								($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2)=($input_chr_breakpoint_1, $input_pos_breakpoint_1, $input_strand_breakpoint_1, $input_chr_breakpoint_2, $input_pos_breakpoint_2, $input_strand_breakpoint_2);
								($gene_name_breakpoint_1, $gene_type_breakpoint_1, $gene_strand_breakpoint_1, $fusion_pos_location_breakpoint_1, $fusion_pos_region_type_breakpoint_1, $fusion_extra_base_breakpoint_1)=&get_gene_name_and_type($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, "donor");
								($gene_name_breakpoint_2, $gene_type_breakpoint_2, $gene_strand_breakpoint_2, $fusion_pos_location_breakpoint_2, $fusion_pos_region_type_breakpoint_2, $fusion_extra_base_breakpoint_2)=&get_gene_name_and_type($chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, "acceptor");
				}
				
    # whether in the same gene
    my $is_in_the_same_gene="N";
    if($input_chr_breakpoint_1 eq $input_chr_breakpoint_2){
        my @Temp_gene_name_breakpoint_1=();
        my @Temp_gene_name_breakpoint_2=();
        my $temp_key=join ";", ($chr_breakpoint_1, $pos_breakpoint_1);
        my @genes_info=sort keys %{$breakpoint2gene{$temp_key}};
        my %temp_hash;
        foreach my $gene_info (@genes_info){
            my ($transcript_id, $transcript_chr, $transcript_strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames)=split "\t", $select_transcripts{$gene_info};
            next if exists $temp_hash{$name_2};
            push @Temp_gene_name_breakpoint_1, $name_2;
            $temp_hash{$name_2}=0;   
        }
        $temp_key=join ";", ($chr_breakpoint_2, $pos_breakpoint_2);
        @genes_info=sort keys %{$breakpoint2gene{$temp_key}};
        %temp_hash=();
        foreach my $gene_info (@genes_info){
            my ($transcript_id, $transcript_chr, $transcript_strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames)=split "\t", $select_transcripts{$gene_info};
            next if exists $temp_hash{$name_2};
            push @Temp_gene_name_breakpoint_2, $name_2;
            $temp_hash{$name_2}=0;   
        }
        
        foreach my $gene_1 (@Temp_gene_name_breakpoint_1){
            last if $is_in_the_same_gene eq "Y";
            foreach my $gene_2 (@Temp_gene_name_breakpoint_2){
                if($gene_1 eq $gene_2){
                    $is_in_the_same_gene="Y";
                }
                last if $is_in_the_same_gene eq "Y";
            }
        }
    }
    
    next if $filter_in_the_same_gene eq "Y" && $is_in_the_same_gene eq "Y";
    
    # fusion frame
				my $fusion_type="non-coding";
				$fusion_type="protein-coding" if $gene_type_breakpoint_1=~/protein_coding_gene/ && $gene_type_breakpoint_2=~/protein_coding_gene/ ;
				$fusion_type="non-coding" if $fusion_type_must_be_non_coding_as_unable_to_select==1;
    my $fusion_frame="unknown_frame";
    if($fusion_type_must_be_non_coding_as_unable_to_select!=1){
        my @Fusion_extra_base_breakpoint_1=split "/|;", $fusion_extra_base_breakpoint_1;
        my @Fusion_extra_base_breakpoint_2=split "/|;", $fusion_extra_base_breakpoint_2;
        foreach my $temp_extra_base_breakpoint_1 (@Fusion_extra_base_breakpoint_1){
            next if $temp_extra_base_breakpoint_1 eq "N";
            foreach my $temp_extra_base_breakpoint_2 (@Fusion_extra_base_breakpoint_2){
                next if $temp_extra_base_breakpoint_2 eq "N";
                if($temp_extra_base_breakpoint_1 ne $temp_extra_base_breakpoint_2){
                    $fusion_frame="out_of_frame";
                }else{
                    $fusion_frame="in_frame";
                    last;
                }
            }
            last if $fusion_frame eq "in_frame";
        }
    }
    
    # change "protein_coding" to "protein-coding" "non_coding" to "non-coding"
				$gene_type_breakpoint_1=~s/protein_coding/protein-coding/g;
				$gene_type_breakpoint_1=~s/non_coding/non-coding/g;
    $fusion_pos_location_breakpoint_1=~s/protein_coding/protein-coding/g;
    $fusion_pos_location_breakpoint_1=~s/non_coding/non-coding/g;
    $fusion_pos_region_type_breakpoint_1=~s/protein_coding/protein-coding/g;
    $fusion_pos_region_type_breakpoint_1=~s/non_coding/non-coding/g;
				$gene_type_breakpoint_2=~s/protein_coding/protein-coding/g;
				$gene_type_breakpoint_2=~s/non_coding/non-coding/g;
    $fusion_pos_location_breakpoint_2=~s/protein_coding/protein-coding/g;
    $fusion_pos_location_breakpoint_2=~s/non_coding/non-coding/g;
    $fusion_pos_region_type_breakpoint_2=~s/protein_coding/protein-coding/g;
    $fusion_pos_region_type_breakpoint_2=~s/non_coding/non-coding/g;
    
    say join "\t",($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2,
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

system('rm -f temp_annotation.bed temp_fusion_breakpoint.bed temp_overlapped_gene.tsv');

sub get_gene_name_and_type{
				my ($chr_breakpoint, $pos_breakpoint, $strand_breakpoint, $donor_or_acceptor)=@_;
				my $temp_key=join ";", ($chr_breakpoint, $pos_breakpoint);
				my %temp_hash;
				if(exists $breakpoint2gene{$temp_key}){
								my @genes_info=sort keys %{$breakpoint2gene{$temp_key}};
								foreach my $gene_info (@genes_info){
												my ($transcript_id, $transcript_chr, $transcript_strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames)=split "\t", $select_transcripts{$gene_info};
												
												my $is_sense_orientation=0;
												if($donor_or_acceptor eq "donor"){
																$is_sense_orientation=1 if $strand_breakpoint eq $transcript_strand;
												}else{ # acceptor
																$is_sense_orientation=1 if $strand_breakpoint ne $transcript_strand;
												}
												next if $is_sense_orientation==0;
												
												my @Exon_starts=split ",", $exon_starts;
												my @Exon_ends=split ",", $exon_ends;
												my @Exon_frames=split ",", $exon_frames;
												my $fusion_position_location;
												my $fusion_location_index;
												for(my $i=0; $i<$exon_count; $i++){
																last if $Exon_starts[$i]+1>$pos_breakpoint;
																if($i<$exon_count-1){
																				if($pos_breakpoint>$Exon_ends[$i] && $pos_breakpoint<$Exon_starts[$i+1]+1){ # in intron
																								$fusion_position_location="intron";
																								$fusion_location_index=$i;
																								last;
																				}
																}
																next if $pos_breakpoint>$Exon_ends[$i];
																$fusion_location_index=$i;
																if($pos_breakpoint==$Exon_starts[$i]+1){
																				$fusion_position_location="splice_site";
																}elsif($pos_breakpoint==$Exon_ends[$i]){
																				$fusion_position_location="splice_site";
																}else{
																				$fusion_position_location="exon";
																}
												}
												
												my $gene_type=($cds_start==$cds_end) ? "non_coding_gene" : "protein_coding_gene";
												my $fusion_region_type="non_coding";
												if($gene_type eq "protein_coding_gene" && $fusion_position_location ne "intron"){
																if($pos_breakpoint>=$cds_start+1 && $pos_breakpoint<=$cds_end){
																				$fusion_region_type="CDS";
																}else{
																				if($pos_breakpoint<$cds_start+1){
																								$fusion_region_type=($transcript_strand eq '+') ? "5'UTR" : "3'UTR";
																				}else{
																								$fusion_region_type=($transcript_strand eq '+') ? "3'UTR" : "5'UTR";
																				}
																}
												}
												
												my $temp_extra_bases="N";
												if($gene_type eq "protein_coding_gene" && $fusion_position_location ne "intron"){
                if($fusion_region_type eq "CDS"){
                    if($donor_or_acceptor eq "donor"){
                        if($transcript_strand eq '+'){ # gene strand
                            if($Exon_starts[$fusion_location_index]<$cds_start){
                                $temp_extra_bases=$pos_breakpoint-$cds_start+$Exon_frames[$fusion_location_index];
                            }else{
                                $temp_extra_bases=$pos_breakpoint-$Exon_starts[$fusion_location_index]+$Exon_frames[$fusion_location_index];
                            }
                        }else{
                            if($Exon_ends[$fusion_location_index]>$cds_end){
                                $temp_extra_bases=$cds_end-$pos_breakpoint+1+$Exon_frames[$fusion_location_index];
                            }else{
                                $temp_extra_bases=$Exon_ends[$fusion_location_index]-$pos_breakpoint+1+$Exon_frames[$fusion_location_index];
                            }
                        }
                    }else{ # acceptor
                        if($transcript_strand eq '+'){ # gene strand
                            if($Exon_starts[$fusion_location_index]<$cds_start){
                                $temp_extra_bases=$pos_breakpoint-$cds_start-1+$Exon_frames[$fusion_location_index];
                            }else{
                                $temp_extra_bases=$pos_breakpoint-$Exon_starts[$fusion_location_index]-1+$Exon_frames[$fusion_location_index];
                            }
                        }else{
                            if($Exon_ends[$fusion_location_index]>$cds_end){
                                $temp_extra_bases=$cds_end-$pos_breakpoint+$Exon_frames[$fusion_location_index];
                            }else{
                                $temp_extra_bases=$Exon_ends[$fusion_location_index]-$pos_breakpoint+$Exon_frames[$fusion_location_index];
                            } 
                        }
                    }
                    $temp_extra_bases=$temp_extra_bases%3;
                }
												}
												
												# fusion_position_location=("intron", "exon", "splice_site") extra_base=("N", "0", "1", "2")
												# fusion_region_type=("5'UTR", "3'UTR", "CDS", "non_coding")
												$temp_hash{"position"}{$fusion_position_location}{$name_2}{$fusion_region_type}=0;
												$temp_hash{"gene"}{$name_2}{$gene_type}=0;
            $temp_hash{"extra_base"}{$name_2}{$temp_extra_bases}=0;
								}
				}
				
				my ($final_gene_name, $final_gene_type, $final_gene_strand, $final_fusion_position_location, $final_fusion_region_type, $final_extra_base)=("NA") x 6;
				# annotate order: skip intron if there is splice_site or exon	
				if(exists $temp_hash{"position"}{"splice_site"}){
								foreach my $temp_gene_name (sort keys %{$temp_hash{"position"}{"splice_site"}}){
												$final_gene_name.="/".$temp_gene_name;
												#my @Temp_gene_type=sort keys %{$temp_hash{"gene"}{$temp_gene_name}};
												#my $temp_gene_type=join ",", @Temp_gene_type;
												$final_gene_type.="/".$hash_gene_info{"gene_type"}{$temp_gene_name};
												$final_gene_strand.="/".$hash_gene_info{"strand"}{$temp_gene_name};
            my @Temp_extra_base=sort keys %{$temp_hash{"extra_base"}{$temp_gene_name}};
            my $temp_extra_base=join ";", @Temp_extra_base;
            $final_extra_base.="/".$temp_extra_base;
												$final_fusion_position_location.="/"."splice_site";
            
            my @Temp_fusion_region_type=sort keys %{$temp_hash{"position"}{"splice_site"}{$temp_gene_name}};
            my $temp_fusion_region_type=join ",", @Temp_fusion_region_type;
            $final_fusion_region_type.="/".$temp_fusion_region_type;
            
            if(exists $temp_hash{"position"}{"exon"}{$temp_gene_name}){
                $final_fusion_position_location.=";"."exon";
                @Temp_fusion_region_type=sort keys %{$temp_hash{"position"}{"exon"}{$temp_gene_name}};
                $temp_fusion_region_type=join ",", @Temp_fusion_region_type;
                $final_fusion_region_type.=";".$temp_fusion_region_type;
            }
								}
				}
				if(exists $temp_hash{"position"}{"exon"}){
								foreach my $temp_gene_name (sort keys %{$temp_hash{"position"}{"exon"}}){
												next if exists $temp_hash{"position"}{"splice_site"}{$temp_gene_name};
												$final_gene_name.="/".$temp_gene_name;
												#my @Temp_gene_type=sort keys %{$temp_hash{"gene"}{$temp_gene_name}};
												#my $temp_gene_type=join ",", @Temp_gene_type;
												$final_gene_type.="/".$hash_gene_info{"gene_type"}{$temp_gene_name};
												$final_gene_strand.="/".$hash_gene_info{"strand"}{$temp_gene_name};
            my @Temp_extra_base=sort keys %{$temp_hash{"extra_base"}{$temp_gene_name}};
            my $temp_extra_base=join ";", @Temp_extra_base;
            $final_extra_base.="/".$temp_extra_base;
												$final_fusion_position_location.="/"."exon";
            my @Temp_fusion_region_type=sort keys %{$temp_hash{"position"}{"exon"}{$temp_gene_name}};
            my $temp_fusion_region_type=join ",", @Temp_fusion_region_type;
            $final_fusion_region_type.="/".$temp_fusion_region_type;
								}
				}
				
				if(! exists $temp_hash{"position"}{"splice_site"} && ! exists $temp_hash{"position"}{"exon"} && exists $temp_hash{"position"}{"intron"}){
								foreach my $temp_gene_name (sort keys %{$temp_hash{"position"}{"intron"}}){
												$final_gene_name.="/".$temp_gene_name;
												#my @Temp_gene_type=sort keys %{$temp_hash{"gene"}{$temp_gene_name}};
												#my $temp_gene_type=join ",", @Temp_gene_type;
												$final_gene_type.="/".$hash_gene_info{"gene_type"}{$temp_gene_name};
												$final_gene_strand.="/".$hash_gene_info{"strand"}{$temp_gene_name};
            my @Temp_extra_base=sort keys %{$temp_hash{"extra_base"}{$temp_gene_name}};
            my $temp_extra_base=join ";", @Temp_extra_base;
            $final_extra_base.="/".$temp_extra_base;
												$final_fusion_position_location.="/"."intron";
            $final_fusion_region_type.="/"."non_coding";
								}
				}
				
				if($final_gene_name eq "NA"){
								$final_gene_name="Unannotated";
								$final_gene_type="Unknown";
								$final_gene_strand=".";
								$final_fusion_position_location="Unknown";
        $final_fusion_region_type="Unknown";
        $final_extra_base="N";
				}else{
								$final_gene_name=~s/NA\///;
								$final_gene_type=~s/NA\///;
								$final_gene_strand=~s/NA\///;
								$final_fusion_position_location=~s/NA\///;
        $final_fusion_region_type=~s/NA\///;
        $final_extra_base=~s/NA\///;
				}
				
				return ($final_gene_name, $final_gene_type, $final_gene_strand, $final_fusion_position_location, $final_fusion_region_type, $final_extra_base);
}

sub getFlankingSequence{
    my ($chr, $pos, $strand, $length_1_fusion, $length_2_fusion, $upstream_or_downstream_segment)=@_;
    my ($start, $end);
    my $append_N_length=0; # if $start < 1, append N
    if($strand eq '+'){
        $start=$pos-$length_1_fusion+1;
        if($start<1){
            $append_N_length=abs($start)+1;
            $start=1;
        }
        $end=$pos+$length_2_fusion; # end bigger than chr size will generate short sequence
    }else{
        $start=$pos-$length_2_fusion;
        if($start<1){
            $append_N_length=abs($start)+1;
            $start=1;
        }
        $end=$pos+$length_1_fusion-1; # end bigger than chr size will generate short sequence
    }
    
    my $flankingSeq;
    my $artifact_coordinate=$chr.":".$start."-".$end;
    if($upstream_or_downstream_segment eq "upstream"){
								if($strand eq '+'){
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}'`;
								}else{
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}' | rev | tr 'ACGT' 'TGCA'`;
								}
    }else{
								if($strand eq '+'){
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}' | rev | tr 'ACGT' 'TGCA'`;
								}else{
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}'`;
								}
    }
    chomp($flankingSeq);
    
    my $length_artifact_fa=$length_1_fusion+$length_2_fusion;
    if($upstream_or_downstream_segment eq "upstream"){
								if($strand eq '+'){
												if($append_N_length>0){
																my $temp_N_string="N" x $append_N_length;
																$flankingSeq=$temp_N_string.$flankingSeq;
												}
												my $length_flankingSeq=length($flankingSeq);
												if($length_flankingSeq<$length_artifact_fa){
																my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
																$flankingSeq=$flankingSeq.$temp_N_string;
												}
								}else{
												if($append_N_length>0){
																my $temp_N_string="N" x $append_N_length;
																$flankingSeq=$flankingSeq.$temp_N_string;
												}
												my $length_flankingSeq=length($flankingSeq);
												if($length_flankingSeq<$length_artifact_fa){
																my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
																$flankingSeq=$temp_N_string.$flankingSeq;
												}
								}
    }else{
								if($strand eq '+'){
												if($append_N_length>0){
																my $temp_N_string="N" x $append_N_length;
																$flankingSeq=$flankingSeq.$temp_N_string;
												}
												my $length_flankingSeq=length($flankingSeq);
												if($length_flankingSeq<$length_artifact_fa){
																my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
																$flankingSeq=$temp_N_string.$flankingSeq;
												}
								}else{
												if($append_N_length>0){
																my $temp_N_string="N" x $append_N_length;
																$flankingSeq=$temp_N_string.$flankingSeq;
												}
												my $length_flankingSeq=length($flankingSeq);
												if($length_flankingSeq<$length_artifact_fa){
																my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
																$flankingSeq=$flankingSeq.$temp_N_string;
												}
								}
    }
    
    return($flankingSeq);
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
This script was used to annotate fusions
Usage: perl $scriptName -a reference.fa annotation.gpe fusions.tsv >output
Options:

    -f --filter_in_the_same_gene    STR     Filter fusions in the same gene [default: $filter_in_the_same_gene]
    -a --fasta    		                STR	    Reference genome fasta file (must be given)
    -h --help     		                        Print this help information
HELP
    exit(-1);
}


