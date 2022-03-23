#!/usr/bin/env perl

# 2022-03-22

# 1. Function
# Annotate fusions
# Filter fusion locate in the same gene
# Swith left and right segment so fusion is in sense-sense orientation
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

# delete chr*_ by $2!~/_/
system("awk 'BEGIN{FS=OFS=\"\t\"} NR>1 && \$2!~/_/{print \$2,\$4,\$5,\$1,\$3;}' $ARGV[0] | sort -k1,1 -k2,2n | uniq >temp_annotation.bed");
system("awk 'BEGIN{FS=OFS=\"\t\"} NR>1{print \$1,\$2-1,\$2; print \$4,\$5-1,\$5;}' $ARGV[1] | sort -k1,1 -k2,2n | uniq >temp_fusion_breakpoint.bed");
system("bedtools intersect -wa -wb -a temp_fusion_breakpoint.bed -b temp_annotation.bed | cut -f1,3- >temp_overlapped_gene.tsv");

open IN, "temp_overlapped_gene.tsv" or die "Can't open temp_overlapped_gene.tsv:$!";
my %breakpoint2gene;
my %select_transcripts;
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
    my $temp_gene_info=join ";", ($chr, $transcript_start, $transcript_end, $transcript_id, $strand);
    next if ! exists $select_transcripts{$temp_gene_info};
    $select_transcripts{$temp_gene_info}=$_;
}
close(IN);

say join "\t", ("Chr_1", "Breakpoint_1", "Strand_1", "Chr_2", "Breakpoint_2", "Strand_2",
     "Split_read_count_(Tophat_and_Blat)", "Read_pair_count_(Tophat)",
     "Minimum_read_distance_to_left", "Minimum_read_distance_to_right",
     "Identity", "Minimum_blat_distance_to_1", "Minimum_blat_distance_to_2",
     "Total_clusters_1", "Total_clusters_2", "Total_clusters_(merge)",
     "(discordant_reads)%_support_fusion",
     "SD_(discordant_reads)%_in_clusters_1", "SD_(discordant_reads)%_in_clusters_2",
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
        $split_reads_blat, $read_pairs_tophat,
	$cluster_ids_left, $cluster_ids_right,
        $total_reads_in_each_cluster_left, $total_reads_in_each_cluster_right)=split "\t", $_;
    
    my (@Transcripts_left, @Transcripts_right);
    my $temp_key=join ";", ($chr_left, $pos_left);
    if(exists $breakpoint2gene{$temp_key}){
	my $temp_value=&selectLongestExonTranscript($chr_left, $pos_left);
	@Transcripts_left=@$temp_value;
    }else{
	push @Transcripts_left, join "\t",("NA") x 15; # at least one for loop
    }
    $temp_key=join ";", ($chr_right, $pos_right);
    if(exists $breakpoint2gene{$temp_key}){
	my $temp_value=&selectLongestExonTranscript($chr_right, $pos_right);
	@Transcripts_right=@$temp_value;
    }else{
	push @Transcripts_right, join "\t",("NA") x 15; # at least one for loop
    }
    
    
    my (@sense_fusions, @sense_fusions_switch_segment, @fusion_to_unannotated, @fusion_to_unannotated_switch_segment); # @sense_fusions keep current segment order, @sense_fusions_switch_segment swith left and right fusion segment
    my $is_in_the_same_gene=0;
    foreach my $transcript_left (@Transcripts_left){
	my ($transcript_id_left, $chr_transcript_left, $strand_transcript_left, $transcript_start_left, $transcript_end_left, $cds_start_left, $cds_end_left,
	    $exon_count_left, $exon_starts_left, $exon_ends_left, $score_left, $symbol_left, $cds_start_stat_left, $cds_end_stat_left, $exon_frames_left)=split "\t", $transcript_left;
	foreach my $transcript_right (@Transcripts_right){
	    my ($transcript_id_right, $chr_transcript_right, $strand_transcript_right, $transcript_start_right, $transcript_end_right, $cds_start_right, $cds_end_right,
		$exon_count_right, $exon_starts_right, $exon_ends_right, $score_right, $symbol_right, $cds_start_stat_right, $cds_end_stat_right, $exon_frames_right)=split "\t", $transcript_right;
	    
	    if($strand_left eq '+' && $strand_right eq '-' && $symbol_left eq $symbol_right && $symbol_left ne "NA"){
		$is_in_the_same_gene=1;
		goto OUTPUT;
	    }
	    
	    # 1. if $strand_left ne $strand_right, $strand_transcript_left eq $strand_transcript_right will be a valid fusion
	    # 2. if $strand_left eq $strand_right, $strand_transcript_left ne $strand_transcript_right will be a valid fusion
	    # 3. if $strand_left ne $strand_transcript_left, switch left and right symbol
	    # symbol may contain "@", "." and "-"
	    next if $symbol_left eq "NA" && $symbol_right eq "NA";
	    if($symbol_left ne "NA" && $symbol_right ne "NA"){
		my ($gene_type_left, $fusion_pos_location_left, $fusion_region_type_left, $frame_offset_left)=&compareBreakpointAndTranscript($pos_left, $transcript_left);
		my ($gene_type_right, $fusion_pos_location_right, $fusion_region_type_right, $frame_offset_right)=&compareBreakpointAndTranscript($pos_right, $transcript_right);
		
		if((($strand_left ne $strand_right) && ($strand_transcript_left eq $strand_transcript_right)) ||
		    (($strand_left eq $strand_right) && ($strand_transcript_left ne $strand_transcript_right))){ # valid fusion
		    if($strand_left ne $strand_transcript_left){ # swith left and right symbol
			my $fusion_symbols=$symbol_right."--".$symbol_left;
			my $fusion_transcript_ids=$transcript_id_right."--".$transcript_id_left;
			my $fusion_gene_types=$gene_type_right."--".$gene_type_left;
			my $fusion_pos_locations=$fusion_pos_location_right."--".$fusion_pos_location_left;
			my $fusion_region_types=$fusion_region_type_right."--".$fusion_region_type_left;
			if($fusion_region_types eq "cds--cds"){
			    if($frame_offset_left==$frame_offset_right){
				$fusion_region_types.="_(in-frame)";
			    }else{
				$fusion_region_types.="_(frame-shift)";
			    }
			}
			push @sense_fusions_switch_segment, join ",",($fusion_symbols,$fusion_transcript_ids,$fusion_gene_types,$fusion_pos_locations,$fusion_region_types);
		    }else{
			my $fusion_symbols=$symbol_left."--".$symbol_right;
			my $fusion_transcript_ids=$transcript_id_left."--".$transcript_id_right;
			my $fusion_gene_types=$gene_type_left."--".$gene_type_right;
			my $fusion_pos_locations=$fusion_pos_location_left."--".$fusion_pos_location_right;
			my $fusion_region_types=$fusion_region_type_left."--".$fusion_region_type_right;
			if($fusion_region_types eq "cds--cds"){
			    if($frame_offset_left==$frame_offset_right){
				$fusion_region_types.="_(in-frame)";
			    }else{
				$fusion_region_types.="_(frame-shift)";
			    }
			}
			push @sense_fusions, join ",",($fusion_symbols,$fusion_transcript_ids,$fusion_gene_types,$fusion_pos_locations,$fusion_region_types);
		    }
		}else{
		    my $fusion_symbols=$symbol_left."--"."unannotated";
		    my $fusion_transcript_ids=$transcript_id_left."--"."unannotated";
		    my $fusion_gene_types=$gene_type_left."--"."unknown";
		    my $fusion_pos_locations=$fusion_pos_location_left."--"."intergenic";
		    my $fusion_region_types=$fusion_region_type_left."--"."intergenic";
		    push @fusion_to_unannotated, join ",",($fusion_symbols,$fusion_transcript_ids,$fusion_gene_types,$fusion_pos_locations,$fusion_region_types);
		    
		    $fusion_symbols=$symbol_right."--"."unannotated";
		    $fusion_transcript_ids=$transcript_id_right."--"."unannotated";
		    $fusion_gene_types=$gene_type_right."--"."unknown";
		    $fusion_pos_locations=$fusion_pos_location_right."--"."intergenic";
		    $fusion_region_types=$fusion_region_type_right."--"."intergenic";
		    push @fusion_to_unannotated_switch_segment, join ",",($fusion_symbols,$fusion_transcript_ids,$fusion_gene_types,$fusion_pos_locations,$fusion_region_types);
		}
	    }else{
		if($symbol_left ne "NA"){
		    my ($gene_type_left, $fusion_pos_location_left, $fusion_region_type_left, $frame_offset_left)=&compareBreakpointAndTranscript($pos_left, $transcript_left);
		    my $fusion_symbols=$symbol_left."--"."unannotated";
		    my $fusion_transcript_ids=$transcript_id_left."--"."unannotated";
		    my $fusion_gene_types=$gene_type_left."--"."unknown";
		    my $fusion_pos_locations=$fusion_pos_location_left."--"."intergenic";
		    my $fusion_region_types=$fusion_region_type_left."--"."intergenic";
		    push @fusion_to_unannotated, join ",",($fusion_symbols,$fusion_transcript_ids,$fusion_gene_types,$fusion_pos_locations,$fusion_region_types);
		}else{
		    my ($gene_type_right, $fusion_pos_location_right, $fusion_region_type_right, $frame_offset_right)=&compareBreakpointAndTranscript($pos_right, $transcript_right);
		    my $fusion_symbols=$symbol_right."--"."unannotated";
		    my $fusion_transcript_ids=$transcript_id_right."--"."unannotated";
		    my $fusion_gene_types=$gene_type_right."--"."unknown";
		    my $fusion_pos_locations=$fusion_pos_location_right."--"."intergenic";
		    my $fusion_region_types=$fusion_region_type_right."--"."intergenic";
		    push @fusion_to_unannotated_switch_segment, join ",",($fusion_symbols,$fusion_transcript_ids,$fusion_gene_types,$fusion_pos_locations,$fusion_region_types);
		}
	    }
	}
    }
    
    my $fusion_annotation="NA";
    my $sense_fusion_count=@sense_fusions;
    my $sense_fusion_switch_count=@sense_fusions_switch_segment;
    my $fusion_to_unannotated_count=@fusion_to_unannotated;
    my $fusion_to_unannotated_swith_count=@fusion_to_unannotated_switch_segment;
    
    if($sense_fusion_count+$sense_fusion_switch_count+$fusion_to_unannotated_count+$fusion_to_unannotated_swith_count==0){
	$fusion_annotation="unannotated--unannotated,unannotated--unannotated,unknown--unknown,intergenic--intergenic,intergenic--intergenic";
    }elsif($sense_fusion_count+$sense_fusion_switch_count==0){
	if($fusion_to_unannotated_count>0 && $fusion_to_unannotated_swith_count>0){
	    foreach my $temp_value (@fusion_to_unannotated){
		$fusion_annotation.=";".$temp_value;
	    }
	    foreach my $temp_value (@fusion_to_unannotated_switch_segment){
		$fusion_annotation.=";".$chr_right.":".$pos_right.":".$strand_right."--".$chr_left.":".$pos_left.":".$strand_left.",".$temp_value;
	    }
	}elsif($fusion_to_unannotated_swith_count==0){
	    foreach my $temp_value (@fusion_to_unannotated){
		$fusion_annotation.=";".$temp_value;
	    }
	}else{ # $fusion_to_unannotated_count==0
	    foreach my $temp_value (@fusion_to_unannotated_switch_segment){
		$fusion_annotation.=";".$temp_value;
	    }
	    ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right,
	    $min_read_pair_distance_left, $min_read_pair_distance_right,
	    $minimum_blat_distance_left, $minimum_blat_distance_right,
	    $total_clusters_left, $total_clusters_right,
	    $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right)=(
		$chr_right, $pos_right, $strand_right, $chr_left, $pos_left, $strand_left,
		$min_read_pair_distance_right, $min_read_pair_distance_left,
		$minimum_blat_distance_right, $minimum_blat_distance_left,
		$total_clusters_right, $total_clusters_left,
		$sd_percentage_discorant_read_in_each_cluster_right, $sd_percentage_discorant_read_in_each_cluster_left);
	}
	$fusion_annotation=~s/NA;//;
    }else{ # $sense_fusion_count+$sense_fusion_switch_count>0
	if($sense_fusion_count>0 && $sense_fusion_switch_count>0){
	    foreach my $temp_value (@sense_fusions){
		$fusion_annotation.=";".$temp_value;
	    }
	    foreach my $temp_value (@sense_fusions_switch_segment){
		$fusion_annotation.=";".$chr_right.":".$pos_right.":".$strand_right."--".$chr_left.":".$pos_left.":".$strand_left.",".$temp_value;
	    }
	}elsif($sense_fusion_switch_count==0){
	    foreach my $temp_value (@sense_fusions){
		$fusion_annotation.=";".$temp_value;
	    }
	}else{ # $sense_fusion_count==0
	    foreach my $temp_value (@sense_fusions_switch_segment){
		$fusion_annotation.=";".$temp_value;
	    }
	    ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right,
	    $min_read_pair_distance_left, $min_read_pair_distance_right,
	    $minimum_blat_distance_left, $minimum_blat_distance_right,
	    $total_clusters_left, $total_clusters_right,
	    $sd_percentage_discorant_read_in_each_cluster_left, $sd_percentage_discorant_read_in_each_cluster_right)=(
		$chr_right, $pos_right, $strand_right, $chr_left, $pos_left, $strand_left,
		$min_read_pair_distance_right, $min_read_pair_distance_left,
		$minimum_blat_distance_right, $minimum_blat_distance_left,
		$total_clusters_right, $total_clusters_left,
		$sd_percentage_discorant_read_in_each_cluster_right, $sd_percentage_discorant_read_in_each_cluster_left);
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

sub selectLongestExonTranscript{
    my ($chr, $pos)=@_;
    my $temp_key=join ";", ($chr, $pos);
    my @genes_info=sort keys %{$breakpoint2gene{$temp_key}};
    
    my %temp_select_transcripts;
    foreach my $gene_info (@genes_info){
	my ($transcript_id, $chr, $strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames)=split "\t", $select_transcripts{$gene_info};
	my @Exon_starts=split ",", $exon_starts;
	my @Exon_ends=split ",", $exon_ends;
	my $exon_length=0;
	for(my $i=0; $i<$exon_count; $i++){
	    $exon_length+=$Exon_ends[$i]-$Exon_starts[$i];
	}
	if(exists $temp_select_transcripts{$name_2}){
	    if($exon_length>$temp_select_transcripts{$name_2}{exon_length}){
		$temp_select_transcripts{$name_2}{exon_length}=$exon_length;
		$temp_select_transcripts{$name_2}{transcript}=join "\t", ($transcript_id, $chr, $strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames);
	    }
	}else{
	    $temp_select_transcripts{$name_2}{exon_length}=$exon_length;
	    $temp_select_transcripts{$name_2}{transcript}=join "\t", ($transcript_id, $chr, $strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames);
	}
    }
    
    my @Select_transcripts;
    foreach my $gene (sort keys %temp_select_transcripts){
	push @Select_transcripts, $temp_select_transcripts{$gene}{transcript};
    }
    return(\@Select_transcripts);
}

sub compareBreakpointAndTranscript{
    my ($fusion_pos, $transcript)=@_;
    my ($transcript_id, $chr, $strand, $transcript_start, $transcript_end, $cds_start, $cds_end, $exon_count, $exon_starts, $exon_ends, $score, $name_2, $cds_start_stat, $cds_end_stat, $exon_frames)=split "\t", $transcript;
    my @Exon_starts=split ",", $exon_starts;
    my @Exon_ends=split ",", $exon_ends;
    my @Exon_frames=split ",", $exon_frames;
    my $gene_type=($cds_start==$cds_end) ? "non_coding_gene" : "protein_coding_gene";
    
    my $fusion_region_type="non_coding";
    if($gene_type eq "protein_coding_gene"){
	if($fusion_pos>=$cds_start+1 && $fusion_pos<=$cds_end){
	    $fusion_region_type="cds";
	}else{
	    if($fusion_pos<$cds_start+1){
		$fusion_region_type=($strand eq '+') ? "5'UTR" : "3'UTR";
	    }else{
		$fusion_region_type=($strand eq '+') ? "3'UTR" : "5'UTR";
	    }
	}
    }
    
    my $fusion_position_location;
    my $fusion_location_index;
    for(my $i=0; $i<$exon_count; $i++){
	last if $Exon_starts[$i]+1>$fusion_pos;
	if($i<$exon_count-1){
	    if($fusion_pos>$Exon_ends[$i] && $fusion_pos<$Exon_starts[$i+1]+1){ # in intron
		$fusion_position_location="intron";
		$fusion_location_index=$i;
		last;
	    }
	}
	next if $fusion_pos>$Exon_ends[$i];
	$fusion_location_index=$i;
	if($fusion_pos==$Exon_starts[$i]+1){
	    $fusion_position_location="split_site";
	}elsif($fusion_pos==$Exon_ends[$i]){
	    $fusion_position_location="split_site";
	}else{
	    $fusion_position_location="exon";
	}
    }
    
    my $frame_shift="N";
    if($gene_type eq "protein_coding_gene"){
	if($fusion_region_type eq "cds"){
	    if($strand eq '+'){ # gene strand
		if($fusion_pos>=$Exon_ends[$fusion_location_index]){
		    $frame_shift=$Exon_frames[$fusion_location_index+1];
		}elsif($fusion_pos==$Exon_starts[$fusion_location_index]+1){
		    $frame_shift=$Exon_frames[$fusion_location_index];
		}else{
		    my $temp_distance=$fusion_pos-$Exon_starts[$fusion_location_index];
		    $frame_shift=($temp_distance+$Exon_frames[$fusion_location_index])%3;
		}
	    }else{
		if($fusion_pos>=$Exon_ends[$fusion_location_index]){
		    $frame_shift=$Exon_frames[$fusion_location_index];
		}elsif($fusion_pos==$Exon_starts[$fusion_location_index]+1){
		    $frame_shift=$Exon_frames[$fusion_location_index-1];
		}else{
		    my $temp_distance=$Exon_ends[$fusion_location_index]-$fusion_pos+1;
		    $frame_shift=($temp_distance+$Exon_frames[$fusion_location_index])%3;
		}
	    }
	}
    }
    
    return ($gene_type, $fusion_position_location, $fusion_region_type, $frame_shift);
}

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


