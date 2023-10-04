#!/usr/bin/env perl

# 2023-03-08

# 1. Function 
# scan for canonical split motif in region defined by --motif_searching_length


# 2. Input
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
# column 13: presence of canonical splice site
# column 14: split reads (e.g., read_23)
# column 15: read pairs (e.g., read_38,read_62,read_70,read_8)
# column 16: distance of read pair to breakpoint 1 
# column 17: distance of read pair to breakpoint 2
# column 18: discordant read cluster IDs around breakpoint 1 (e.g., 1)
# column 19: discordant read cluster IDs around breakpoint 2 (e.g., 1)
# column 20: supporting reads count in each cluster around breakpoint 1 (e.g., 4)
# column 21: supporting reads count in each cluster around breakpoint 2 (e.g., 3)



use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $motif_searching_length=5;
my $fasta;
my $filter_by_canonical_split_motif="Y";

GetOptions(
    'm|motif_searching_length=i'												=> \$motif_searching_length,
    'f|fasta=s'     																								=> \$fasta,
    'c|filter_by_canonical_split_motif=s'  	=> \$filter_by_canonical_split_motif,
    'h|help'    																												=> sub{usage()}
)||usage();

die "Please set -f or --fasta" if ! defined $fasta;

say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2", 
    "Split_read_count", "Read_pair_count", 
    "Minimum_read_pair_distance_to_breakpoint_1", "Minimum_read_pair_distance_to_breakpoint_2",
    "SD_(discordant_reads)%_in_breakpoint_1", "SD_(discordant_reads)%_in_breakpoint_2",
    "Have_canonical_motif",
    "Split_reads", "Read_pairs",
    "Read_pair_distance_to_breakpoint_1", "Read_pair_distance_to_breakpoint_2",
    "Cluster_ids_breakpoint_1", "Cluster_ids_breakpoint_2",
    "Read_count_in_each_cluster_breakpoint_1", "Read_count_in_each_cluster_breakpoint_2");


open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2,
        $split_read_count, $read_pair_count,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
        $split_reads, $read_pairs,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2,
        $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
        $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2)=split "\t", $_;
        
    # is contain canonical split motif
    my $flanking_search_sequence_breakpoint_1=&getFlankingSequence($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $motif_searching_length, $motif_searching_length, "upstream");
    my $is_breakpoint_1_contain_donor_motif=($flanking_search_sequence_breakpoint_1=~/GT/) ? "Y" : "N";
    $flanking_search_sequence_breakpoint_1=~tr/[ACTG]/[TGAC]/;
    $flanking_search_sequence_breakpoint_1=reverse($flanking_search_sequence_breakpoint_1);
    my $is_breakpoint_1_contain_acceptor_motif=($flanking_search_sequence_breakpoint_1=~/[CTA]AG/) ? "Y" : "N";
    
    my $flanking_search_sequence_breakpoint_2=&getFlankingSequence($chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $motif_searching_length, $motif_searching_length, "downstream");
    my $is_breakpoint_2_contain_acceptor_motif=($flanking_search_sequence_breakpoint_2=~/[CTA]AG/) ? "Y" : "N";
    $flanking_search_sequence_breakpoint_2=~tr/[ACTG]/[TGAC]/;
    $flanking_search_sequence_breakpoint_2=reverse($flanking_search_sequence_breakpoint_2);
    my $is_breakpoint_2_contain_donor_motif=($flanking_search_sequence_breakpoint_2=~/GT/) ? "Y" : "N";
    
    my $is_split_site_contain_canonical_motif="N";
    $is_split_site_contain_canonical_motif="Y" if ($is_breakpoint_1_contain_donor_motif eq "Y" && $is_breakpoint_2_contain_acceptor_motif eq "Y") ||
								($is_breakpoint_2_contain_donor_motif eq "Y" && $is_breakpoint_1_contain_acceptor_motif eq "Y");
    
    next if $filter_by_canonical_split_motif eq "Y" && $is_split_site_contain_canonical_motif ne "Y";
    
    say join "\t", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2,
        $split_read_count, $read_pair_count,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2,
        $sd_percentage_discorant_read_in_each_cluster_breakpoint_1, $sd_percentage_discorant_read_in_each_cluster_breakpoint_2,
        $is_split_site_contain_canonical_motif,
        $split_reads, $read_pairs,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2,
        $cluster_ids_breakpoint_1, $cluster_ids_breakpoint_2,
        $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2);
}
close(IN);


sub getFlankingSequence{
    my ($chr, $pos, $strand, $outside_length_fusion, $inside_length_fusion, $upstream_or_downstream_segment)=@_;
    my ($start, $end);
    my $append_N_length=0; # if $start < 1, append N
    if($strand eq '+'){
        $start=$pos-$outside_length_fusion+1;
        if($start<1){
            $append_N_length=abs($start)+1;
            $start=1;
        }
        $end=$pos+$inside_length_fusion; # end bigger than chr size will generate short sequence
    }else{
        $start=$pos-$inside_length_fusion;
        if($start<1){
            $append_N_length=abs($start)+1;
            $start=1;
        }
        $end=$pos+$outside_length_fusion-1; # end bigger than chr size will generate short sequence
    }
    
    my $flankingSeq;
    my $artifact_coordinate=$chr.":".$start."-".$end;
    if($upstream_or_downstream_segment eq "upstream"){
								if($strand eq '+'){
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print toupper(output);}'`;
								}else{
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print toupper(output);}' | rev | tr 'ACGT' 'TGCA'`;
								}
    }else{
								if($strand eq '+'){
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print toupper(output);}' | rev | tr 'ACGT' 'TGCA'`;
								}else{
												$flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print toupper(output);}'`;
								}
    }
    chomp($flankingSeq);
    
    my $length_artifact_fa=$outside_length_fusion+$inside_length_fusion;
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
This script was used to get blat supported statistics
Usage: perl $scriptName input >output
Options:

    -m --motif_searching_length											INT				Splice site motifs (GT in the donor, AAG/CAG/TAG in the acceptor) are searched within this window size of breakpoints [default: $motif_searching_length]
    -f --fasta    																								STR				Reference genome fasta file (must be given)
    -c --filter_by_canonical_split_motif		STR				Filter by canonical splice site motif [default: $filter_by_canonical_split_motif]
    -h --help     																															Print this help information
HELP
    exit(-1);
}



