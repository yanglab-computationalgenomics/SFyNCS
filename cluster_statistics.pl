#!/usr/bin/env perl

# 2022-11-04

# 1. Function
# Get totoal discordant reads, total clusters, (discordant reads)% in each cluster, (discordant reads)% in fusion cluster, (discordant reads)% support ratio within --flanking_length bp
# Drop cluster id

# to do in next version:
# Change chr23 to chrX, and chr24 to chrY

# 2. Input 
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand 
# column 4: chromosome of breakpoint 2
# column 5: breakpoint 2 position (1-based)
# column 6: breakpoint 2 strand
# column 7: cluster id
# column 8: split read count reported by STAR (e.g., 1)
# column 9: read pair count reported by STAR (e.g., 4)
# column 10: split read count supported by TopHat (e.g., 1)
# column 11: potential split read count not supported by TopHat (e.g., 0)
# column 12: read pair count reported by TopHat
# column 13: split read count supported by Blat (e.g., 1)
# column 14: split read count supported by Blat (considering split read supported by TopHat only, e.g., 1)
# column 15: minimal distance between read pair and breakpoint 1 after aligning by TopHat (e.g., 13)
# column 16: minimal distance between read pair and breakpoint 2 after aligning by TopHat (e.g., 22)
# column 17: sequence identity (e.g., 0.52)
# column 18: minimal distance between split read and breakpoint 1 after aligning by Blat (e.g., 0, use mean if both read 1 and read 2 are split read)
# column 19: minimal distance between split read and breakpoint 2 after aligning by Blat (e.g., 0, use mean if both read 1 and read 2 are split read)
# column 20: minimal distance between split read and breakpoint 1 after aligning by Blat (considering split read supported by TopHat only, e.g., 0, use mean if both read 1 and read 2 are split read)
# column 21: minimal distance between split read and breakpoint 2 after aligning by Blat (considering split read supported by TopHat only, e.g., 0, use mean if both read 1 and read 2 are split read)
# column 22: presence of canonical splice site
# column 23: split reads reported by STAR (e.g., read_23)
# column 24: read pairs reported by STAR (e.g., read_38,read_62,read_70,read_8)
# column 25: split reads supported by TopHat (e.g., read_23)
# column 26: potential split reads not supported by TopHat (e.g., NA)
# column 27: read pairs supported by TopHat (e.g., read_38,read_70,read_8)
# column 28: split reads supported by Blat(blat tophat split read and tophat potential split read)
# column 29: distance between each read pair and breakpoint 1 after aligning by TopHat (e.g., 740,23,13)
# column 30: distance between each read pair and breakpoint 2 after aligning by TopHat (e.g., 292,2592,22)
# column 31: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 32: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 33: distance between each split read and breakpoint 1 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 34: distance between each split read and breakpoint 2 after aligning to artifact reference by Blat (considering split read supported by TopHat only, use mean if both read 1 and read 2 are split read)
# column 35: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 1 (e.g., --AGTGGGCCAGGTAG-GGCTGG)
# column 36: Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 2 (e.g., CCACT--GCCAGG-AGAACCTCA)

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
    "Split_read_count_(star)", "Read_pair_count_(star)", "Split_read_count_(tophat)", "Potential_split_read_count_(tophat)", "Read_pair_count_(tophat)",
    "Split_read_count_(blat_tophat_split_and_tophat_potential_split_reads)", "Split_read_count_(blat_tophat_split_reads)",
    "Minimum_read_pair_distance_to_breakpoint_1", "Minimum_read_pair_distance_to_breakpoint_2",
    "Sequence_identity",
    "Minimum_blat_distance_to_breakpoint_1_(tophat_split_and_potential_split_reads)", "Minimum_blat_distance_to_breakpoint_2_(tophat_split_and_potential_split_reads)",
    "Minimum_blat_distance_to_breakpoint_1_(tophat_split_reads)", "Minimum_blat_distance_to_breakpoint_2_(tophat_split_reads)",
    "Have_canonical_motif",
    "Cluster_count_breakpoint_1", "Cluster_count_breakpoint_2", "Cluster_count_breakpoint_1_and_2",
    "(discordant_reads)%_support_fusion",
    "SD_(discordant_reads)%_in_breakpoint_1", "SD_(discordant_reads)%_in_breakpoint_2",
    "Split_reads_(star)", "Read_pairs_(star)",
    "Split_reads_(tophat)", "Potential_split_reads_(tophat)", "Read_pairs_(tophat)",
    "Split_reads_(blat_tophat_split_and_tophat_potential_split_reads)",
    "Read_pair_distance_to_breakpoint_1", "Read_pair_distance_to_breakpoint_2",
    "Blat_distance_to_breakpoint_1_(tophat_split_and_tophat_potential_split_reads)", "Blat_distance_to_breakpoint_2_(tophat_split_and_tophat_potential_split_reads)",
    "Blat_distance_to_breakpoint_1_(tophat_split_reads)", "Blat_distance_to_breakpoint_2_(tophat_split_eads)",
    "Sequence_identity_alignment_breakpoint_1", "Sequence_identity_alignment_breakpoint_2",
    "Cluster_ids_breakpoint_1", "Cluster_ids_breakpoint_2",
    "Read_count_in_each_cluster_breakpoint_1", "Read_count_in_each_cluster_breakpoint_2");

open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cluster_id,
        $split_read_count_star, $read_pair_count_star, $split_read_count_tophat, $potential_split_read_count_tophat, $read_pair_count_tophat, $split_read_count_blat, $split_read_count_blat_base_on_tophat_split_read,
        $min_read_pair_distance_breakpoint_1, $min_read_pair_distance_breakpoint_2, $identity_output,
        $minimum_blat_distance_breakpoint_1, $minimum_blat_distance_breakpoint_2, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_1, $minimum_blat_distance_base_on_tophat_split_read_breakpoint_2,
        $is_split_site_contain_canonical_motif,
        $split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat, $split_reads_blat,
        $read_pair_distance_to_breakpoint_1, $read_pair_distance_to_breakpoint_2, 
        $blat_distance_breakpoint_1, $blat_distance_breakpoint_2, $blat_distance_base_on_tophat_split_read_breakpoint_1, $blat_distance_base_on_tophat_split_read_breakpoint_2,
        $identity_seq_breakpoint_1, $identity_seq_breakpoint_2)=split "\t", $_;
    my ($total_clusters_breakpoint_1, $percentage_discordant_read_in_each_cluster_breakpoint_1, $all_discordant_reads_breakpoint_1, $clusters_breakpoint_1, $cluster_ids_breakpoint_1, $total_reads_in_each_cluster_breakpoint_1)=&get_cluster_and_discordant_statistics($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cluster_id);
    my ($total_clusters_breakpoint_2, $percentage_discordant_read_in_each_cluster_breakpoint_2, $all_discordant_reads_breakpoint_2, $clusters_breakpoint_2, $cluster_ids_breakpoint_2, $total_reads_in_each_cluster_breakpoint_2)=&get_cluster_and_discordant_statistics($chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cluster_id);
    
    my $key_breakpoint_1=join ":", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1);
    my $key_breakpoint_2=join ":", ($chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
    my %processed_reads;
    my ($total_discordant_reads, $fusion_supported_reads)=(0) x 2;
    foreach my $read (split ",", $all_discordant_reads_breakpoint_1){
        next if exists $processed_reads{$read};
        $processed_reads{$read}=0;
        $total_discordant_reads++;
        $fusion_supported_reads++ if (exists $overlapped_discordant_reads{$key_breakpoint_1}{"split_reads"}{$read} && exists $overlapped_discordant_reads{$key_breakpoint_2}{"split_reads"}{$read}) ||
            (exists $overlapped_discordant_reads{$key_breakpoint_1}{"read_pairs"}{$read} && exists $overlapped_discordant_reads{$key_breakpoint_2}{"read_pairs"}{$read});
    }
    foreach my $read (split ",", $all_discordant_reads_breakpoint_2){
        next if exists $processed_reads{$read};
        $processed_reads{$read}=0;
        $total_discordant_reads++;
    }
    my $percentage_support_fusion=sprintf("%.2f", $fusion_supported_reads/$total_discordant_reads);
    
    my %prcessed_cluster_ids;
    my $total_clusters=0;
    foreach my $cluster_id(split ",", $clusters_breakpoint_1){
        next if exists $prcessed_cluster_ids{$cluster_id};
        $prcessed_cluster_ids{$cluster_id}=0;
        $total_clusters++;
    }
    foreach my $cluster_id(split ",", $clusters_breakpoint_2){
        next if exists $prcessed_cluster_ids{$cluster_id};
        $prcessed_cluster_ids{$cluster_id}=0;
        $total_clusters++;
    }
    
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
        $total_reads_in_each_cluster_breakpoint_1, $total_reads_in_each_cluster_breakpoint_2);
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

    -f --flanking_length    INT     Breakpoint flanking size to calculate standard deviation [default: $flanking_length]
    -s --sd_cutoff          FLOAT   Standard deviation cutoff to filter fusions  [default: $sd_cutoff]
    -h --help     	                 Print this help information
HELP
    exit(-1);
}