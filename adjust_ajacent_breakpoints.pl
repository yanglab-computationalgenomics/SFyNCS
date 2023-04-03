#!/usr/bin/env perl

# 2022-11-03

# 1. Function
# Merge split reads' adjacent breakpoints (junction type is 1)
# 1.1. junction ChrA:PosA1:-__ChrB:PosB1:+ should merge with ChrB:PosB2:-__chrA:PosA2:+ if they are within --adjacent_distance
# 1.2. change junction type 0/2 to 1
# 1.3. when merging junction, selecte jucntion type first (1>2>0), then base on read count, select most left junction if junction type and read count are same

# 2. Input
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 4: chromosome of breakpoint 2
# column 5: breakpoint 2 position (1-based)
# column 6: breakpoint 2 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 7: junction type: -1=encompassing junction (between the mates), 1=split reads
# column 8: read name

# 3. Output
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 4: chromosome of breakpoint 2
# column 5: breakpoint 2 position (1-based)
# column 6: breakpoint 2 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 7: junction type: -1=encompassing junction (between the mates), 1=split reads
# column 8: read name


use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $adjacent_distance=5;
GetOptions(
    'a|adjacent_distance=i'=>	\$adjacent_distance,
    'h|help'    => sub{usage()}
)||usage();


# sort split read by chr_breakpoint_1, strand_breakpoint_1, chr_breakpoint_2, strand_breakpoint_2, pos_breakpoint_1, pos_breakpoint_2
# input have header
system("awk 'BEGIN{FS=OFS=\"\\t\"} NR>1 && \$7!=-1' $ARGV[0] | sort -k1,1 -k3,3 -k4,4 -k6,6 -k2,2n -k5,5n >temp_discordant_breakpoint_sorted.tsv");
open IN, "temp_discordant_breakpoint_sorted.tsv" or die "Can't open temp_discordant_breakpoint_sorted.tsv:$!";

my @Cluster_array; # reads belong to the same cluster
my @Unprocessed_array; # reads taken from file but not belong to the current cluster
my %fusion_type; # store fusion type
my %fusion_count; # store the number of read that support fusion
my %merge_fusion; # store select fusion

chomp(my $input_line=<IN>);

my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split "\t", $input_line;
$fusion_type{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}=$junc_type;
$fusion_count{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}=1;
push @Cluster_array, join ",", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);

# only one position needed 
my ($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
    $chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);


while(<IN>){
    chomp;
    ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split "\t", $_;
    $fusion_type{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}=$junc_type;
    if(exists $fusion_count{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}){
								$fusion_count{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}++;
    }else{
								$fusion_count{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}=1;
    }
    
    if($chr_breakpoint_1 ne $compare_chr_breakpoint_1 || $strand_breakpoint_1 ne $compare_strand_breakpoint_1 || $chr_breakpoint_2 ne $compare_chr_breakpoint_2 ||
							$strand_breakpoint_2 ne $compare_strand_breakpoint_2){ # belong to different cluster from this reads
								&process_cluster;
								push @Cluster_array, join ",", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
								($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
								$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
    }else{
								my $distance_breakpoint_1=abs($pos_breakpoint_1-$compare_pos_breakpoint_1);
								my $distance_breakpoint_2=abs($pos_breakpoint_2-$compare_pos_breakpoint_2);
								if($distance_breakpoint_1>$adjacent_distance && $distance_breakpoint_2>$adjacent_distance){ # belong to different cluster from this reads
												&process_cluster;
												push @Cluster_array, join ",", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
												($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
												$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
								}else{
												if($distance_breakpoint_1<=$adjacent_distance && $distance_breakpoint_2<=$adjacent_distance){
																push @Cluster_array, join ",", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
																# comment code: try to change min and max cluster position, didn't do it in current version
																#$compare_pos_breakpoint_1=$pos_breakpoint_1 if $pos_breakpoint_1 > $compare_pos_breakpoint_1;
																#$compare_start_breakpoint_2=$pos_breakpoint_2 if $pos_breakpoint_2 < $compare_start_breakpoint_2;
																#$compare_end_breakpoint_2=$pos_breakpoint_2 if $pos_breakpoint_2 > $compare_end_breakpoint_2;
												}else{
																push  @Unprocessed_array, join ",", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
												}
								}
    }
}

close(IN);
system("rm temp_discordant_breakpoint_sorted.tsv");

&process_cluster;

say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2", "Junction_type", "Read_name");
open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split "\t", $_;
				#$junc_type=1 if $junc_type!=-1;
    if($junc_type==-1 || !exists $merge_fusion{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}){
								say join "\t", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name);
    }else{
								my ($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2)=split ",", $merge_fusion{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2};
								say join "\t", ($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2, $junc_type, $read_name);
    }
}



sub process_cluster{
    &select_fusion_from_adjacent(\@Cluster_array) if @Cluster_array>1; # do nothing if no nearby fusion
    
    @Cluster_array=();
    until(@Unprocessed_array==0){
								my $process_read_info=shift @Unprocessed_array;
								push @Cluster_array, $process_read_info;
								my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split ",", $process_read_info;
								my ($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
												$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
	
								my @Temp_unprocessed_array;
								foreach $process_read_info (@Unprocessed_array){
												my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split ",", $process_read_info;
												my $distance_breakpoint_1=abs($pos_breakpoint_1-$compare_pos_breakpoint_1);
												my $distance_breakpoint_2=abs($pos_breakpoint_2-$compare_pos_breakpoint_2);
	    
												if($distance_breakpoint_1<=$adjacent_distance && $distance_breakpoint_2<=$adjacent_distance){
																push @Cluster_array, $process_read_info;
												}else{
																push @Temp_unprocessed_array, $process_read_info;
												}
								}
	
								&select_fusion_from_adjacent(\@Cluster_array) if @Cluster_array>1; # do nothing if no nearby fusion
	
								@Cluster_array=();
								@Unprocessed_array=@Temp_unprocessed_array;
    }
}

sub select_fusion_from_adjacent{
    my ($temp)=@_;
    my @input_array=@$temp;
    my ($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2);
    my ($selected_junction_type, $selected_read_count)=("NA", 0);
    my %temp_hash; # store different fusion
    # selection base on jucntion type first (1>2>0), then base on read count, select most left junction if all same
    foreach my $read_info (@Cluster_array){
								my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2)=split ",", $read_info;
								$temp_hash{$read_info}=0;
								my $junc_type=$fusion_type{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2};
								my $read_count=$fusion_count{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2};
								if($selected_junction_type eq "NA"){
												($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2)=(
												$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
												$selected_junction_type=$junc_type;
												$selected_read_count=$read_count;
								}else{
												if($junc_type eq $selected_junction_type){
																if($read_count > $selected_read_count){
																				($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2)=(
																				$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
																				$selected_read_count=$read_count;
																}
												}else{
																if(($junc_type==1) || ($junc_type==2 && $selected_junction_type==0)){
																				($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2)=(
																				$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
																				$selected_read_count=$read_count;
																}
												}
								}
    }
    
    if(keys %temp_hash>1){
								foreach my $read_info (@Cluster_array){
												my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2)=split ",", $read_info;
												$merge_fusion{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2}=join ",",($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2);
												# debug
												#say STDERR join "\t", ($selected_chr_breakpoint_1, $selected_pos_breakpoint_1, $selected_strand_breakpoint_1, $selected_chr_breakpoint_2, $selected_pos_breakpoint_2, $selected_strand_breakpoint_2, $chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $fusion_type{$chr_breakpoint_1}{$pos_breakpoint_1}{$strand_breakpoint_1}{$chr_breakpoint_2}{$strand_breakpoint_2}{$pos_breakpoint_2});
								}
    }
}



sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to merge adjacent fusion breakpoints
Usage: perl $scriptName input >output

    -a	--adjacent_distance	INT		Breakpoints within this distance will be adjusted [default: $adjacent_distance]
    -h	--help																			Print this help information screen

HELP
    exit(-1);
}

