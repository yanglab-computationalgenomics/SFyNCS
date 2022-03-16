#!/usr/bin/env perl

# 2022-03-15

# 1. Function
# Merge split reads' adjacent breakpoints (junction type is 0, 1, or 2)
# 1.1. junction ChrA:PosA1:-__ChrB:PosB1:+ should merge with ChrB:PosB2:-__chrA:PosA2:+ if they are within --adjacent_distance
# 1.2. change junction type 0/2 to 1
# 1.3. when merging junction, selecte jucntion type first (1>2>0), then base on read count, select most left junction if junction type and read count are same

# 2. Input
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC, 0=any other motif
# column 8: read name

# 3. Output
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: junction type: -1=encompassing junction (between the mates), 1=split read
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


# sort split read by chr_left, strand_left, chr_right, strand_right, pos_left, pos_right
# input have header
system("awk 'BEGIN{FS=OFS=\"\\t\"} NR>1 && \$7!=-1' $ARGV[0] | sort -k1,1 -k3,3 -k4,4 -k6,6 -k2,2n -k5,5n >temp_discordant_breakpoint_sorted.tsv");
open IN, "temp_discordant_breakpoint_sorted.tsv" or die "Can't open temp_discordant_breakpoint_sorted.tsv:$!";

my @Cluster_array; # reads belong to the same cluster
my @Unprocessed_array; # reads taken from file but not belong to the current cluster
my %fusion_type; # store fusion type
my %fusion_count; # store the number of read that support fusion
my %merge_fusion; # store select fusion

chomp(my $input_line=<IN>);

my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split "\t", $input_line;
$fusion_type{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}=$junc_type;
$fusion_count{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}=1;
push @Cluster_array, join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);

# only one position needed 
my ($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
    $chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);


while(<IN>){
    chomp;
    ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split "\t", $_;
    $fusion_type{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}=$junc_type;
    if(exists $fusion_count{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}){
	$fusion_count{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}++;
    }else{
	$fusion_count{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}=1;
    }
    
    if($chr_left ne $compare_chr_left || $strand_left ne $compare_strand_left || $chr_right ne $compare_chr_right || $strand_right ne $compare_strand_right){ # belong to different cluster from this reads
	&process_cluster;
	push @Cluster_array, join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
	    $chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    }else{
	my $distance_left=abs($pos_left-$compare_pos_left);
	my $distance_right=abs($pos_right-$compare_pos_right);
	if($distance_left>$adjacent_distance && $distance_right>$adjacent_distance){ # belong to different cluster from this reads
	    &process_cluster;
	    push @Cluster_array, join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	    ($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
		$chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	}else{
	    if($distance_left<=$adjacent_distance && $distance_right<=$adjacent_distance){
		push @Cluster_array, join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
		#$compare_pos_left=$pos_left if $pos_left > $compare_pos_left;
		#$compare_start_right=$pos_right if $pos_right < $compare_start_right;
		#$compare_end_right=$pos_right if $pos_right > $compare_end_right;
	    }else{
		push  @Unprocessed_array, join ",", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	    }
	}
    }
}

close(IN);
system("rm temp_discordant_breakpoint_sorted.tsv");

&process_cluster;

say join "\t", ("Chr_left", "Pos_left", "Strand_left", "Chr_right", "Pos_right", "Strand_right", "Junction_type", "Read_name");
open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split "\t", $_;
    $junc_type=1 if $junc_type!=-1;
    if($junc_type==-1 || !exists $merge_fusion{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}){
	say join "\t", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name);
    }else{
	my ($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right)=split ",", $merge_fusion{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right};
	say join "\t", ($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right, $junc_type, $read_name);
    }
}



sub process_cluster{
    &select_fusion_from_adjacent(\@Cluster_array) if @Cluster_array>1; # do nothing if no nearby fusion
    
    @Cluster_array=();
    until(@Unprocessed_array==0){
	my $process_read_info=shift @Unprocessed_array;
	push @Cluster_array, $process_read_info;
	my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split ",", $process_read_info;
	my ($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
	    $chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	
	my @Temp_unprocessed_array;
	foreach $process_read_info (@Unprocessed_array){
	    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split ",", $process_read_info;
	    my $distance_left=abs($pos_left-$compare_pos_left);
	    my $distance_right=abs($pos_right-$compare_pos_right);
	    
	    if($distance_left<=$adjacent_distance && $distance_right<=$adjacent_distance){
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
    my ($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right);
    my ($selected_junction_type, $selected_read_count)=("NA", 0);
    my %temp_hash; # store different fusion
    # selection base on jucntion type first (1>2>0), then base on read count, select most left junction if all same
    foreach my $read_info (@Cluster_array){
	my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right)=split ",", $read_info;
	$temp_hash{$read_info}=0;
	my $junc_type=$fusion_type{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right};
	my $read_count=$fusion_count{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right};
	if($selected_junction_type eq "NA"){
	    ($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right)=(
		$chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	    $selected_junction_type=$junc_type;
	    $selected_read_count=$read_count;
	}else{
	    if($junc_type eq $selected_junction_type){
		if($read_count > $selected_read_count){
		    ($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right)=(
			$chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
		    $selected_read_count=$read_count;
		}
	    }else{
		if(($junc_type==1) || ($junc_type==2 && $selected_junction_type==0)){
		    ($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right)=(
			$chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
		    $selected_read_count=$read_count;
		}
	    }
	}
    }
    
    if(keys %temp_hash>1){
	foreach my $read_info (@Cluster_array){
	    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right)=split ",", $read_info;
	    $merge_fusion{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right}=join ",",($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right);
	    # debug
	    #say STDERR join "\t", ($selected_chr_left, $selected_pos_left, $selected_strand_left, $selected_chr_right, $selected_pos_right, $selected_strand_right, $chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $fusion_type{$chr_left}{$pos_left}{$strand_left}{$chr_right}{$strand_right}{$pos_right});
	}
    }
}



sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to merge adjacent fusion breakpoints
Usage: perl $scriptName input >output

    -a	--adjacent_distance	breakpoints within this range will be merged [default: $adjacent_distance]
    -h	--help			print this help information screen

HELP
    exit(-1);
}

