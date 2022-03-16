#!/usr/bin/env perl

# 2022-03-15

# 1. Function
# Cluster reads within defined --window_size

# 2. Input
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: junction type:  -1=read pair (between the mates), 1=split reads
# column 8: read name

# 3. Output
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: junction type:  -1=read pair (between the mates), 1=split reads
# column 8: read name
# column 9: cluster id



use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $help;
my $window_size=1000000;
GetOptions(
    'w|window_size=i'=>	\$window_size,
    'h|help'=>        	\$help
)||usage(); 
usage () if defined $help;

$ARGV[0]='-' unless defined $ARGV[0];
# sort fusion breakpoint input by chr_left, strand_left, chr_right, strand_right, pos_left, pos_right
# input have header
system("sed -n '2,\$p' $ARGV[0] | sort -k1,1 -k3,3 -k4,4 -k6,6 -k2,2n -k5,5n >temp_discordant_breakpoint_sorted.tsv");

open IN, "temp_discordant_breakpoint_sorted.tsv" or die "Can't open temp_discordant_breakpoint_sorted.tsv:$!";
say join "\t", ("Chr_left", "Pos_left", "Strand_left", "Chr_right", "Pos_right", "Strand_right", "Junction_type", "Read_name", "Cluster_id");

my @Cluster_array; # reads belong to the same cluster
my @Unprocessed_array; # reads taken from file but not belong the current cluster
my $cluster_id=1;

chomp(my $input_line=<IN>);
push @Cluster_array, $input_line;
my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split "\t", $input_line;
my ($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
    $chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);

while(<IN>){
    chomp;
    $input_line=$_;
    ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split "\t", $input_line;
    if($chr_left ne $compare_chr_left || $strand_left ne $compare_strand_left || $chr_right ne $compare_chr_right || $strand_right ne $compare_strand_right){ # belong to different cluster from this reads
	&process_cluster;
	push @Cluster_array, $input_line;
	($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
	    $chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    }else{
	my $distance_left=abs($pos_left-$compare_pos_left);
	my $distance_right=abs($pos_right-$compare_pos_right);
	if($distance_left>$window_size && $distance_right>$window_size){ # belong to different cluster from this reads
	    &process_cluster;
	    push @Cluster_array, $input_line;
	    ($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
		$chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	}else{
	    if($distance_left<=$window_size && $distance_right<=$window_size){
		push @Cluster_array, $input_line;
	    }else{
		push @Unprocessed_array, $input_line;
	    }
	}
    }
}
close(IN);
system("rm temp_discordant_breakpoint_sorted.tsv");

&process_cluster;

sub process_cluster{
    foreach my $read_info (@Cluster_array){
	say join "\t", ($read_info, $cluster_id);
    }
    
    $cluster_id++;
    @Cluster_array=();
    until(@Unprocessed_array==0){
	my $process_read_info=shift @Unprocessed_array;
	push @Cluster_array, $process_read_info;
	my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split "\t", $process_read_info;
	my ($compare_chr_left, $compare_pos_left, $compare_strand_left, $compare_chr_right, $compare_pos_right, $compare_strand_right)=(
	    $chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
	
	my @Temp_unprocessed_array;
	foreach $process_read_info (@Unprocessed_array){
	    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name)=split "\t", $process_read_info;
	    my $distance_left=abs($pos_left-$compare_pos_left);
	    my $distance_right=abs($pos_right-$compare_pos_right);
	    
	    if($distance_left<=$window_size && $distance_right<=$window_size){
		push @Cluster_array, $process_read_info;
	    }else{
		push @Temp_unprocessed_array, $process_read_info;
	    }
	}
	
	foreach my $read_info (@Cluster_array){
	    say join "\t", ($read_info, $cluster_id);
	}
	
	$cluster_id++;
	@Cluster_array=();
	@Unprocessed_array=@Temp_unprocessed_array;
    }
}

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to cluster read within defined window size
Usage: perl $scriptName input >output

    -w	--window_size	discordant reads within this window size will be clustered [default: $window_size]
    -h	--help		print this help information screen

HELP
    exit(-1);
}