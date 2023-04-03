#!/usr/bin/env perl

# 2022-11-03

# 1. Function
# Cluster reads within defined --window_size

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
# sort fusion breakpoint input by chr_breakpoint_1, strand_breakpoint_1, chr_breakpoint_2, strand_breakpoint_2, pos_breakpoint_1, pos_breakpoint_2
# input have header
system("sed -n '2,\$p' $ARGV[0] | sort -k1,1 -k3,3 -k4,4 -k6,6 -k2,2n -k5,5n >temp_discordant_breakpoint_sorted.tsv");

open IN, "temp_discordant_breakpoint_sorted.tsv" or die "Can't open temp_discordant_breakpoint_sorted.tsv:$!";
say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2", "Junction_type", "Read_name", "Cluster_id");

my @Cluster_array; # reads belong to the same cluster
my @Unprocessed_array; # reads taken from file but not belong the current cluster
my $cluster_id=1;

chomp(my $input_line=<IN>);
push @Cluster_array, $input_line;
my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split "\t", $input_line;
my ($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
    $chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);

while(<IN>){
    chomp;
    $input_line=$_;
    ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split "\t", $input_line;
    if($chr_breakpoint_1 ne $compare_chr_breakpoint_1 || $strand_breakpoint_1 ne $compare_strand_breakpoint_1 ||
							$chr_breakpoint_2 ne $compare_chr_breakpoint_2 || $strand_breakpoint_2 ne $compare_strand_breakpoint_2){ # belong to different cluster from this reads
								&process_cluster;
								push @Cluster_array, $input_line;
								($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
												$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
    }else{
								my $distance_breakpoint_1=abs($pos_breakpoint_1-$compare_pos_breakpoint_1);
								my $distance_breakpoint_2=abs($pos_breakpoint_2-$compare_pos_breakpoint_2);
								if($distance_breakpoint_1>$window_size && $distance_breakpoint_2>$window_size){ # belong to different cluster from this reads
												&process_cluster;
												push @Cluster_array, $input_line;
												($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
																$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
								}else{
												if($distance_breakpoint_1<=$window_size && $distance_breakpoint_2<=$window_size){
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
								my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split "\t", $process_read_info;
								my ($compare_chr_breakpoint_1, $compare_pos_breakpoint_1, $compare_strand_breakpoint_1, $compare_chr_breakpoint_2, $compare_pos_breakpoint_2, $compare_strand_breakpoint_2)=(
								$chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
	
								my @Temp_unprocessed_array;
								foreach $process_read_info (@Unprocessed_array){
												my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $junc_type, $read_name)=split "\t", $process_read_info;
												my $distance_breakpoint_1=abs($pos_breakpoint_1-$compare_pos_breakpoint_1);
												my $distance_breakpoint_2=abs($pos_breakpoint_2-$compare_pos_breakpoint_2);
	    
												if($distance_breakpoint_1<=$window_size && $distance_breakpoint_2<=$window_size){
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

    -w	--window_size	INT	Split reads and read pairs within this distance will be clustered together [default: $window_size]
    -h	--help			Print this help information screen

HELP
    exit(-1);
}