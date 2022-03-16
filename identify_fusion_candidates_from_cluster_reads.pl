#!/usr/bin/env perl

# 2022-03-15

# 1. Function
# Identify fusion candidates from clustered reads 

# 2. Input
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: junction type: -1=read pair (between the mates), 1=split reads
# column 8: read name
# column 9: cluster id

# 3. Output
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: junction type: -1=read pair (between the mates), 1=split reads
# column 8: cluster id
# column 9: split read count
# column 10: read pair count
# column 11: split reads
# column 12: read pairs



use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $help;
GetOptions(
    'h|help'=>        	\$help
)||usage(); 
usage () if defined $help;


say join "\t", ("Chr_left", "Pos_left", "Strand_left", "Chr_right", "Pos_right", "Strand_right", "Junction_type", "Cluster_id", "Split_read_count", "Read_pair_count", "Split_reads", "Read_pairs");
$ARGV[0]='-' unless defined $ARGV[0];
open IN,  $ARGV[0] or die "Can't open  $ARGV[0]:$!";
# input have header
<IN>;

my @Split_reads_info; # split reads
my @Span_reads_info; # read pairs

chomp(my $input_line=<IN>);
my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name, $cluster_id)=split "\t", $input_line;
if($junc_type == "-1"){
    push @Span_reads_info, $input_line;
}else{
    push @Split_reads_info, $input_line;
}
my $current_cluster_id=$cluster_id;

while(<IN>){
    chomp;
    $input_line=$_;
    ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $junc_type, $read_name, $cluster_id)=split "\t", $input_line;
    if($cluster_id != $current_cluster_id){ # different cluster
        &process_cluster;
        $current_cluster_id=$cluster_id;
    }
    
    if($junc_type == "-1"){
        push @Span_reads_info, $input_line;
    }else{
        push @Split_reads_info, $input_line;
    }
}
close(IN);

&process_cluster;

sub process_cluster{
    my %number_of_span_read_belong_to_split_read;
    until(@Split_reads_info==0){		
        my $process_line=shift @Split_reads_info;
        my ($output_chr_left, $output_pos_left, $output_strand_left, $output_chr_right, $output_pos_right, $output_strand_right, $output_junc_type, $output_read_name, $output_cluster_id)=split "\t", $process_line;
        my $output_split_reads=$output_read_name;
        my $output_split_read_count=1;
        
	my @Temp_unprocessed_reads;
        foreach $process_line (@Split_reads_info){
            my ($process_chr_left, $process_pos_left, $process_strand_left, $process_chr_right, $process_pos_right, $process_strand_right, $process_junc_type, $process_read_name, $process_cluster_id)=split "\t", $process_line;
            if($output_pos_left==$process_pos_left && $output_pos_right==$process_pos_right){
                $output_split_reads.=",".$process_read_name;
                $output_split_read_count++;
            }else{
		push @Temp_unprocessed_reads, $process_line;
	    }
        }
	@Split_reads_info=@Temp_unprocessed_reads;
        
        my $output_read_pairs="NA";
        my $output_read_pair_count=0;
        foreach $process_line (@Span_reads_info){
            my ($process_chr_left, $process_pos_left, $process_strand_left, $process_chr_right, $process_pos_right, $process_strand_right, $process_junc_type, $process_read_name, $process_cluster_id)=split "\t", $process_line;
            $number_of_span_read_belong_to_split_read{$process_read_name}=0 if !exists $number_of_span_read_belong_to_split_read{$process_read_name};
	    my ($is_left_support, $is_right_support)=("0") x 2;
            $is_left_support=1 if ($output_strand_left eq "+" && $process_pos_left<=$output_pos_left) || ($output_strand_left eq "-" && $process_pos_left>=$output_pos_left);
            $is_right_support=1 if ($output_strand_right eq "+" && $process_pos_right<=$output_pos_right) || ($output_strand_right eq "-" && $process_pos_right>=$output_pos_right);
            if($is_left_support==1 && $is_right_support==1){
                $output_read_pairs.=",".$process_read_name;
                $output_read_pair_count++;
                $number_of_span_read_belong_to_split_read{$process_read_name}+=1;
            }
        }
        
        $output_read_pairs=~s/^NA,//;
        say join "\t", ($output_chr_left, $output_pos_left, $output_strand_left, $output_chr_right, $output_pos_right, $output_strand_right, $output_junc_type, $output_cluster_id, $output_split_read_count, $output_read_pair_count, $output_split_reads, $output_read_pairs);
    }
    
    # for read pairs not belong to split reads, select most likely position
    @Span_reads_info=();
    foreach my $read (sort keys %number_of_span_read_belong_to_split_read){
        push @Span_reads_info, $read if $number_of_span_read_belong_to_split_read{$read}==0;
    }
    if(@Span_reads_info>0){
        my $output_read_pair_count=@Span_reads_info;
        
        my $process_line=shift @Span_reads_info;
        my ($output_chr_left, $output_pos_left, $output_strand_left, $output_chr_right, $output_pos_right, $output_strand_right, $output_junc_type, $output_read_name, $output_cluster_id)=split "\t", $process_line;
        my $output_read_pairs=$output_read_name;
        
        foreach $process_line (@Span_reads_info){
            my ($process_chr_left, $process_pos_left, $process_strand_left, $process_chr_right, $process_pos_right, $process_strand_right, $process_junc_type, $process_read_name, $process_cluster_id)=split "\t", $process_line;
            $output_pos_left=$process_pos_left if ($process_strand_left eq '+' && $process_pos_left>$output_pos_left) || ($process_strand_left eq '-' && $process_pos_left<$output_pos_left);
            $output_pos_right=$process_pos_right if ($process_strand_right eq '+' && $process_pos_right>$output_pos_right) || ($process_strand_right eq '-' && $process_pos_right<$output_pos_right);
            $output_read_pairs.=",".$process_read_name;
        }
        
        say join "\t", ($output_chr_left, $output_pos_left, $output_strand_left, $output_chr_right, $output_pos_right, $output_strand_right, $output_junc_type, $output_cluster_id, "0", $output_read_pair_count, "NA", $output_read_pairs);
        @Span_reads_info=();
    }
} 

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to identify fusion candidates from cluster reads
Usage: perl $scriptName input >fusion_candidates

    -h	--help	print this help information screen

HELP
    exit(-1);
}