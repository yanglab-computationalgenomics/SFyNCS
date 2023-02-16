#!/usr/bin/env perl

# 2023-02-14

# 1. function
# discard duplicate reads
# get the number of reads that supportes junctions

# 2. input
#column 1: chromosome of breakpoint 1
#column 2: position of breakpoint 1 (1-base)
#column 3: strand of the breakpoint 1 (different from input strand, + means left of the site will be used, while - means right of site will be used)
#column 4: CIGAR of breakpoint 1
#column 5: chromosome of breakpoint 2
#column 6: position of breakpoint 2 (1-base)
#column 7: strand of the breakpoint 2 (different from input strand, + means left of the site will be used, while - means right of site will be used)
#column 8: CIGAR of breakpoint 2
#column 9: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC, 0=any other motif
#column 10: read name

# 3. output
#column 1: chromosome of breakpoint 1
#column 2: position of breakpoint 1 (1-base)
#column 3: strand of the breakpoint 1 (different from input strand, + means left of the site will be used, while - means right of site will be used)
#column 4: chromosome of breakpoint 2
#column 5: position of breakpoint 2 (1-base)
#column 6: strand of the breakpoint 2 (different from input strand, + means left of the site will be used, while - means right of site will be used)
#column 7: read count


use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $help;
GetOptions(
    'h|help'=>        \$help
)||usage(); 
usage () if defined $help;

$ARGV[0]='-' unless defined $ARGV[0];
open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";

my %duplication;
my %read_count;
# input didn't have header
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right, $junc_type, $read_name)=split "\t", $_;
    my $temp_key=join "_", ($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right);
    next if exists $duplication{$temp_key};
    $duplication{$temp_key}=0;
    $read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}{$pos_left}{$pos_right}{$read_name}=0;
}
close(IN);

say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2", "Read_count");
foreach my $strand_left (sort keys %read_count){
    foreach my $strand_right (sort keys %{$read_count{$strand_left}}){
        foreach my $chr_left (sort keys %{$read_count{$strand_left}{$strand_right}}){
            foreach my $chr_right (sort keys %{$read_count{$strand_left}{$strand_right}{$chr_left}}){
                foreach my $pos_left (sort keys %{$read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}}){
                    foreach my $pos_right (sort keys %{$read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}{$pos_left}}){
                       my $temp_read_count=keys %{$read_count{$strand_left}{$strand_right}{$chr_left}{$chr_right}{$pos_left}{$pos_right}};
                       say join "\t", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $temp_read_count);
                    }
                }
            }
        }
    }
}

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to calculate the number of reads that supports junction
Usage: perl $scriptName format_chimeric.tsv >output

    -h	--help				print this help information screen

HELP
    exit(-1);
}
