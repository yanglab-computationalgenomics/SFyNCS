#!/usr/bin/env perl

# 2022-03-15

# 1. Function
# Remove multiple-mapped reads

# 2. Input
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: CIGAR of left breakpoint
# column 5: chromosome of right segment
# column 6: right segment breakpoint
# column 7: strand of the right segment
# column 8: CIGAR of right segment
# column 9: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC, 0=any other motif
# column 10: read name

# 3. Output
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: CIGAR of left segment
# column 5: chromosome of right segment
# column 6: right segment breakpoint
# column 7: strand of the right segment
# column 8: CIGAR of right segment
# column 9: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC, 0=any other motif
# column 10: read name


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

my %hash;
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right, $junc_type, $read_name)=split "\t",$_;
    my $temp_key=join "_", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    $hash{$read_name}{$temp_key}=0;
}

seek(IN, 0, 0);
# input have header
<IN>;
say join "\t", ("Chr_left", "Pos_left", "Strand_left", "CIGAR_left", "Chr_right", "Pos_right", "Strand_right", "CIGAR_righ", "Junction_type", "Read_name");
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right, $junc_type, $read_name)=split "\t",$_;
    my $temp_key=join "_", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right);
    my $read_align_count=keys %{$hash{$read_name}};
    say $_ if $read_align_count==1;
}

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to remove multiple-mapped discordant reads
Usage: perl $scriptName input >output

    -h	--help	print this help information screen

HELP
    exit(-1);
}