#!/usr/bin/env perl

# 2022-11-03

# 1. Function
# Remove multiple-mapped reads

# 2. Input
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 4: CIGAR of breakpoint 1
# column 5: chromosome of breakpoint 2
# column 6: breakpoint 2 position (1-based)
# column 7: breakpoint 2 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 8: CIGAR of breakpoint 2
# column 9: junction type: -1=encompassing junction (between the mates), 1=split reads
# column 10: read name

# 3. Output
# column 1: chromosome of breakpoint 1
# column 2: breakpoint 1 position (1-based)
# column 3: breakpoint 1 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 4: CIGAR of breakpoint 1
# column 5: chromosome of breakpoint 2
# column 6: breakpoint 2 position (1-based)
# column 7: breakpoint 2 strand (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 8: CIGAR of breakpoint 2
# column 9: junction type: -1=encompassing junction (between the mates), 1=split reads
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
    my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cigar_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cigar_breakpoint_2, $junc_type, $read_name)=split "\t",$_;
    my $temp_key=join "_", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
    $hash{$read_name}{$temp_key}=0;
}

seek(IN, 0, 0);
# input have header
<IN>;
say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "CIGAR_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2", "CIGAR_breakpoint_2", "Junction_type", "Read_name");
while(<IN>){
    chomp;
    my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cigar_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cigar_breakpoint_2, $junc_type, $read_name)=split "\t",$_;
    my $temp_key=join "_", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2);
    my $read_align_count=keys %{$hash{$read_name}};
    say $_ if $read_align_count==1;
}

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to remove multiple-aligned discordant reads
Usage: perl $scriptName input >output

    -h	--help	Print this help information screen

HELP
    exit(-1);
}