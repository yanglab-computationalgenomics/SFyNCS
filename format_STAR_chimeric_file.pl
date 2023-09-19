#!/usr/bin/env perl

# 2022-11-03

# 1. Function
# Add "chr" if chromosome is numeral
# Change upstream/downstream position in STAR chimeirc file to breakpoint position
# Format breakpoint strand: + means left or upstream of the breakpoint will be used in fusion, while - means right or downstream of breakpoint
# For a discordant read, sort breakpoints by chromosome if breakpoints locate in different chromosomes
# For a discordant read, sort breakpoints by breakpoint position if breakpoints locate in same chromosome

# to do in next version:
# maybe to do in next version: delete $junc_type=1 if $junc_type!=-1;
# maybe change chrX/chrx to chr23 and chrY/chry to chr24

# 2. Input from chimeric file (taking from star document):
# Every line contains one chimerically aligned read, e.g.:
# chr22 23632601 + chr9 133729450 + 1 0 0 SINATRA-0006:3:3:6387:56650 23632554 47M29S 133729451 47S29M40p76M
# The first 9 columns give information about the chimeric junction. Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like andlignments are given with respect to the (+) strand:
# column 1: chromosome of the donor
# column 2: rst base of the intron of the donor (1-based)
# column 3: strand of the donor (dornor read aligned strand)
# column 4: chromosome of the acceptor
# column 5: rst base of the intron of the acceptor (1-based)
# column 6: strand of the acceptor (acceptor read aligned strand)
# column 7: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC, 0=any other motif
# column 8: repeat length to the left of the junction
# column 9: repeat length to the right of the junction
# column 10: read name
# column 11: rst base of the rst segment (on the + strand)
# column 12: CIGAR of the rst segment
# column 13: rst base of the second segment
# column 14: CIGAR of the second segment

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

say join "\t", ("Chr_breakpoint_1", "Pos_breakpoint_1", "Strand_breakpoint_1", "CIGAR_breakpoint_1", "Chr_breakpoint_2", "Pos_breakpoint_2", "Strand_breakpoint_2", "CIGAR_breakpoint_2", "Junction_type", "Read_name");
while(<IN>){
    chomp;
    next if $_=~/^chr_donorA/; # in case have header 
    next if $_=~/^#/;
    my ($chr_d, $rst_d, $strand_d, $chr_a, $rst_a, $strand_a, $junc_type, $rep_len_d, $rep_len_a, $read_name, $align_d, $cigar_d, $align_a, $cigar_a)=split "\t",$_;
    $chr_a="chr".$chr_a if $chr_a!~/^chr/;
    $chr_a=~s/Chr/chr/;
    #$chr_a=~s/chrX/chr23/;
    #$chr_a=~s/chrx/chr23/;
    #$chr_a=~s/chrY/chr24/;
    #$chr_a=~s/chry/chr24/;
    $chr_d="chr".$chr_d if $chr_d!~/^chr/;
    $chr_d=~s/Chr/chr/;
    #$chr_d=~s/chrX/chr23/;
    #$chr_d=~s/chrx/chr23/;
    #$chr_d=~s/chrY/chr24/;
    #$chr_d=~s/chry/chr24/;
    
    # maybe delete in next version
    $junc_type=1 if $junc_type!=-1;
    
    # change upstream/downstream site to breakpoint
    $rst_d = ($strand_d eq '+') ? --$rst_d : ++$rst_d; 
    $rst_a = ($strand_a  eq '+') ? ++$rst_a : --$rst_a;
    
    # change breakpoint strand: + means left or upstream of the breakpoint will be used in fusion, while - means right or downstream of breakpoint
    # no need to change strand_d
    $strand_a = ($strand_a eq '+') ? "-" : "+";
    
    # order breakpoint
    my ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cigar_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cigar_breakpoint_2)=($chr_d, $rst_d, $strand_d, $cigar_d, $chr_a, $rst_a, $strand_a, $cigar_a);
    if($chr_d eq $chr_a){
								($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cigar_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cigar_breakpoint_2)=($chr_a, $rst_a, $strand_a, $cigar_a, $chr_d, $rst_d, $strand_d, $cigar_d) if $rst_a < $rst_d;
    }else{
								($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cigar_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cigar_breakpoint_2)=($chr_a, $rst_a, $strand_a, $cigar_a, $chr_d, $rst_d, $strand_d, $cigar_d) if substr($chr_a, 3) < substr($chr_d, 3);
    }
    say join "\t", ($chr_breakpoint_1, $pos_breakpoint_1, $strand_breakpoint_1, $cigar_breakpoint_1, $chr_breakpoint_2, $pos_breakpoint_2, $strand_breakpoint_2, $cigar_breakpoint_2, $junc_type, $read_name);
}
close(IN);

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to format STAR chimeric input
Usage: perl $scriptName Chimeric.out.junction >output

    -h	--help	Print this help information screen

HELP
    exit(-1);
}
