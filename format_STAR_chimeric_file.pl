#!/usr/bin/env perl

# 2022-03-15

# 1. Function
# Add "chr" if chromosome is numeral
# Change chrX/chrx to chr23 and chrY/chry to chr24
# Change upstream/downstream position in STAR chimeirc file to breakpoint position
# Format breakpoint strand: + means left or upstream of the breakpoint will be used in fusion, while - means right or downstream of breakpoint
# For a discordant read, sort breakpoints by chromosome if breakpoints locate in different chromosomes
# For a discordant read, sort breakpoints by breakpoint position if breakpoints locate in same chromosome

# to do: delete $junc_type=1 if $junc_type!=-1;

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
# column 1: chromosome of the left segment
# column 2: left segment breakpoint (1-based)
# column 3: strand of the left segment (different from input strand, + means left of the site will be used, while - means right of site will be used)
# column 4: CIGAR of left segment
# column 5: chromosome of right segment
# column 6: right segment breakpoint (1-based)
# column 7: strand of the right segment (different from input strand, + means left of the site will be used, while - means right of site will be used)
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

say join "\t", ("Chr_left", "Pos_left", "Strand_left", "CIGAR_left", "Chr_right", "Pos_right", "Strand_right", "CIGAR_righ", "Junction_type", "Read_name");
# input didn't have header
while(<IN>){
    chomp;
    next if $_=~/^#/;
    my ($chr_d, $rst_d, $strand_d, $chr_a, $rst_a, $strand_a, $junc_type, $rep_len_d, $rep_len_a, $read_name, $align_d, $cigar_d, $align_a, $cigar_a)=split "\t",$_;
    #$chr_a="chr".$chr_a if $chr_a!~/^chr/;
    #$chr_a=~s/Chr/chr/;
    #$chr_a=~s/chrX/chr23/;
    #$chr_a=~s/chrx/chr23/;
    #$chr_a=~s/chrY/chr24/;
    #$chr_a=~s/chry/chr24/;
    #$chr_d="chr".$chr_d if $chr_d!~/^chr/;
    #$chr_d=~s/Chr/chr/;
    #$chr_d=~s/chrX/chr23/;
    #$chr_d=~s/chrx/chr23/;
    #$chr_d=~s/chrY/chr24/;
    #$chr_d=~s/chry/chr24/;
    
    # need to delete
    $junc_type=1 if $junc_type!=-1;
    
    # change upstream/downstream site do breakpoint
    $rst_d = ($strand_d eq '+') ? --$rst_d : ++$rst_d; 
    $rst_a = ($strand_a  eq '+') ? ++$rst_a : --$rst_a;
    
    # change breakpoint strand: + means left or upstream of the breakpoint will be used in fusion, while - means right or downstream of breakpoint
    # no need to change strand_d
    $strand_a = ($strand_a eq '+') ? "-" : "+";
    
    # order breakpoint
    my ($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right)=($chr_d, $rst_d, $strand_d, $cigar_d, $chr_a, $rst_a, $strand_a, $cigar_a);
    if($chr_d eq $chr_a){
	($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right)=($chr_a, $rst_a, $strand_a, $cigar_a, $chr_d, $rst_d, $strand_d, $cigar_d) if $rst_a < $rst_d;
    }else{
	($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right)=($chr_a, $rst_a, $strand_a, $cigar_a, $chr_d, $rst_d, $strand_d, $cigar_d) if substr($chr_a, 3) < substr($chr_d, 3);
    }
    say join "\t", ($chr_left, $pos_left, $strand_left, $cigar_left, $chr_right, $pos_right, $strand_right, $cigar_right, $junc_type, $read_name);
}
close(IN);

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to format STAR chimeric input
Usage: perl $scriptName Chimeric.out.junction >output

    -h	--help				print this help information screen

HELP
    exit(-1);
}