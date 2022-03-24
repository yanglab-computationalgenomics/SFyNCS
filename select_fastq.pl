#!/usr/bin/env perl

# 2022-03-15

# 1. Function
# Select potential candidates' fastq
# 1.1. will delete strings after space if read name contain space

# 2. Input
# 2.1. potential candiates' read name
# 2.2. fastq

# 3. Output
# Potential candidates' fastq


use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $have_space;
GetOptions(
    's|have_space'  => \$have_space,
    'h|help'	    => sub{usage()}
)||usage();

my %hash;
open IN, $ARGV[0];
while(<IN>){
    chomp;
    $hash{$_}=0;
}
close(IN);

open IN, $ARGV[1];
while(<IN>){
    chomp;
    my $name=$_;
    my $name2=$name;
    $name2=~s/^@//;
    $name2=~s/\s.*// if defined $have_space;
    my $name3=$name;
    $name3=~s/\s.*// if defined $have_space;
    chomp(my $seq=<IN>);
    chomp(my $name_2=<IN>);
    chomp(my $qual=<IN>);
    if(exists $hash{$name2}){
        say $name3;
        say $seq;
        say $name_2;
        say $qual;
    }
}
close(IN);

open IN, $ARGV[2];
while(<IN>){
    chomp;
    my $name=$_;
    my $name2=$name;
    $name2=~s/^@//;
    $name2=~s/\s.*// if defined $have_space;
    my $name3=$name;
    $name3=~s/\s.*// if defined $have_space;
    chomp(my $seq=<IN>);
    chomp(my $name_2=<IN>);
    chomp(my $qual=<IN>);
    if(exists $hash{$name2}){
	say STDERR $name3;
        say STDERR $seq;
        say STDERR $name_2;
        say STDERR $qual;
    }
}
close(IN);


sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName reads_name.tsv 1.fastq 2.fastq >output_1.fastq 2>output_2.fastq
Options:

    -h --help		print this help information
    -s --have_space 	read name have space. String after space will be deleted in STAR, use this parameter to get discordant reads
HELP
    exit(-1);
}
