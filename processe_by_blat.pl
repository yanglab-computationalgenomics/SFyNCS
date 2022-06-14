#!/usr/bin/env perl

# 2022-06-14

# 1. Function
# Add blat supported statistics
# 1.1. split read alignment should span breakpoints in the artifact reference 
# 1.2. for each junction, take upstream and downstream sequence to construct artifact reference, sequence length is defiend by $length_in_fusion and $length_out_fusion
# 1.3. for each read, at least align_percentage of read can be blat
# 1.4. for each read, select the best pslScore as alignment. If there are more than two alignments have best pslScore, this read will not support breakpoint
# 1.5. calculate sequence identity after taking --length_for_identity bp flanking sequence and align them with Needleman Wunsch algorithm 
# 1.6. calculate the distance between segment and artifact breakpoints 
# 1.7. scan for canonical split motif in region defined by --motif_searching_length

# 2. Input
# column 1: chromosome of the left segment
# column 2: left segment breakpoint
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment breakpoint
# column 6: strand of the right segment
# column 7: cluster id
# column 8: split read count (star)
# column 9: read pair count (star)
# column 10: split read count (processed by tophat)
# column 11: potential split read count (processed by tophat)
# column 12: read pair count (processed by tophat)
# column 13: minimum distance of read pair to left breakpoint 
# column 14: minimum distance of read pair to right breakpoint 
# column 15: split reads (star)
# column 16: read pairs (star)
# column 17: split reads (processed by tophat)
# column 18: potential split reads (processed by tophat)
# column 19: read pairs (processed by tophat)
# column 20: distance of read pair to left breakpoint
# column 21: distance of read pair to right breakpoint


# 3. Output
# column 1: chromosome of the left segment
# column 2: left segment site
# column 3: strand of the left segment
# column 4: chromosome of right segment
# column 5: right segment site
# column 6: strand of the right segment
# column 7: cluster id
# column 8: split read count (star)
# column 9: read pair count (star)
# column 10: split read count (processed by tophat)
# column 11: potential split read count (processed by tophat)
# column 12: read pair count (processed by tophat)
# column 13: split read count (blat tophat split reads and tophat potential split reads)
# column 14: split read count (blat tophat split reads)
# column 15: minimum distance of read pair to left breakpoint 
# column 16: minimum distance of read pair to right breakpoint 
# column 17: identity
# column 18: minimum blat distace of left breakpoint (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 19: minimum blat distace of right breakpoint (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 20: minimum blat distace of left breakpoint (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 21: minimum blat distace of right breakpoint (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 22: is canonical split site
# column 23: split reads (star)
# column 24: read pairs (star)
# column 25: split reads (processed by tophat)
# column 26: potential split reads (processed by tophat)
# column 27: read pairs (processed by tophat)
# column 28: split reads (blat tophat split read and tophat potential split read)
# column 29: distance of read pair to left breakpoint
# column 30: distance of read pair to right breakpoint
# column 31: distace of split read to left breakpoint when blating to artifact reference (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 32: distace of split read to right breakpoint when blating to artifact reference (blat align tophat split read and tophat potential split read, use mean if both read 1 and read 2 are split read)
# column 33: distace of split read to left breakpoint when blating to artifact reference (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 34: distace of split read to right breakpoint when blating to artifact reference (blat tophat split read, use mean if both read 1 and read 2 are split read)
# column 35: alignment of left segment
# column 36: alignment of right segment


use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $slop_length=5;
my $motif_searching_length=5;
my $length_in_fusion=1000000; # 1Mb
my $length_out_fusion=100;
my $length_for_identity=10;
my $fasta;
my $align_percentage=0.9; # blat alignment >=read_length*$align_percentage will be kept

GetOptions(
    's|slop_length=i'  			=> \$slop_length,
    'm|motif_searching_length=i'	=>\$motif_searching_length,
    'f|fasta=s'     			=> \$fasta,
    'a|align_percentage=f'     		=> \$align_percentage,
    'i|length_in_fusion=i'     		=> \$length_in_fusion,
    'o|length_out_fusion=i'     	=> \$length_out_fusion,
    'h|help'    			=> sub{usage()}
)||usage();

die "Please set -f or --fasta" if ! defined $fasta;
my $left_pos_in_artifact_seq=$length_in_fusion;
my $right_pos_in_artifact_seq=$length_in_fusion+2*$length_out_fusion+1;

say join "\t", ("Chr_left", "Pos_left", "Strand_left", "Chr_right", "Pos_right", "Strand_right", "Cluster_id",
    "Split_read_count_(star)", "Read_pair_count_(star)", "Split_read_count_(tophat)", "Potential_split_read_count_(tophat)", "Read_pair_count_(tophat)", "Split_read_count_(blat_tophat_split_and_tophat_potential_split_reads)", "Split_read_count_(blat_tophat_split_reads)",
    "Minimum_read_pair_distance_to_left", "Minimum_read_pair_distance_to_right",
    "Identity", "Minimum_blat_distance_to_left_(tophat_split_and_potential_split_reads)", "Minimum_blat_distance_to_right_(tophat_split_and_potential_split_reads)", "Minimum_blat_distance_to_left_(tophat_split_reads)", "Minimum_blat_distance_to_right_(tophat_split_reads)",
    "Is_canonical_motif",
    "Split_reads_(star)", "Read_pairs_(star)", "Split_reads_(tophat)", "Potential_split_reads_(tophat)", "Read_pairs_(tophat)", "Split_reads_(blat_tophat_split_and_tophat_potential_split_reads)",
    "Read_pair_distance_to_left", "Read_pair_distance_to_right",
    "Blat_distance_to_left_(tophat_split_and_tophat_potential_split_reads)", "Blat_distance_to_right_(tophat_split_and_tophat_potential_split_reads)", "Blat_distance_to_left_(tophat_split_reads)", "Blat_distance_to_right_(tophat_split_eads)",
    "Identity_align_left", "Identity_align_right");


open IN, $ARGV[0] or die "Can't open $ARGV[0]:$!";
# input have header
<IN>;
while(<IN>){
    chomp;
    my ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $cluster_id,
        $split_read_count_star, $read_pair_count_star, $split_read_count_tophat, $potential_split_read_count_tophat, $read_pair_count_tophat,
        $min_read_pair_distance_left, $min_read_pair_distance_right,
        $split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat,
        $read_pair_distance_to_left, $read_pair_distance_to_right)=split "\t", $_;
    my $split_reads_blat="NA";
    my $split_read_count_blat=0;
    
    my @Split_reads_tophat=split ",", $split_reads_tophat;
    my @Split_reads_poential_tophat=split ",", $potential_split_reads_tophat;
    
    my $artifact_seq_left=&getFlankingSequence($chr_left, $pos_left, $strand_left, $length_in_fusion, $length_out_fusion, "left");
    my $artifact_seq_right=&getFlankingSequence($chr_right, $pos_right, $strand_right, $length_in_fusion, $length_out_fusion, "right");
    
    my $artifact_seq=$artifact_seq_left.$artifact_seq_right;
    # system failed to ouput long $artifact_seq, use handle instead
    # system("echo '>artifact_ref' >temp_artifact.fa");
    # system("echo $artifact_seq >>temp_artifact.fa");
    open ARTIFACT, ">temp_artifact.fa";
    say ARTIFACT ">artifact_ref";
    say ARTIFACT $artifact_seq;
    close(ARTIFACT);
    
    my ($blat_distance_left_output, $blat_distance_right_output)=("NA") x 2;
    foreach my $read (@Split_reads_tophat, @Split_reads_poential_tophat){
        next if $read eq "NA";
        system("echo $read >temp_reads.tsv");
        system('grep -A 1 -f temp_reads.tsv selected_discordant_reads_1.fastq | grep -v "^--$" | sed "/^@/ s#^@#>#;/^>/ s#\$#_1#" >temp_reads.fa');
        system('grep -A 1 -f temp_reads.tsv selected_discordant_reads_2.fastq | grep -v "^--$" | sed "/^@/ s#^@#>#;/^>/ s#\$#_2#" >>temp_reads.fa');
        
        my @BLAT=`blat -t=dna -q=dna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -noHead temp_artifact.fa temp_reads.fa /dev/stdout | grep -v "^[a-zA-Z]"`;
        next if @BLAT==0;
        
        my $temp_processed_Alignment=&processAlignment(\@BLAT);
        my @Processed_Alignment=@{$temp_processed_Alignment};
	
	if(@Processed_Alignment==2 || @Processed_Alignment==3){ # one mate may align to two segments, so it would be 2 or 3
	    my ($name_1, $start_1, $end_1, $name_2, $start_2, $end_2);
	    my ($blockSizes_1, $tStarts_1, $blockSizes_2, $tStarts_2);
	    
	    if(@Processed_Alignment==2){
		($name_1, $start_1, $end_1)=($Processed_Alignment[0][7], $Processed_Alignment[0][3], $Processed_Alignment[0][4]);  # $qSize, $qStart, $qEnd, $tStart, $tEnd, $pslScore, $percentIdentity, $qName, $blockSize, $tStarts in Array
		($name_2, $start_2, $end_2)=($Processed_Alignment[1][7], $Processed_Alignment[1][3], $Processed_Alignment[1][4]);
		($blockSizes_1, $tStarts_1, $blockSizes_2, $tStarts_2)=($Processed_Alignment[0][8], $Processed_Alignment[0][9], $Processed_Alignment[1][8], $Processed_Alignment[1][9]);
		next if $name_1 eq $name_2; # only have one reads
	    }else{ # one mate align to different region
		my ($temp_name_1, $temp_start_1, $temp_end_1, $temp_blockSizes_1, $temp_tStarts_1)=($Processed_Alignment[0][7], $Processed_Alignment[0][3], $Processed_Alignment[0][4], $Processed_Alignment[0][8], $Processed_Alignment[0][9]);
		my ($temp_name_2, $temp_start_2, $temp_end_2, $temp_blockSizes_2, $temp_tStarts_2)=($Processed_Alignment[1][7], $Processed_Alignment[1][3], $Processed_Alignment[1][4], $Processed_Alignment[1][8], $Processed_Alignment[1][9]);
		my ($temp_name_3, $temp_start_3, $temp_end_3, $temp_blockSizes_3, $temp_tStarts_3)=($Processed_Alignment[2][7], $Processed_Alignment[2][3], $Processed_Alignment[2][4], $Processed_Alignment[2][8], $Processed_Alignment[2][9]);
		next if $temp_name_1 eq $temp_name_2 && $temp_name_1 eq $temp_name_3; # only have one reads
		my ($temp_multiple_start_1, $temp_multiple_end_1, $temp_multiple_start_2, $temp_multiple_end_2, $temp_multiple_name);
		my ($temp_multiple_blockSizes_1, $temp_multiple_tStarts_1, $temp_multiple_blockSizes_2, $temp_multiple_tStarts_2);
		if($temp_name_1 eq $temp_name_2){
		    ($name_2, $start_2, $end_2)=($temp_name_3, $temp_start_3, $temp_end_3);
		    $temp_multiple_name=$temp_name_1;
		    ($temp_multiple_start_1, $temp_multiple_end_1, $temp_multiple_start_2, $temp_multiple_end_2)=($temp_start_1, $temp_end_1, $temp_start_2, $temp_end_2);
		    ($blockSizes_2, $tStarts_2)=($temp_blockSizes_3, $temp_tStarts_3);
		    ($temp_multiple_blockSizes_1, $temp_multiple_tStarts_1, $temp_multiple_blockSizes_2, $temp_multiple_tStarts_2)=($temp_blockSizes_1, $temp_tStarts_1, $temp_blockSizes_2, $temp_tStarts_2);
		}elsif($temp_name_1 eq $temp_name_3){
		    ($name_2, $start_2, $end_2)=($temp_name_2, $temp_start_2, $temp_end_2);
		    $temp_multiple_name=$temp_name_1;
		    ($temp_multiple_start_1, $temp_multiple_end_1, $temp_multiple_start_2, $temp_multiple_end_2)=($temp_start_1, $temp_end_1, $temp_start_3, $temp_end_3);
		    ($blockSizes_2, $tStarts_2)=($temp_blockSizes_2, $temp_tStarts_2);
		    ($temp_multiple_blockSizes_1, $temp_multiple_tStarts_1, $temp_multiple_blockSizes_2, $temp_multiple_tStarts_2)=($temp_blockSizes_1, $temp_tStarts_1, $temp_blockSizes_3, $temp_tStarts_3);
		}else{ # $temp_name_2 eq $temp_name_3
		    ($name_2, $start_2, $end_2)=($temp_name_1, $temp_start_1, $temp_end_1);
		    $temp_multiple_name=$temp_name_2;
		    ($temp_multiple_start_1, $temp_multiple_end_1, $temp_multiple_start_2, $temp_multiple_end_2)=($temp_start_2, $temp_end_2, $temp_start_3, $temp_end_3);
		    ($blockSizes_2, $tStarts_2)=($temp_blockSizes_1, $temp_tStarts_1);
		    ($temp_multiple_blockSizes_1, $temp_multiple_tStarts_1, $temp_multiple_blockSizes_2, $temp_multiple_tStarts_2)=($temp_blockSizes_2, $temp_tStarts_2, $temp_blockSizes_3, $temp_tStarts_3);
		}
		
		my $is_multiple_aligned_read_support_fusion=0;
		#$is_multiple_aligned_read_support_fusion=1 if ($temp_multiple_end_1<=$left_pos_in_artifact_seq && $temp_start_2>=$right_pos_in_artifact_seq) || ($temp_multiple_end_2<=$left_pos_in_artifact_seq && $temp_start_1>=$right_pos_in_artifact_seq);
		$is_multiple_aligned_read_support_fusion=1 if ($temp_multiple_end_1<=$left_pos_in_artifact_seq && ($temp_multiple_start_2+1)>=$right_pos_in_artifact_seq) || ($temp_multiple_end_2<=$left_pos_in_artifact_seq && ($temp_multiple_start_1+1)>=$right_pos_in_artifact_seq);
		next if $is_multiple_aligned_read_support_fusion==0;
		($name_1, $start_1, $end_1)=($temp_multiple_name, $temp_multiple_start_1, $temp_multiple_end_1); # randomly select one alignment
		($blockSizes_1, $tStarts_1)=($temp_multiple_start_1, $temp_multiple_end_1);
	    }
	    
	    my ($is_support_split_1, $is_support_split_2)=(0) x 2;
	    $is_support_split_1=1 if $left_pos_in_artifact_seq>($start_1+1) && $right_pos_in_artifact_seq<$end_1;
	    $is_support_split_2=1 if $left_pos_in_artifact_seq>($start_2+1) && $right_pos_in_artifact_seq<$end_2;
	    if($is_support_split_1==1 || $is_support_split_2==1){
		$split_reads_blat.=",".$read;
		$split_read_count_blat++;
	    }
	    
	    # get distance to fusion breakpoint, can include this to above if($is_support_split_1==1 || $is_support_split_2==1) later 
	    my ($blat_distance_left_1, $blat_distance_right_1)=&getDistanceInBlat($blockSizes_1, $tStarts_1);
	    my ($blat_distance_left_2, $blat_distance_right_2)=&getDistanceInBlat($blockSizes_2, $tStarts_2);
	    
	    my ($temp_blat_distance_left, $temp_blat_distance_right)=("NA") x 2;
	    if($blat_distance_left_1 ne "NA" && $blat_distance_left_2 ne "NA"){ # take average if both are segment
		$temp_blat_distance_left=sprintf("%.2f", (abs($blat_distance_left_1)+abs($blat_distance_left_2))/2);
		$temp_blat_distance_right=sprintf("%.2f", (abs($blat_distance_right_1)+abs($blat_distance_right_2))/2);
	    }elsif($blat_distance_left_1 ne "NA"){
		($temp_blat_distance_left, $temp_blat_distance_right)=($blat_distance_left_1, $blat_distance_right_1);
	    }elsif($blat_distance_left_2 ne "NA"){
		($temp_blat_distance_left, $temp_blat_distance_right)=($blat_distance_left_2, $blat_distance_right_2);
	    }
	    
	    $blat_distance_left_output.=",".$temp_blat_distance_left if $is_support_split_1==1 || $is_support_split_2==1;
	    $blat_distance_right_output.=",".$temp_blat_distance_right if $is_support_split_1==1 || $is_support_split_2==1;
	}
    }
    
    # get identity
    my $identity_seq_left=&getFlankingSequence($chr_left, $pos_left, $strand_left, $length_for_identity, $length_for_identity, "left");
    my $identity_seq_right=&getFlankingSequence($chr_right, $pos_right, $strand_right, $length_for_identity, $length_for_identity, "right");
    my ($identity_output, $identity_align_left_output, $identity_align_right_output)=&getIdentity($identity_seq_left, $identity_seq_right);

    # is contain canonical split motif
    my $flanking_search_sequence_left=$identity_seq_left;
    $flanking_search_sequence_left=substr($flanking_search_sequence_left, $length_for_identity-$motif_searching_length, 2*$motif_searching_length);
    my $is_left_contain_donor_motif=($flanking_search_sequence_left=~/GT/) ? "Y" : "N";
    $flanking_search_sequence_left=~tr/[ACTG]/[TGAC]/;
    $flanking_search_sequence_left=reverse($flanking_search_sequence_left);
    my $is_left_contain_acceptor_motif=($flanking_search_sequence_left=~/[CTA]AG/) ? "Y" : "N";
    
    my $flanking_search_sequence_right=$identity_seq_right;
    $flanking_search_sequence_right=substr($flanking_search_sequence_right, $length_for_identity-$motif_searching_length, 2*$motif_searching_length);
    my $is_right_contain_acceptor_motif=($flanking_search_sequence_right=~/[CTA]AG/) ? "Y" : "N";
    $flanking_search_sequence_right=~tr/[ACTG]/[TGAC]/;
    $flanking_search_sequence_right=reverse($flanking_search_sequence_right);
    my $is_right_contain_donor_motif=($flanking_search_sequence_right=~/GT/) ? "Y" : "N";
    
    my $is_split_site_contain_canonical_motif="N";
    $is_split_site_contain_canonical_motif="Y" if ($is_left_contain_donor_motif eq "Y" && $is_right_contain_acceptor_motif eq "Y") ||
	($is_right_contain_donor_motif eq "Y" && $is_left_contain_acceptor_motif eq "Y");
    
    
    $split_reads_blat=~s/NA,//;
    $blat_distance_left_output=~s/NA,//;
    $blat_distance_right_output=~s/NA,//;
    # minimum distance
    my ($minimum_blat_distance_left, $minimum_blat_distance_right, $sum_minimum_blat_distance)=("NA") x 3;
    my @Blat_distance_left=split ",", $blat_distance_left_output;
    my @Blat_distance_right=split ",", $blat_distance_right_output;
    for(my $i=0; $i<@Blat_distance_left; $i++){
        next if $Blat_distance_left[$i] eq "NA" || $Blat_distance_right[$i] eq "NA";
        my $temp_sum_blat_distance=abs($Blat_distance_left[$i]+$Blat_distance_right[$i]);
        if($sum_minimum_blat_distance eq "NA" || $sum_minimum_blat_distance>$temp_sum_blat_distance){
            ($minimum_blat_distance_left, $minimum_blat_distance_right)=(abs($Blat_distance_left[$i]), abs($Blat_distance_right[$i]));
            $sum_minimum_blat_distance=$temp_sum_blat_distance;
        }
    }
    
    # statistics base on tophat split reads
    my %Temp_tophat_split_reads;
    foreach my $read (split ",", $split_reads_tophat){
	next if $read eq "NA";
	$Temp_tophat_split_reads{$read}=0;
    }
    my $split_read_count_blat_base_on_tophat_split_read=0;
    my ($minimum_blat_distance_base_on_tophat_split_read_left, $minimum_blat_distance_base_on_tophat_split_read_right, $sum_minimum_blat_base_on_tophat_split_read_distance)=("NA") x 3;
    my ($blat_distance_base_on_tophat_split_read_left, $blat_distance_base_on_tophat_split_read_right)=("NA") x 2;
    my @Split_reads_blat=split ",", $split_reads_blat;
    for(my $i=0; $i<@Blat_distance_left; $i++){
        next if $Blat_distance_left[$i] eq "NA" || $Blat_distance_right[$i] eq "NA";
	my $read=$Split_reads_blat[$i];
	next if $read eq "NA" || !exists $Temp_tophat_split_reads{$read};
	$split_read_count_blat_base_on_tophat_split_read++;
	$blat_distance_base_on_tophat_split_read_left.=",".$Blat_distance_left[$i];
	$blat_distance_base_on_tophat_split_read_right.=",".$Blat_distance_right[$i];
        my $temp_sum_blat_distance=abs($Blat_distance_left[$i]+$Blat_distance_right[$i]);
        if($sum_minimum_blat_base_on_tophat_split_read_distance eq "NA" || $sum_minimum_blat_base_on_tophat_split_read_distance>$temp_sum_blat_distance){
            ($minimum_blat_distance_base_on_tophat_split_read_left, $minimum_blat_distance_base_on_tophat_split_read_right)=(abs($Blat_distance_left[$i]), abs($Blat_distance_right[$i]));
            $sum_minimum_blat_base_on_tophat_split_read_distance=$temp_sum_blat_distance;
        }
    }
    $blat_distance_base_on_tophat_split_read_left=~s/NA,//;
    $blat_distance_base_on_tophat_split_read_right=~s/NA,//;
    
    
    say join "\t", ($chr_left, $pos_left, $strand_left, $chr_right, $pos_right, $strand_right, $cluster_id,
        $split_read_count_star, $read_pair_count_star, $split_read_count_tophat, $potential_split_read_count_tophat, $read_pair_count_tophat, $split_read_count_blat, $split_read_count_blat_base_on_tophat_split_read,
        $min_read_pair_distance_left, $min_read_pair_distance_right, $identity_output,
	$minimum_blat_distance_left, $minimum_blat_distance_right, $minimum_blat_distance_base_on_tophat_split_read_left, $minimum_blat_distance_base_on_tophat_split_read_right,
        $is_split_site_contain_canonical_motif,
	$split_reads_star, $read_pairs_star, $split_reads_tophat, $potential_split_reads_tophat, $read_pairs_tophat, $split_reads_blat,
        $read_pair_distance_to_left, $read_pair_distance_to_right, 
	$blat_distance_left_output, $blat_distance_right_output, $blat_distance_base_on_tophat_split_read_left, $blat_distance_base_on_tophat_split_read_right,
	$identity_align_left_output, $identity_align_right_output);
}
close(IN);
system('rm -f temp_artifact.fa temp_reads.tsv temp_reads.fa');

sub pslCalcMilliBad{
    my ($sizeMul, $qEnd, $qStart, $tEnd, $tStart, $qNumInsert, $tNumInsert, $matches, $repMatches, $misMatches, $isMrna) = @_;
    my $milliBad = 0;
    my $qAliSize = $sizeMul * ($qEnd - $qStart);
    my $tAliSize = $tEnd - $tStart;
    my $aliSize = $qAliSize;
    $aliSize = $tAliSize if ($tAliSize < $qAliSize);
    if($aliSize<=0){
        return $milliBad;
    }
    my $sizeDif = $qAliSize - $tAliSize;
    if($sizeDif<0){
        if($isMrna){
            $sizeDif=0;
        }else{
            $sizeDif=-$sizeDif;
        }
    }
    my $insertFactor = $qNumInsert;
    if(0==$isMrna){
        $insertFactor+=$tNumInsert;
    }
    my $total=($sizeMul * ($matches + $repMatches + $misMatches));
    if($total!=0){
        my $roundAwayFromZero = 3*log(1+$sizeDif);
        if($roundAwayFromZero<0){
            $roundAwayFromZero = int($roundAwayFromZero - 0.5);
        }else{
            $roundAwayFromZero = int($roundAwayFromZero + 0.5);
        }
        $milliBad = (1000 * ($misMatches*$sizeMul + $insertFactor + $roundAwayFromZero)) / $total;
    }
    return $milliBad;
} # sub pslCalcMilliBad()

sub processAlignment{
    my ($input_variable)=@_;
    my @readAlign=@{$input_variable};
    my (@Processed_BLAT, @Compare_BLAT);
    my $sizeMul=1; # used in processing blat alignment
    my $pre_qName="NA"; 
    
    foreach my $alignment (@readAlign){
        chomp;
        my ($matches, $misMatches, $repMatches, $nCount, $qNumInsert, $qBaseInsert, $tNumInsert, $tBaseInsert, $strand,
            $qName, $qSize, $qStart, $qEnd, $tName, $tSize, $tStart, $tEnd, $blockCount, $blockSizes, $qStarts, $tStarts) = split('\t', $alignment);
        next if ($qEnd-$qStart)<$qSize*$align_percentage; # filter by align length
        
        my $pslScore = $sizeMul * ($matches + ( $repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumInsert;
        my $milliBad = int(&pslCalcMilliBad($sizeMul, $qEnd, $qStart, $tEnd, $tStart, $qNumInsert, $tNumInsert, $matches, $repMatches, $misMatches, 1));
        my $percentIdentity = 100.0 - $milliBad * 0.1;
        
        if($qName eq $pre_qName){
            push @Compare_BLAT, [$qSize, $qStart, $qEnd, $tStart, $tEnd, $pslScore, $percentIdentity, $blockSizes, $tStarts];
        }else{
            if($pre_qName ne "NA"){
                my $processed_alignment=&compareAlignment(\@Compare_BLAT) if $pre_qName ne "NA";
                foreach my $processed_line (@{$processed_alignment}){
                    my ($temp_qSize, $temp_qStart, $temp_qEnd, $temp_tStart, $temp_tEnd, $temp_pslScore, $temp_percentIdentity, $temp_blockSizes, $temp_tStarts)=@{$processed_line};
                    push @Processed_BLAT, [$temp_qSize, $temp_qStart, $temp_qEnd, $temp_tStart, $temp_tEnd, $temp_pslScore, $temp_percentIdentity, $pre_qName, $temp_blockSizes, $temp_tStarts];
                }
                @Compare_BLAT=();
            }
            $pre_qName=$qName;
            push @Compare_BLAT, [$qSize, $qStart, $qEnd, $tStart, $tEnd, $pslScore, $percentIdentity, $blockSizes, $tStarts];
        }
    }
    if($pre_qName ne "NA"){
        my $processed_alignment=&compareAlignment(\@Compare_BLAT);
        foreach my $processed_line (@{$processed_alignment}){
            my ($temp_qSize, $temp_qStart, $temp_qEnd, $temp_tStart, $temp_tEnd, $temp_pslScore, $temp_percentIdentity, $temp_blockSizes, $temp_tStarts)=@{$processed_line};
            push @Processed_BLAT, [$temp_qSize, $temp_qStart, $temp_qEnd, $temp_tStart, $temp_tEnd, $temp_pslScore, $temp_percentIdentity, $pre_qName, $temp_blockSizes, $temp_tStarts];
        }
    }
    
    return(\@Processed_BLAT);
}


sub compareAlignment{
    my ($input_variable)=@_;
    my @readAlign2=@{$input_variable};
    my @Select_alignment;
    
    @readAlign2=sort {$b->[5] <=> $a->[5] || $b->[6] <=> $a->[6]} @readAlign2;
    my ($topQSize, $topQStart, $topQEnd, $topTStart, $topTEnd, $topPslScore, $topPercentIdentity, $topBlockSizes, $topTStarts)=@{shift @readAlign2};
    push @Select_alignment, [$topQSize, $topQStart, $topQEnd, $topTStart, $topTEnd, $topPslScore, $topPercentIdentity, $topBlockSizes, $topTStarts];
    foreach my $processed_line (@readAlign2){
        my ($temp_qSize, $temp_qStart, $temp_qEnd, $temp_tStart, $temp_tEnd, $temp_pslScore, $temp_percentIdentity, $temp_blockSizes, $temp_tStarts)=@{$processed_line};
        next if $temp_pslScore<$topPslScore;
        push @Select_alignment, [$temp_qSize, $temp_qStart, $temp_qEnd, $temp_tStart, $temp_tEnd, $temp_pslScore, $temp_percentIdentity, $temp_blockSizes, $temp_tStarts];
    }
    return(\@Select_alignment);
}


sub getDistanceInBlat{
    my ($blockSizes, $tStarts)=@_;
    my ($distance_left, $distance_right)=("NA") x 2;
    if($blockSizes ne "NA"){
	EXIT_get_distance_block:{
	    my @BlockSizes=split ",", $blockSizes;
	    my @TStarts=split ",", $tStarts;
	    my $blockCount=@BlockSizes;
	    my $start=$TStarts[0]+1; # change to 1-base
	    my $end=$TStarts[$blockCount-1]+$BlockSizes[$blockCount-1];
	    last EXIT_get_distance_block if $end<=$right_pos_in_artifact_seq || $start>=$left_pos_in_artifact_seq;
	    
	    for(my $i=0; $i<$blockCount-1; $i++){
		next if $TStarts[$i+1]+1<$left_pos_in_artifact_seq; # +1 change to 1-base # $i will indicates the last left segment, $i+1 indicates the first right segment
		$distance_left=$TStarts[$i]+$BlockSizes[$i]-$left_pos_in_artifact_seq;
		$distance_right=$TStarts[$i+1]+1-$right_pos_in_artifact_seq;
		last;
	    }
	}
    }
    return ($distance_left, $distance_right);
}


sub getFlankingSequence{
    my ($chr, $pos, $strand, $length_in_fusion, $length_out_fusion, $left_or_right_segment)=@_;
    my ($start, $end);
    my $append_N_length=0; # if $start < 1, append N
    if($strand eq '+'){
        $start=$pos-$length_in_fusion+1;
        if($start<1){
            $append_N_length=abs($start)+1;
            $start=1;
        }
        $end=$pos+$length_out_fusion; # end bigger than chr size will generate short sequence
    }else{
        $start=$pos-$length_out_fusion;
        if($start<1){
            $append_N_length=abs($start)+1;
            $start=1;
        }
        $end=$pos+$length_in_fusion-1; # end bigger than chr size will generate short sequence
    }
    
    my $flankingSeq;
    my $artifact_coordinate=$chr.":".$start."-".$end;
    if($left_or_right_segment eq "left"){
	if($strand eq '+'){
	    $flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}'`;
	}else{
	    $flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}' | rev | tr 'ACGT' 'TGCA'`;
	}
    }else{
	if($strand eq '+'){
	    $flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}' | rev | tr 'ACGT' 'TGCA'`;
	}else{
	    $flankingSeq=`samtools faidx $fasta $artifact_coordinate | awk 'NR==2{output=\$0;} NR>2{output=output""\$0;} END{print output;}'`;
	}
    }
    chomp($flankingSeq);
    
    my $length_artifact_fa=$length_in_fusion+$length_out_fusion;
    if($left_or_right_segment eq "left"){
	if($strand eq '+'){
	    if($append_N_length>0){
		my $temp_N_string="N" x $append_N_length;
		$flankingSeq=$temp_N_string.$flankingSeq;
	    }
	    my $length_flankingSeq=length($flankingSeq);
	    if($length_flankingSeq<$length_artifact_fa){
		my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
		$flankingSeq=$flankingSeq.$temp_N_string;
	    }
	}else{
	    if($append_N_length>0){
		my $temp_N_string="N" x $append_N_length;
		$flankingSeq=$flankingSeq.$temp_N_string;
	    }
	    my $length_flankingSeq=length($flankingSeq);
	    if($length_flankingSeq<$length_artifact_fa){
		my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
		$flankingSeq=$temp_N_string.$flankingSeq;
	    }
	}
    }else{
	if($strand eq '+'){
	    if($append_N_length>0){
		my $temp_N_string="N" x $append_N_length;
		$flankingSeq=$flankingSeq.$temp_N_string;
	    }
	    my $length_flankingSeq=length($flankingSeq);
	    if($length_flankingSeq<$length_artifact_fa){
		my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
		$flankingSeq=$temp_N_string.$flankingSeq;
	    }
	}else{
	    if($append_N_length>0){
		my $temp_N_string="N" x $append_N_length;
		$flankingSeq=$temp_N_string.$flankingSeq;
	    }
	    my $length_flankingSeq=length($flankingSeq);
	    if($length_flankingSeq<$length_artifact_fa){
		my $temp_N_string="N" x ($length_artifact_fa-$length_flankingSeq);
		$flankingSeq=$flankingSeq.$temp_N_string;
	    }
	}
    }
    
    return($flankingSeq);
}


sub getIdentity{
    my ($seq1, $seq2)=@_;

    # scoring scheme
    my $MATCH=1; # +1 for letters that match
    my $MISMATCH=-1; # -1 for letters that mismatch
    my $GAP=-1; # -1 for any gap
    
    # initialization
    my @matrix;
    $matrix[0][0]{score}=0;
    $matrix[0][0]{pointer}="none";
    
    
    for(my $j=1; $j<=length($seq1); $j++) {
	$matrix[0][$j]{score}=$GAP*$j;
	$matrix[0][$j]{pointer}="left";
    }
    for(my $i=1; $i<=length($seq2); $i++) {
	$matrix[$i][0]{score}=$GAP*$i;
	$matrix[$i][0]{pointer}="up";
    }
    
    # fill
    for(my $i=1; $i<=length($seq2); $i++) {
	for(my $j=1; $j<=length($seq1); $j++) {
	    my ($diagonal_score, $left_score, $up_score);
    
	    # calculate match score
	    my $letter1=substr($seq1, $j-1, 1);
	    my $letter2=substr($seq2, $i-1, 1);                         
	    
	    if($letter1 eq $letter2) {
		$diagonal_score=$matrix[$i-1][$j-1]{score} + $MATCH;
	    }
	    else{
		$diagonal_score=$matrix[$i-1][$j-1]{score} + $MISMATCH;
	    }
    
	    # calculate gap scores
	    $up_score=$matrix[$i-1][$j]{score}+$GAP;
	    $left_score=$matrix[$i][$j-1]{score}+$GAP;
    
	    # choose best score
	    if($diagonal_score>=$up_score) {
		if($diagonal_score>=$left_score) {
		    $matrix[$i][$j]{score}=$diagonal_score;
		    $matrix[$i][$j]{pointer}="diagonal";
		}
	    else{
		    $matrix[$i][$j]{score}=$left_score;
		    $matrix[$i][$j]{pointer}="left";
		}
	    }else{
		if($up_score>=$left_score) {
		    $matrix[$i][$j]{score}=$up_score;
		    $matrix[$i][$j]{pointer}="up";
		}
		else{
		    $matrix[$i][$j]{score}=$left_score;
		    $matrix[$i][$j]{pointer}="left";
		}
	    }
	}
    }
    # trace-back
    
    my $align1="";
    my $align2="";
    
    # start at last cell of matrix
    my $j=length($seq1);
    my $i=length($seq2);
    
    while(1){
	last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix
    
	if($matrix[$i][$j]{pointer} eq "diagonal"){
	    $align1.=substr($seq1, $j-1, 1);
	    $align2.=substr($seq2, $i-1, 1);
	    $i--;
	    $j--;
	}elsif($matrix[$i][$j]{pointer} eq "left"){
	    $align1.=substr($seq1, $j-1, 1);
	    $align2.="-";
	    $j--;
	}elsif($matrix[$i][$j]{pointer} eq "up") {
	    $align1.="-";
	    $align2.=substr($seq2, $i-1, 1);
	    $i--;
	}    
    }
    $align1=reverse $align1;
    $align2=reverse $align2;
    
    # use blast identity here, refer to https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity for more
    my $align_length=length($align1);
    my @Align_1=split "", $align1;
    my @Align_2=split "", $align2;
    my $align_match=0;
    for(my $i=0; $i<$align_length; $i++){
	$align_match++ if $Align_1[$i] eq $Align_2[$i];
    }
    my $identity=sprintf("%.2f", $align_match/$align_length);
    return ($identity, $align1, $align2);
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
This script was used to get blat supported statistics
Usage: perl $scriptName input >output
Options:

    -s --slop_length		Slop length, read align within this slop_length of breakpoint will support the breakpoint [default: $slop_length]
    -m --motif_searching_length	Searching motif in upstream and downstream sequence of breakpoint [default: $motif_searching_length]
    -f --fasta    		reference fasta file (must be given)
    -a --align_percentage   	at least align_percentage of read should be aligned by blat [default: $align_percentage]
    -i --length_in_fusion   	take this value bp flanking breakpoint in fusion region to contructe artifact reference [default: $length_in_fusion]
    -o --length_out_fusion   	take this value bp flanking breakpoint out fusion region to contructe artifact reference[default: $length_out_fusion]
    -h --help     		print this help information
HELP
    exit(-1);
}


