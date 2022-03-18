#!/bin/bash

# 2022-03-17

usage(){
 cat << EOF
Description:
    This script was used to identify fusion transcripts from pair-end RNA-seq data.
        
Usage:
    $0 [-p 1 -o output] -a annotation_file -g genome_fasta -s star_index -t tophat_index read_1.fastq read_2.fastq
    
Options:
    -a Gene annotation file.
    -g Genome fasta.
    -o Output directory [default: current directory].
    -p Number of threads [default: 1].
    -s Path to the genome directory where STAR genome indices where generated.
    -t Path to Tophat index.
    -h Print this help menu.

Gene annotation format (must have header, start and end are 1-base):
    Chr     Start       End         Strand      Symbol      Gene_type
    chr7    55019021    55211628    +           EGFR        protein_coding
EOF
    exit 0
}

[ $# -eq 0 ] && usage

thread_number=1
output_directory=$PWD

while getopts "a:g:o:p:s:t:h" OPTION
 do
  case $OPTION in
        a) annotation_file=$OPTARG;;
        g) genome_fasta=$OPTARG;;
        o) output_directory=$OPTARG;;
        p) thread_number=$OPTARG;;
        s) star_index=$OPTARG;;
        t) tophat_index=$OPTARG;;
        h) usage;;
  esac
 done
shift $((OPTIND - 1))

[ -z $annotation_file ] && echo "Please provide gene annotation with -a" && usage
[ -z $genome_fasta ] && echo "Please provide genome fasta with -g" && usage
[ -z $star_index ] && echo "Please provide STAR index with -s" && usage 
[ -z $tophat_index ] && echo "Please provide Tophat index with -t" && usage

function runningTime(){
  startS=$(echo $1 | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }') 
  endS=$(echo $2 | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }') 
  runTime=$(( $endS - $startS ))
  date -d@${runTime} -u '+%H:%M:%S'
}

startTime=$( date +%H:%M:%S)
echo -e "START:\t$startTime"


# Generate output directory
[ -e ${output_directory}/temp_output_SFyNCS ] && echo "Have ${output_directory}/temp_output_SFyNCS, please rename or delete the directory" && exit
mkdir ${output_directory}/temp_output_SFyNCS && cd ${output_directory}/temp_output_SFyNCS
fastq_1=$(readlink -f $1)
fastq_2=$(readlink -f $2)
ln -s $fastq_1 1.fastq
ln -s $fastq_2 2.fastq


# 1. Align fastq with STAR
echo -e "Step 1: Getting discordant information file with STAR"
STAR --genomeDir $star_index  \
    --readFilesIn 1.fastq 2.fastq  \
    --runThreadN $thread_number \
    --outFileNamePrefix star_output \
    --outReadsUnmapped None  \
    --twopassMode Basic \    
    --outSAMstrandField intronMotif  \
    --outSAMunmapped Within  \
    --chimSegmentMin 12  \
    --chimJunctionOverhangMin 12  \
    --chimOutJunctionFormat 1  \
    --alignSJDBoverhangMin 10  \
    --alignMatesGapMax 100000  \
    --alignIntronMax 100000  \
    --alignSJstitchMismatchNmax 5 -1 5 5  \
    --outSAMattrRGline ID:GRPundef  \
    --chimMultimapScoreRange 10 \
    --chimMultimapNmax 10 \
    --chimNonchimScoreDropMin 10  \
    --peOverlapNbasesMin 12 \
    --peOverlapMMp 0.1  \
    --outSAMtype BAM Unsorted  \
    --genomeLoad NoSharedMemory

discordant_count=$(grep -v "^#" star_output/Chimeric.out.junction | wc -l)
if [ $discordant_count -eq 0 ]; then
  echo -e "Chr_left\tPos_left\tStrand_left\tChr_right\tPos_right\tStrand_right\t \
     Split_read_count_(Tophat_and_Blat)\tRead_pair_count_(Tophat)\t \
     Minimum_read_distance_to_left\tMinimum_read_distance_to_right\t \
     Identity\tMinimum_blat_distance_to_left\tMinimum_blat_distance_to_right\t \
     Total_clusters_(left)\tTotal_clusters_(right)\tTotal_clusters_(merge)\t \
     (discordant_reads)%_support_fusion\t 
     SD_(discordant_reads)%_in_clusters_(left)\tSD_(discordant_reads)%_in_clusters_(right)\t \
     Fusion_annotations\tSplit_reads_(tophat_and_blat)\tRead_pairs_(tophat)" | sed "s# ##g" >../fusions.tsv
  echo "Finished!" && cd .. && rm -rf ${output_directory}/temp_output_SFyNCS && exit
fi


# 2. Format and generate preliminary fusions
echo "Step 2: Generating preliminary fusions"
toolDir=$(readlink $0 | xargs dirname)
perl $toolDir/format_STAR_chimeric_file.pl  star_output/Chimeric.out.junction >format_chimeric.tsv
perl $toolDir/remove_multiple_mapped_reads.pl format_chimeric.tsv >no_multiple-mapped.tsv
perl $toolDir/remove_duplicate_reads.pl no_multiple-mapped.tsv | sort | uniq >no_duplicate.tsv
perl $toolDir/merge_ajacent_breakpoints.pl no_duplicate.tsv >merge_adjacent.tsv
perl $toolDir/cluster_discordant_reads.pl merge_adjacent.tsv >cluster.tsv
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2-1,$2,$3,$7,$8,$9; print $4,$5-1,$5,$6,$7,$8,$9;}' cluster.tsv | sort -k1,1 -k2,2n | uniq >temp_cluster.bed
perl $toolDir/identify_fusion_candidates_from_cluster_reads.pl cluster.tsv >preliminary_candidates.tsv
# split_read >=1 && read_pair>=1 && split_read+read_pair>=3, remove chrM related fusions
awk 'NR==1 || ($9>0 && $10>0 && $9+$10>=3)' preliminary_candidates.tsv | grep -v chrM >temp.tsv
mv temp.tsv preliminary_candidates.tsv
rm cluster.tsv format_chimeric.tsv merge_adjacent.tsv no_duplicate.tsv no_multiple-mapped.tsv


# 3. process preliminary fusions with tophat
echo "Step 3: Processing with Tophat"
cut -f11,12 preliminary_candidates.tsv | sed "s#\t#\n#;s#,#\n#g" | grep -vw "NA" | grep -vP "Split_reads|Read_pairs" | sort | uniq >selected_discordant_reads.tsv
perl $toolDir/select_fastq.pl selected_discordant_reads.tsv 1.fastq 2.fastq >selected_discordant_reads_1.fastq 2>selected_discordant_reads_2.fastq
rm selected_discordant_reads.tsv

tophat --no-coverage-search \
    --fusion-search \
    --fusion-anchor-length 12 \
    --fusion-min-dist 100000 \
    --read-mismatches 4 \
    --read-gap-length 4 \
    --read-edit-dist 4 \
    --splice-mismatches 2 \
    --max-insertion-length 4 \
    --max-deletion-length4 \
    --segment-mismatches 3 \
    --fusion-read-mismatches 4 \
    -o tophat_out \
    -p $thread_number \
    $tophat_index \
    selected_discordant_reads_1.fastq selected_discordant_reads_2.fastq
ln -s tophat_out/accepted_hits.bam tophat.bam
samtools index tophat.bam
perl $toolDir/processe_read_supported_by_tophat.pl tophat.bam preliminary_candidates.tsv >processed_with_tophat.tsv
# split_read_tophat >=1 && read_pair_tophat>=1 &&  split_read_tophat+read_pair_tophat>=3
awk 'NR==1 || ($8>0 && $9>0 && $8+$9>=3)' processed_with_tophat.tsv >temp.tsv
mv temp.tsv processed_with_tophat.tsv
rm -rf tophat* preliminary_candidates.tsv 


# 4. process preliminary fusions with blat
echo "Step 4: Processing with Blat"
# split file to 50 chunk to redunce memory and speed up processing
lineCount=$( wc -l processed_with_tophat.tsv | cut -f1 -d " " )
chunk=$( echo "$lineCount/50" | bc )
for((i=1; i<50; i++))
  do
    mkdir chunk_${i}
    head -n 1 processed_with_tophat.tsv >chunk_${i}/processed_with_tophat.tsv
    start=$[ ($i-1)*$chunk+1 ]
    end=$[ $i*chunk ]
    sed -n "${start},${end}p" processed_with_tophat.tsv >>chunk_${i}/processed_with_tophat.tsv
  done
sed -i '1d' chunk_1/processed_with_tophat.tsv
start=$[ 49*$chunk+1 ]
end=$lineCount
mkdir chunk_50
head -n 1 processed_with_tophat.tsv >chunk_50/processed_with_tophat.tsv
sed -n "${start},${end}p" processed_with_tophat.tsv >>chunk_50/processed_with_tophat.tsv

for i in {1..50}
  do
    cd chunk_${i}
    ln -s ../selected_discordant_reads_1.fastq .
    ln -s ../selected_discordant_reads_1.fastq .
    perl $toolDir/processe_read_supported_by_blat.pl -f $genome_fasta processed_with_tophat.tsv >processed_with_blat.tsv
    # split_read_blat >=1 && split_read_blat+read_pair_tophat>=3 && (already fulfillment: read_pair_tophat>=1)
    awk 'NR==1 || ($8>0 && $8+$9>=3)' processed_with_blat.tsv >temp.tsv
    mv temp.tsv processed_with_blat.tsv
    cd ..
  done

cp chunk_1/processed_with_blat.tsv .
for((i=2; i<51; i++))
  do
    sed -n '2,$p' chunk_${i}/processed_with_blat.tsv >>processed_with_blat.tsv
  done
rm -rf chunk*
rm selected_discordant_reads_1.fastq selected_discordant_reads_2.fastq


# 5. generating fusion statistics
echo "Step 5: Generating fusion statistics"
# split file to 50 chunk to redunce memory and speed up processing
lineCount=$( wc -l processed_with_blat.tsv | cut -f1 -d " " )
chunk=$( echo "$lineCount/50" | bc )
for((i=1; i<50; i++))
  do
    mkdir chunk_${i}
    head -n 1 processed_with_blat.tsv >chunk_${i}/processed_with_blat.tsv
    start=$[ ($i-1)*$chunk+1 ]
    end=$[ $i*chunk ]
    sed -n "${start},${end}p" processed_with_blat.tsv >>chunk_${i}/processed_with_blat.tsv
  done
sed -i '1d' chunk_1/processed_with_blat.tsv
start=$[ 49*$chunk+1 ]
end=$lineCount
mkdir chunk_50
head -n 1 processed_with_blat.tsv >chunk_50/processed_with_blat.tsv
sed -n "${start},${end}p" processed_with_blat.tsv >>chunk_50/processed_with_blat.tsv

for i in {1..50}
  do
    cd chunk_${i}
    ln -s ../temp_cluster.bed .
    perl $toolDir/cluster_statistics.pl processed_with_blat.tsv >fusion_statistics.tsv
    cd ..
  done

cp chunk_1/fusion_statistics.tsv .
for((i=2; i<51; i++))
  do
    sed -n '2,$p' chunk_${i}/fusion_statistics.tsv >>fusion_statistics.tsv
  done
rm -rf chunk* temp_cluster.bed


# 6. annotating and getting final fusions
echo "Step 6: Generating final fusions"
# fusion_distance>=500kbp && read_pair_distance<=500bp && identity<=0.8 && sd>=0.1 && (already fulfillment: split_read_blat >=1 && read_pair_tophat>=1 && split_read_blat+read_pair_tophat>=3)
awk -v fusion_distance=500000 'NR==1 || ($9<=500 && $10<=500 && $11<=0.8 && $18>=0.1 && $19>=0.1 && (($1!=$4) || ($1==$4 && $5-$2>=fusion_distance)))' fusion_statistics.tsv >filtered_fusions_statistics.tsv
perl $toolDir/annotate_fusions.pl $annotation_file filtered_fusions_statistics.tsv >../fusions.tsv
rm fusion_statistics.tsv filtered_fusions_statistics.tsv

# 7. to do: cocordant statistics

# 8. delete temp directory
echo "Finished!"
cd .. && rm -rf ${output_directory}/temp_output_SFyNCS


endTime=$( date +%H:%M:%S)
echo -e "END:\t$endTime"
runTime=$( runningTime $startTime $endTime )
echo -e "RUNNING:\t$runTime"
