#!/bin/bash

# 2. before blat
discordant_count=$(grep -v "^#" Chimeric.out.junction | wc -l)
if [ $discordant_count -eq 0 ]; then
  echo "No discordant reads" >note_no_discordant.tsv
  exit
fi

perl $toolDir/format_STAR_chimeric_file.pl  Chimeric.out.junction >format_chimeric.tsv
perl $toolDir/remove_multiple_mapped_reads.pl format_chimeric.tsv >no_multiple-mapped.tsv
perl $toolDir/remove_duplicate_reads.pl no_multiple-mapped.tsv | sort | uniq >no_duplicate.tsv
perl $toolDir/merge_ajacent_breakpoints.pl no_duplicate.tsv >merge_adjacent.tsv
perl $toolDir/cluster_discordant_reads.pl merge_adjacent.tsv >cluster.tsv
perl $toolDir/identify_fusion_candidates_from_cluster_reads.pl cluster.tsv | awk '$9>0' | grep -v chrM >potential_candidates.tsv
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2-1,$2,$3,$7,$8,$9; print $4,$5-1,$5,$6,$7,$8,$9;}' cluster.tsv | sort -k1,1 -k2,2n | uniq >temp_cluster.bed
rm Chimeric.out.junction cluster.tsv format_chimeric.tsv merge_adjacent.tsv no_duplicate.tsv no_multiple-mapped.tsv

cut -f11,12 potential_candidates.tsv | sed "s#\t#\n#;s#,#\n#g" | grep -vw "NA" | grep -vP "Split_reads|Read_pairs" | sort | uniq >potential_candidates_reads.tsv
perl $toolDir/select_fq.pl potential_candidates_reads.tsv discordent-1.fq discordent-2.fq >potential_candidates_reads_1.fq 2>potential_candidates_reads_2.fq
rm discordent-{1,2}.fq  potential_candidates_reads.tsv

tophat -p 3 --no-coverage-search --fusion-search --fusion-anchor-length 12 --fusion-min-dist 100000 -o tophat_out --read-mismatches 4 --read-gap-length 4 --read-edit-dist 4 --splice-mismatches 2 --max-insertion-length 4 --max-deletion-length4 --segment-mismatches 3 --fusion-read-mismatches 4 /scratch/xiaomingzhong/fusion-trans/tools/tophatIndex/ref potential_candidates_reads_*fq >log/tophat.log 2>log/tophat.err
ln -s tophat_out/accepted_hits.bam tophat.bam
samtools index tophat.bam
perl $toolDir/processe_read_supported_by_tophat.pl tophat.bam  potential_candidates.tsv | awk 'NR==1 || $8>0' >processed_with_tophat.tsv
rm -rf tophat*

# 3. split tophat processed
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

# 4. blat
for i in {1..50}
  do
    cd chunk_${i}
    ln -s ../potential_candidates_reads_1.fq .
    ln -s ../potential_candidates_reads_2.fq .
    perl $toolDir/processe_read_supported_by_blat.pl -f /scratch/xiaomingzhong/fusion-trans/tools/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa processed_with_tophat.tsv >processed_with_blat.tsv
    cd ..
  done


# 5. merge split blat result
cp chunk_1/processed_with_blat.tsv .
for((i=2; i<51; i++))
  do
    sed -n '2,$p' chunk_${i}/processed_with_blat.tsv >>processed_with_blat.tsv
  done
rm -rf chunk*
rm potential_candidates*

# 6. split blat processed
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

# 6. generate flanking statistics
for i in {1..50}
  do
    cd chunk_${i}
    ln -s ../temp_cluster.bed .
    perl $toolDir/cluster_statistics.pl processed_with_blat.tsv >fusion_statistics.tsv
    cd ..
  done

# 7. merge split flanking statistics
cp chunk_1/fusion_statistics.tsv .
for((i=2; i<51; i++))
  do
    sed -n '2,$p' chunk_${i}/fusion_statistics.tsv >>fusion_statistics.tsv
  done
rm -rf chunk*
