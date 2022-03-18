## About
SFyNCS detect fusion transcripts from pair-end RNA-seq data

## Prerequisites
- Unix/Linux system
- Blat (testing with v35)
- Bedtools (testing with v2.27.1)
- Perl (5.010 or above)
- STAR (testing with 2.7.0f)
- Tophat (testing with v2.1.0)

## How to use SFyNCS
#### Installation
```
wget https://github.com/yanglab-computationalgenomics/SFyNCS/archive/refs/heads/main.zip
tar -zxvf main.zip
```
or
```
git clone https://github.com/yanglab-computationalgenomics/SFyNCS.git
```

#### Workflow
SFyNCS requir fastq, STAR index, Tophat index, genome fasta and gene annotation to identified fusion.
1. STAR index can be builded by following command and please refer to [STAR](https://github.com/alexdobin/STAR) for more detail. Although SFyNCS can run with unlocalized/unplaced sequences, we sugget use chr1-chr22, chrX, chrY and chrM fasta.
```
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir /path/to/directory/to/store/STAR/index/files --genomeFastaFiles /path/to/genome/fasta
```
2. Tophat index can be builded by following command.
```
bowtie-build /path/to/genome/fasta /path/to/directory/to/store/Tophat/index/files
```
3. Gene annotation contain 6 columns and must have header. Start and End are 1-base. X chromosome should be chrX and Y chromosome should be chrY.
```
Chr     Start       End         Strand      Symbol      Gene_type
chr7    55019021    55211628    +           EGFR        protein_coding
```
4. SFyNCS can be run by following command. User can further provide "-p thread_numbers" to speed the step of STAR and Tophat.
```
/path/to/run_SFyNCS.sh -o /path/to/output/direcotry -a /path/to/gene_annotation_file -g /path/to/fasta -s /path/to/directory/to/store/STAR/index/files -t /path/to/directory/to/store/Tophat/index/files
```

#### Output
The output is a tab-delimited file named "fusions.tsv" with the following format:
| Chr_left | Pos_left | Strand_left | Chr_right | Pos_right | Strand_right | Split_read_count_(Tophat_and_Blat) | Read_pair_count_(Tophat) | Minimum_read_distance_to_left | Minimum_read_distance_to_right | Identity | Minimum_blat_distance_to_left | Minimum_blat_distance_to_right | Total_clusters_(left) | Total_clusters_(right) | Total_clusters_(merge) | (discordant_reads)%_support_fusion | SD_(discordant_reads)%_in_clusters_(left) | SD_(discordant_reads)%_in_clusters_(right) | Fusion_annotations | Split_reads_(tophat_and_blat) | Read_pairs_(tophat) |
| :-- | :-- | :-- | :-- | :-- | :--  | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| chr9 | 2876236 | + | chr17 | 7314362 | -  | 1 | 2 | 99 | 0 | 0.59 | 0 | 0 | 2 | 6 | 7 | 0.05 | 0.36 | 0.29 | AC026954.2_(protein_coding)--GPS2P1_(processed_pseudogene);GPS2_(protein_coding)--GPS2P1_(processed_pseudogene) | UNC11-SN627:278:D1LY0ACXX:7:1308:7484:19502 | UNC11-SN627:278:D1LY0ACXX:7:1205:12470:97470,UNC11-SN627:278:D1LY0ACXX:7:1305:18043:40271 |

**Chr_left:** chromosome of the left segment <br>
**Pos_left:** left segment site <br>
**Strand_left:** strand of the left segment <br>
**Chr_right:** chromosome of right segment <br>
**Pos_right:** right segment site <br>
**Strand_right:** strand of the right segment <br>
**Split_read_count_(Tophat_and_Blat):** split read count (processed by tophat and blat) <br>
**Read_pair_count_(Tophat):** read pair count (processed by tophat) <br>
**Minimum_read_distance_to_left:** minimum distance of read pair to left breakpoint <br>
**Minimum_read_distance_to_right:** minimum distance of read pair to right breakpoint <br>
**Identity:** identity <br>
**Minimum_blat_distance_to_left:** minimum blat distace of left breakpoint (use mean if both read 1 and read 2 are split read) <br>
**Minimum_blat_distance_to_right:** minimum blat distace of right breakpoint (use mean if both read 1 and read 2 are split read) <br>
**Total_clusters_(left):** total clusters in cluster (left) <br>
**Total_clusters_(right):** total clusters in cluster (right) <br>
**Total_clusters_(merge):** total clusters in cluster (merge) <br>
**(discordant_reads)%\_support_fusion:** (discordant read)% support fusion <br>
**SD_(discordant_reads)%\_in_clusters_(left):** sd of (discordant read)% in each of cluster (left) <br>
**SD_(discordant_reads)%\_in_clusters_(right):** sd of (discordant read)% in each of cluster (right) <br>
**Fusion_annotations:** fusion annotations <br>
**Split_reads_(tophat_and_blat):** split reads (processed by tophat and blat) <br>
**Read_pairs_(tophat):** read pairs (processed by tophat) <br>

