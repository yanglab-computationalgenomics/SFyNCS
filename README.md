## About
SFyNCS detect fusion transcripts from pair-end RNA-seq data

## Prerequisites
- Unix/Linux system
- Blat (testing with v35)
- Bedtools (testing with v2.27.1)
- Bowtie2 (testing with v2.2.5)
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
3. Gene annotation is in [Gene Predictions (Extended)](https://genome.ucsc.edu/FAQ/FAQformat.html#format9) format and header is needed (can be any artifact header). X chromosome should be chrX and Y chromosome should be chrY.
```
Transcript_id           Chr     Strand   Transcript_start_position  Transcript_end_position CDS_start_position  CDS_end_position  Exon_count    Exon_starts                     Exon_ends                     Score   Symbol  CDS_start_stat  CDS_end_stat  Exon_frame
ENST00000485503.1       chr7    +         55192810                  55200802                55200802            55200802          3             55192810,55198716,55200315,     55192841,55198863,55200802,   0       EGFR    none            none          -1,-1,-1,
```
4. SFyNCS can be run by following command. User can further provide "-p thread_numbers" to speed the step of STAR and Tophat.
```
/path/to/run_SFyNCS.sh -o /path/to/output/direcotry -a /path/to/gene_annotation_file -g /path/to/fasta -s /path/to/directory/to/store/STAR/index/files -t /path/to/directory/to/store/Tophat/index/files
```

#### Output
The output is a tab-delimited file named "fusions.tsv" with the following format:
| Chr_1 | Breakpoint_1 | Strand_1 | Chr_2 | Breakpoint_2 | Strand_2 | Split_read_count_(Tophat_and_Blat) | Read_pair_count_(Tophat) | Minimum_read_distance_to_left | Minimum_read_distance_to_right | Identity | Minimum_blat_distance_to_1 | Minimum_blat_distance_to_2 | Total_clusters_1 | Total_clusters_2 | Total_clusters_(merge) | (discordant_reads)%\_support_fusion | SD_(discordant_reads)%\_in_clusters_1 | SD_(discordant_reads)%\_in_clusters_2 | Fusion_annotations | Split_reads_(tophat_and_blat) | Read_pairs_(tophat) |
| :-- | :-- | :-- | :-- | :-- | :--  | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| chr17 | 7314362 | - | chr9 | 2876236 | +  | 1 | 2 | 99 | 0 | 0.59 | 0 | 0 | 2 | 6 | 7 | 0.05 | 0.36 | 0.29 | AC026954.2--GPS2P1,ENST00000575474.1--ENST00000411761.2,protein_coding_gene--non_coding_gene,exon--exon,3'UTR--non_coding;GPS2--GPS2P1,ENST00000571697.5--ENST00000411761.2,non_coding_gene--non_coding_gene,exon--exon,non_coding--non_coding | UNC11-SN627:278:D1LY0ACXX:7:1308:7484:19502 | UNC11-SN627:278:D1LY0ACXX:7:1205:12470:97470,UNC11-SN627:278:D1LY0ACXX:7:1305:18043:40271 |

1. **Chr_1:** chromosome of breakpoint 1  
2. **Breakpoint_1:** breakpoint 1 coordinates <br>
3. **Strand_1:** strand of breakpoint 1 <br>
4. **Chr_2:** chromosome of breakpoint 2  
5. **Breakpoint_2:** breakpoint 2 coordinates <br>
6. **Strand_2:** strand of breakpoint 2 <br>
7. **Split_read_count_(Tophat_and_Blat):** split read count (processed by tophat and blat) <br>
8. **Read_pair_count_(Tophat):** read pair count (processed by tophat) <br>
9. **Minimum_read_distance_to_1:** minimum distance of read pair to breakpoint 1 <br>
10. **Minimum_read_distance_to_2:** minimum distance of read pair to breakpoint 2 <br>
11. **Identity:** identity <br>
12. **Minimum_blat_distance_to_1:** minimum blat distace to breakpoint 1 (use mean if both read 1 and read 2 are split read) <br>
13. **Minimum_blat_distance_to_2:** minimum blat distace to breakpoint 2 (use mean if both read 1 and read 2 are split read) <br>
14. **Total_clusters_1:** total clusters in flanking region of breakpoint 1 <br>
15. **Total_clusters_2:** total clusters in flanking region of breakpoint 2 <br>
16. **Total_clusters_(merge):** total clusters in flanking region of breakpoint 1 and breakpint 2 <br>
17. **(discordant_reads)%\_support_fusion:** (discordant read)% support fusion <br>
18. **SD_(discordant_reads)%\_in_clusters_1:** sd of (discordant read)% in each of cluster in flanking region of breakpoint 1 <br>
19. **SD_(discordant_reads)%\_in_clusters_2:** sd of (discordant read)% in each of cluster in flanking region of breakpoint 2 <br>
20. **Fusion_annotations:** fusion annotations <br>
21. **Split_reads_(tophat_and_blat):** split reads (processed by tophat and blat) <br>
22. **Read_pairs_(tophat):** read pairs (processed by tophat) <br>

For gene with multiple transcripts, the longest transcript defined as total length of exons was used in generating column **Fusion_annotations**. If there are more than two genes overlapped with fusion, they are seperated by ";". Each pair gene are seperated by "," and in following format:
- symbol--symbol (e.g. AC026954.2--GPS2P1, will be "unannotated" if no annotated genes)
- transcript_id--transcript_id (e.g. ENST00000575474.1--ENST00000411761.2, will be "unannotated" if no annotated genes)
- gene_type--gene_type (e.g. protein_coding_gene--non_coding_gene, will be "unknown" if no annotated genes)
- breakpoint_position--breakpoint_position (e.g. exon--exon, it can be one of "intergenic", "intron", "exon", "split_site")
- fusion_region_type--fusion_region_type (e.g. 3'UTR--non_coding, it can be one of "intergenic", "non_coding", "5'UTR", "3'UTR", "cds"). If the type is "cds--cds", "in-frame" or "frame-shift" will be added (e.g. cds--cds_(in-frame)) 
- fusion position (e.g. chr2:200:---chr1:110:+, this information will be provided in scenario below only):
```
           chr1:110+         chr2:200:-
    -----------|                 |-----------
       gene_A --->             ---> gene_B    (---> is genes' orientation)
       gene_C <---             <--- gene_D    (<--- is genes' orientation)

There are two potential fusions: gene_A--gene_B fusion and gene_D--gene_C fusion, fusion position will added to gene_D--gene_C to keep one row only in output and it would be like: gene_D--gene_C,transcript_id_D--transcript_id_C,gene_type_D--gene_type_C,breakpoint_position_D--breakpoint_position_C,fusion_region_type_D--fusion_region_type_C,chr2:200:---chr1:110:+
```

