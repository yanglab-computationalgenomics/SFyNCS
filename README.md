## I. About
Somatic Fusions involving Non-Coding Sequences (SFyNCS) detects both coding and non-coding fusion transcripts from paired-end RNA-seq data. It uses discordant read pairs and split reads to detect the fusions and refine fusion breakpoints to basepair resolution.

## II. Prerequisites
### 1. Software environment
- Unix/Linux system
- Blat (tested v35, v36), can be downloaded at http://hgdownload.soe.ucsc.edu/admin/exe or installed with conda (https://anaconda.org/bioconda/blat).
- Bedtools (tested v2.25.0, v2.26.0, v2.27.1), can be downloaded at https://github.com/arq5x/bedtools2/releases or installed with conda (https://anaconda.org/bioconda/bedtools).
- Bowtie2 (tested v2.1.0, v2.2.5, v2.2.9, v2.3.0, v2.3.2, v2.3.4.3), can be downloaded at https://sourceforge.net/projects/bowtie-bio/files/bowtie2/ or installed with conda (https://anaconda.org/bioconda/bowtie2).
- Perl (v5.010 or above), can be downloaded at https://www.perl.org/get.html or installed with conda (https://anaconda.org/conda-forge/perl).
- STAR (v2.6.1a or above, tested v2.6.1a, v2.6.1d, v2.7.0f), can be downloaded at https://github.com/alexdobin/STAR/releases or installed with conda (https://anaconda.org/bioconda/star).
- Samtools (tested v0.1.19, v1.1, v1.3.1, v1.5), can be downloaded at https://github.com/samtools/samtools/releases or installed with conda (https://anaconda.org/bioconda/samtools).
- TopHat (must be v2.1.0), can be downloaded at http://ccb.jhu.edu/software/tophat/downloads or installed with conda (https://anaconda.org/bioconda/tophat, please specify version when installing).  
Note that TopHat needs to be v2.1.0 (2.1.1 or 2.0.13 won't work).            
Follow these steps to install:
   ```
   wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz
   tar -zxvf tophat-2.1.0.Linux_x86_64.tar.gz
   export PATH=$PWD/tophat-2.1.0.Linux_x86_64:$PATH
   chmod u+x tophat-2.1.0.Linux_x86_64/*
   ```
### 2. Reference, index and annotation files
- Reference genome sequence in fasta format
- STAR index of the reference genome (can be skipped if Chimeric.out.junction is available). The STAR index can be built with the following command. We suggest using chr1-chr22, chrX, chrY and chrM only. Please specify the number_of_thread in the command below.
   ```
   mkdir -p /path_to/star_index_dir
   STAR --runThreadN number_of_thread --runMode genomeGenerate --genomeDir /path_to/star_index_dir --genomeFastaFiles /path_to/genome.fasta
   ```
- TopHat/Bowtie2 index of the reference genome. The TopHat/Bowtie2 index can be built with the following command. Please specify file_prefix in the commands below.
   ```
   mkdir -p /path_to/tophat_index_dir
   bowtie2-build /path_to/genome.fasta /path_to/tophat_index_dir/file_prefix
   ln -s /path_to/genome.fasta /path_to/tophat_index_dir/file_prefix.fa
   ``` 
 - Fai index of the reference genome (automatically generated by samtools)
 - Gene annotation of the reference genome in [Gene Predictions Extended (gpe) format](https://genome.ucsc.edu/FAQ/FAQformat.html#format9). A header is needed (can be any artifact header). X chromosome should be chrX and Y chromosome should be chrY. Gencode_v29_hg38.gpe is provided with the package.        
 Example of gene annotation:
   ```
   Column 1: Transcript ID (e.g., ENST00000485503.1)
   Column 2: Chromosome (e.g., chr7)
   Column 3: Strand (e.g., +)
   Column 4: Transcript start position (e.g., 55192810)
   Column 5: Transcript end position (e.g., 55200802)
   Column 6: Coding region (CDS) start position (e.g., 55200802)
   Column 7: Coding region end position (e.g., 55200802)
   Column 8: Number of exons (e.g., 3)
   Column 9: Exon start positions (e.g., 55192810,55198716,55200315)
   Column 10: Exon end positions (e.g., 55192841,55198863,55200802)
   Column 11: Score (e.g., 0)
   Column 12: Gene symbol (e.g., EGFR)
   Column 13: Status of CDS start annotation (none, unknown, incomplete, or complete. e.g., none)
   Column 14: Status of CDS end annotation (none, unknown, incomplete, or complete. e.g., none)
   Column 15: Exon frame offsets (e.g., -1,-1,-1)
   ```
### 3. Input files
  - 3.1. Two paired-end fastq.gz files or fastq files. SFyNCS uses STAR to align reads and use Chimeric.out.junction generated by STAR to call fusion candidates. If users have Chimeric.out.junction prior to running SFyNCS, the STAR alignment step can be skipped, and it can speed up the process significantly.    
    Chimeric.out.junction file can be produced by STAR with following command, please specify /path_to/star_index_dir, number_of_thread and /path_to/star_output_dir in the command:
    ```
    STAR   --genomeDir /path_to/star_index_dir \
           --readFilesIn 1.fastq.gz 2.fastq.gz \
           --readFilesCommand zcat \
           --runThreadN number_of_thread \
           --outFileNamePrefix /path_to/star_output_dir \
           --outReadsUnmapped None \
           --twopassMode Basic \
           --outSAMstrandField intronMotif \
           --outSAMunmapped Within \
           --chimSegmentMin 12 \
           --chimJunctionOverhangMin 12 \
           --chimOutJunctionFormat 1 \
           --alignSJDBoverhangMin 10 \
           --alignMatesGapMax 100000 \
           --alignIntronMax 100000 \
           --alignSJstitchMismatchNmax 5 -1 5 5 \
           --outSAMattrRGline ID:GRPundef  \
           --chimMultimapScoreRange 10 \
           --chimMultimapNmax 10 \
           --chimNonchimScoreDropMin 10 \
           --peOverlapNbasesMin 12 \
           --peOverlapMMp 0.1 \
           --outSAMtype BAM SortedByCoordinate \
           --genomeLoad NoSharedMemory
    ```
  - 3.2. Fusion breakpoints from normal samples. They are used to filter out germline events as well as artifacts. We provide fusion breakpoint files from TCGA normal samples. Users can generate their own if RNAseq data from normal samples are available:
     - Generate Chimeric.out.junction for each normal sample (refer to 3.1/II).
     - Format split reads and discordant read pairs information in Chimeric.out.junction.
       ```
       perl /path/to/format_STAR_chimeric_file.pl Chimeric.out.junction >format_chimeric.tsv
       ```
     - Remove duplicated reads and get the number of supporting reads for each breakpoint.
       ```
       perl /path/to/get_juntions_in_normal_sample.pl format_chimeric.tsv >no_duplication_junction_read_count.tsv
       ```
       Example of gene annotation:
       ```
       Column 1: Chromosome of breakpoint 1 (e.g., chr1)
       Column 2: Position of breakpoint 1 (e.g., 10005917)
       Column 3: Strand of breakpoint 1 (e.g., +)
       Column 4: Chromosome of breakpoint 2 (e.g., chr1)
       Column 5: Position of breakpoint 2 (e.g., 10005936)
       Column 6: Strand of breakpoint 2 (e.g., +)
       Column 7: Breakpoint supported read count (e.g., 1)
       ```
     - Rename each sample’s no_duplication_junction_read_count.tsv and put them under the same directory. Replace “sample_1” with each sample id.
       ```
       mv no_duplication_junction_read_count.tsv sample_1.tsv
       ```
     - Optional: each normal junction file can be compressed.
       ```
       gzip sample_1.tsv
       ```

## III. Example
Note: all files under the example directory are for testing only.
1. Download and decompress SFyNCS.
```
wget https://github.com/yanglab-computationalgenomics/SFyNCS/
mv master.zip SFyNCS.zip
unzip SFyNCS.zip
```
2. Decompress example.tar.gz and enter example directory.
```
cd SFyNCS
tar -zxvf example.tar.gz
cd example
```
3. Build STAR index.
```
mkdir star_index
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles toy_reference_genome_sequence.fasta
```
4. Build TopHat/Bowtie2 index.
```
mkdir -p tophat_index
bowtie2-build toy_reference_genome_sequence.fasta tophat_index/tophat
ln -s $PWD/toy_reference_genome_sequence.fasta $PWD/tophat_index/tophat.fa
```
5. Run SFyNCS. Users should get fusions.tsv.gz and fusions_abridged.tsv.gz under demo_output, please make sure the content in these two output files is the same as the files under example.  
  5.1. Start from fastq.gz files.
  ```
  bash ../run_SFyNCS.sh -p 1 -o demo_output -a toy_gene_annotation.gpe -g toy_reference_genome_sequence.fasta -s star_index -t tophat_index/tophat -d toy_normal_directory  toy_pair_end_reads_1.fastq.gz toy_pair_end_reads_2.fastq.gz
  ```  
   5.2.  Start from Chimeric.out.junction produced by STAR and fastq.gz files.
  ```
  bash ../run_SFyNCS.sh -p 1 -c Chimeric.out.junction -o demo_output -a toy_gene_annotation.gpe -g toy_reference_genome_sequence.fasta -t tophat_index/tophat -d toy_normal_directory toy_pair_end_reads_1.fastq.gz toy_pair_end_reads_2.fastq.gz
  ```
  

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
mkdir -p /path/to/directory/to/store/STAR/index/files
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir /path/to/directory/to/store/STAR/index/files --genomeFastaFiles /path/to/genome/fasta
```
2. Tophat index can be builded by following command.
```
bowtie2-build /path/to/genome/fasta /path/to/directory/to/store/Tophat/index/files/file_prefix
# tophat need fasta sequence, make a soft link
ln -s /path/to/genome/fasta /path/to/directory/to/store/Tophat/index/files/file_prefix.fasta
```
3. Gene annotation is in [Gene Predictions (Extended)](https://genome.ucsc.edu/FAQ/FAQformat.html#format9) format and header is needed (can be any artifact header). X chromosome should be chrX and Y chromosome should be chrY.
```
Transcript_id           Chr     Strand   Transcript_start_position  Transcript_end_position CDS_start_position  CDS_end_position  Exon_count    Exon_starts                     Exon_ends                     Score   Symbol  CDS_start_stat  CDS_end_stat  Exon_frame
ENST00000485503.1       chr7    +         55192810                  55200802                55200802            55200802          3             55192810,55198716,55200315,     55192841,55198863,55200802,   0       EGFR    none            none          -1,-1,-1,
```
4. SFyNCS can be run by following command. User can further provide "-p thread_numbers" to speed the step of STAR and Tophat.
```
/path/to/run_SFyNCS.sh -o /path/to/output/direcotry -a /path/to/gene_annotation_file -g /path/to/fasta -s /path/to/directory/to/store/STAR/index/files -t /path/to/directory/to/store/Tophat/index/files 1.fastq 2.fastq
```

#### Example
```
# 1. Build STAR index as above
# 2. Build Tophat index as above
# 3. tar -zxvf example_data.tar.gz
# 4. bash SFyNCS_sherlock/run_SFyNCS.sh -p 1 -o demo_output -a gencode_v29.gpe -g reference/ref_genome.fa -s reference/star -t reference/tophat/tophat example_1.fastq example_2.fastq
# 5. You should get fusions.tsv.gz and fusions_abridged.tsv.gz under demo_output, please make sure these two output files are same as files under example_output 
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

For gene with multiple transcripts, the longest transcript defined as total length of exons was used in generating column **Fusion_annotations**. If there are more than two genes overlapped with fusion, they are seperated by ";". Each paired gene are seperated by "," and in following format:
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

#### Issue
use tophat 2.1.0 if having following issue:
"./SeqAn-1.4.2/seqan/basic/basic_exception.h:236 FAILED!  (Uncaught exception of type St12out_of_range: basic_string::substr: __pos (which is 40) > this->size() (which is 0))"
