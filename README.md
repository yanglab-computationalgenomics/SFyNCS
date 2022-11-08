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
- 2.1. Reference genome sequence in fasta format
- 2.2. STAR index of the reference genome (can be skipped if Chimeric.out.junction is available). The STAR index can be built with the following command. We suggest using chr1-chr22, chrX, chrY and chrM only. Please specify the number_of_thread in the command below.
   ```
   mkdir -p /path_to/star_index_dir
   STAR --runThreadN number_of_thread --runMode genomeGenerate --genomeDir /path_to/star_index_dir --genomeFastaFiles /path_to/genome.fasta
   ```
- 2.3. TopHat/Bowtie2 index of the reference genome. The TopHat/Bowtie2 index can be built with the following command. Please specify file_prefix in the commands below.
   ```
   mkdir -p /path_to/tophat_index_dir
   bowtie2-build /path_to/genome.fasta /path_to/tophat_index_dir/file_prefix
   ln -s /path_to/genome.fasta /path_to/tophat_index_dir/file_prefix.fa
   ``` 
 - 2.4. Fai index of the reference genome (automatically generated by samtools)
 - 2.5. Gene annotation of the reference genome in [Gene Predictions Extended (gpe) format](https://genome.ucsc.edu/FAQ/FAQformat.html#format9). A header is needed (can be any artifact header). X chromosome should be chrX and Y chromosome should be chrY. Gencode_v29_hg38.gpe is provided with the package.        
 Example of gene annotation:
   ```
   Column 1:   Transcript ID (e.g., ENST00000485503.1)
   Column 2:   Chromosome (e.g., chr7)
   Column 3:   Strand (e.g., +)
   Column 4:   Transcript start position (e.g., 55192810)
   Column 5:   Transcript end position (e.g., 55200802)
   Column 6:   Coding region (CDS) start position (e.g., 55200802)
   Column 7:   Coding region end position (e.g., 55200802)
   Column 8:   Number of exons (e.g., 3)
   Column 9:   Exon start positions (e.g., 55192810,55198716,55200315)
   Column 10:  Exon end positions (e.g., 55192841,55198863,55200802)
   Column 11:  Score (e.g., 0)
   Column 12:  Gene symbol (e.g., EGFR)
   Column 13:  Status of CDS start annotation (none, unknown, incomplete, or complete. e.g., none)
   Column 14:  Status of CDS end annotation (none, unknown, incomplete, or complete. e.g., none)
   Column 15:  Exon frame offsets (e.g., -1,-1,-1)
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
       Column 1:  Chromosome of breakpoint 1 (e.g., chr1)
       Column 2:  Position of breakpoint 1 (e.g., 10005917)
       Column 3:  Strand of breakpoint 1 (e.g., +)
       Column 4:  Chromosome of breakpoint 2 (e.g., chr1)
       Column 5:  Position of breakpoint 2 (e.g., 10005936)
       Column 6:  Strand of breakpoint 2 (e.g., +)
       Column 7:  Breakpoint supported read count (e.g., 1)
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
  - 5.1. Start from fastq.gz files.
    ```
    bash ../run_SFyNCS.sh -p 1 -o demo_output -a toy_gene_annotation.gpe -g toy_reference_genome_sequence.fasta -s star_index -t tophat_index/tophat -d toy_normal_directory  toy_pair_end_reads_1.fastq.gz toy_pair_end_reads_2.fastq.gz
    ```  
  - 5.2.  Start from Chimeric.out.junction produced by STAR and fastq.gz files.
    ```
    bash ../run_SFyNCS.sh -p 1 -c Chimeric.out.junction -o demo_output -a toy_gene_annotation.gpe -g toy_reference_genome_sequence.fasta -t tophat_index/tophat -d toy_normal_directory toy_pair_end_reads_1.fastq.gz toy_pair_end_reads_2.fastq.gz
    ```

## IV. Workflow
### 1. Running SFyNCS  
SFyNCS can be run by the following commands. Users can further provide "-p thread_numbers" to speed up the steps of STAR and Tophat.
- 1.1. Start from fastq.gz files.
  ```
  /path/to/run_SFyNCS.sh -o /path_to/output_dir -a /path_to/gene_annotation.gpe -g /path_to/genome.fasta -s /path_to/star_index_dir -t /path_to/tophat_index_dir/file_prefix -d /path/to/normal_junction_directory 1.fastq.gz 2.fastq.gz
  ```
- 1.2. Start from Chimeric.out.junction produced by STAR and fastq.gz files.
  ```
  /path/to/run_SFyNCS.sh -c /path_to/Chimeric.out.junction -o /path_to/output_dir -a /path_to/gene_annotation.gpe -g /path_to/genome.fasta -t /path_to/tophat_index_dir/file_prefix -d /path/to/normal_junction_directory 1.fastq.gz 2.fastq.gz
  ```
 - 1.3. run_SFyNCS.sh options. Note that the default parameters are optimized for TCGA data.
   ```
   -a    --annotation_file                   STR      Gene annotation gpe file (section 2.5/II)  
   -g    --genome_fasta                      STR      Reference genome fasta file  
   -o    --output_directory                  STR      Output directory [default: current directory]  
   -p    --thread_number                     INT      Number of threads [default: 1]. Multiple threads can speed up the steps of STAR and TopHat  
   -s    --star_index                        STR      Path to STAR index. This option can be skipped if "-c" is provided  
   -t    --tophat_index                      STR      Path to TopHat index. It includes the name of any of index files up to but not including the first period  
   -c    --chimeric_file                     STR      Chimeric.out.junction file generated by STAR
   -d    --normal_junction_dir               STR      Directory contains normal samples' junctions (section 3.2/II)
         --adjust_adjacent_distance          INT      Breakpoints within this distance will be adjusted [default: 5]
         --cluster_distance                  INT      Split reads and read pairs within this distance will be clustered together [default: 1000000]
         --min_split_reads                   INT      Minimal number of split reads for a fusion transcript to be identified [default: 1]
         --min_read_pairs                    INT      Minimal number of read pairs for a fusion transcript to be identified [default: 1]
         --min_total_reads                   INT      Minimal number of total reads for a fusion transcript to be identified [default: 3]
         --overhang_length                   INT      Breakpoint overhang length for TopHat [default: 5]
         --read_pair_distance                INT      Maximal distance between the TopHat alignment of read pairs and breakpoint [default: 10000]
         --max_mate_split_read_distance      INT      Maximal distance between alignments if both read 1 and read 2 are split reads [default: 200]
         --motif_searching_length_in_blat    INT      Splice site motifs (GT in the donor, AAG/CAG/TAG in the acceptor) are searched within this window size of breakpoints [default: 5]
         --length_1_fusion_in_blat           INT      Specify the window of breakpoint flanking sequence for artificial reference [default: 1000000]
         --length_2_fusion_in_blat           INT      Specify the window of breakpoint flanking sequence for artificial reference [default: 100]
         --length_for_identity_in_blat       INT      Flanking sequences size of both fusion breakpoints to calculate sequence identity [default: 10]
         --align_percentage_in_blat          FLOAT    Minimal percentage of bases of the whole read that is alignable by Blat for split reads [default: 0.9]
         --max_split_read_blat_distance      INT      Maximal distance between Blat alignment and breakpoints [default: 5]
         --max_sequence_identity_in_blat     FLOAT    Maximal sequence identity between flanking sequences of two fusion breakpoints [default: 0.8]
         --filter_by_canonical_splice_motif  STR      Filter by canonical splice site motif [default: Y]
         --length_in_sd                      INT      Breakpoint flanking size to calculate standard deviation [default: 100]
         --sd_cutoff                         FLOAT    Standard deviation cutoff to filter fusions [default: 0.15]
         --filter_in_the_same_gene           STR      Filter fusions in the same gene [default: Y]
         --normal_adjacent_distance          INT      Window size to search for breakpoints in normal samples [default: 500]
         --normal_read_count_cutoff          INT      Filter fusions if numbers of reads (discordant pairs or split reads) in normal samples are equal to or more than the specified value [default: 2]
         --min_distance_plus_minus           INT      Minimal distance allowed between fusion breakpoints if located in the same chromosome and strand is plus (smaller coordinate) minus (larger coordinate) [default: 500000]
         --min_distance_non_plus_minus       INT      Minimal distance allowed between breakpoints if located in the same chromosome and strand is not plus minus [default: 100000]
   ```
 
 ### 2. Output
 There are two output files named fusions.tsv.gz and fusions_abridged.tsv.gz. The second one contains subset columns of the first file and is much smaller. Here are the columns of fusions.tsv.gz.
 ```
Column 1:   Chromosome of breakpoint 1 (e.g., chr1)
Column 2:   Position of breakpoint 1 (e.g., 46389)
Column 3:   Strand of breakpoint 1 (e.g., +)
Column 4:   Chromosome of breakpoint 2 (e.g., chr2)
Column 5:   Position of breakpoint 2 (e.g., 16337)
Column 6:   Strand of breakpoint 2 (e.g., -)
Column 7:   Split read count reported by STAR (e.g., 1)
Column 8:   Read pair count reported by STAR (e.g., 4)
Column 9:   Split read count supported by TopHat (e.g., 1)
Column 10:  Potential split read count not supported by TopHat (e.g., 0)
Column 11:  Read pair count reported by TopHat (e.g., 3)
Column 12:  Split read count supported by Blat (e.g., 1)
Column 13:  Split read count supported by Blat (considering split read supported by TopHat only, e.g., 1)
Column 14:  Minimal distance between read pair and breakpoint 1 after aligning by TopHat (e.g., 13)
Column 15:  Minimal distance between read pair and breakpoint 2 after aligning by TopHat (e.g., 22)
Column 16:  Sequence identity (e.g., 0.52)
Column 17:  Minimal distance between split read and breakpoint 1 after aligning by Blat (e.g., 0)
Column 18:  Minimal distance between split read and breakpoint 2 after aligning by Blat (e.g., 0)
Column 19:  Minimal distance between split read and breakpoint 1 after aligning by Blat (considering split read supported by TopHat only, e.g., 0)
Column 20:  Minimal distance between split read and breakpoint 2 after aligning by Blat (considering split read supported by TopHat only, e.g., 0)
Column 21:  Presence of canonical splice site (Y or N)
Column 22:  Cluster count around breakpoint 1 (e.g., 1)
Column 23:  Cluster count around breakpoint 2 (e.g., 1)
Column 24:  Cluster count around breakpoint 1 and breakpoint 2 (e.g., 1)
Column 25:  Percentage of split reads and read pairs supported fusion transcript’s cluster around breakpoint 1 and breakpoint 2 (e.g., 0.75)
Column 26:  Standard deviation for candidate fusion clusters around breakpoint 1 (e.g., NA, NA for only one cluster)
Column 27:  Standard deviation for candidate fusion clusters around breakpoint 2 (e.g., NA, NA for only one cluster)
Column 28:  Whether fusion transcript locates in the same gene (Y or N)
Column 29:  Fusion type (non-coding or protein-coding)
Column 30:  Overlapped gene in breakpoint 1 (e.g., C9/DAB2, genes were separated by “/”)
Column 31:  Gene type in breakpoint 1 (non-coding_gene, protein-coding_gene or unknown, different genes’ type were separated by “/”)
Column 32:  Gene strand in breakpoint 1 (e.g., +)
Column 33:  Breakpoint 1’s location (intron, exon, splice_site or unknown)
Column 34:  Breakpoint 1’s region type (5’UTR, 3’UTR, CDS, non-coding or unknown)
Column 35:  Exon frame after breakpoint 1 (0, 1, 2 or N, N means unknown)
Column 36:  Overlapped gene in breakpoint 2 (e.g., EML4)
Column 37:  Gene type in breakpoint 2 (non-coding_gene, protein-coding_gene or unknown, different genes’ type were separated by “/”)
Column 38:  Gene strand in breakpoint 2 (e.g., -)
Column 39:  Breakpoint 2’s location (intron, exon, splice_site or unknown)
Column 40:  Breakpoint 2’s region type (5’UTR, 3’UTR, CDS, non-coding or unknown)
Column 41:  Exon frame before breakpoint 2 (0, 1, 2 or N, N means unknown)
Column 42:  Fusion frame (in-frame, out-frame, or unknown)
Column 43:  Split reads reported by STAR (e.g., read_23)
Column 44:  Read pairs reported by STAR (e.g., read_38,read_62,read_70,read_8)
Column 45:  Split reads supported by TopHat (e.g., read_23)
Column 46:  Potential split reads not supported by TopHat (e.g., NA)
Column 47:  Read pairs supported by TopHat (e.g., read_38,read_70,read_8)
Column 48:  Split reads supported by Blat (e.g., read_23)
Column 49:  Distance between each read pair and breakpoint 1 after aligning by TopHat (e.g., 740,23,13)
Column 50:  Distance between each read pair and breakpoint 2 after aligning by TopHat (e.g., 292,2592,22)
Column 51:  Distance between each split read and breakpoint 1 after aligning by Blat (e.g., 0)
Column 52:  Distance between each split read and breakpoint 2 after aligning by Blat (e.g., 0)
Column 53:  Distance between each split read and breakpoint 1 after aligning by Blat (considering split read supported by TopHat only, e.g., 0)
Column 54:  Distance between each split read and breakpoint 2 after aligning by Blat (considering split read supported by TopHat only, e.g., 0)
Column 55:  Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 1 (e.g., --AGTGGGCCAGGTAG-GGCTGG)
Column 56:  Needleman-Wunsch alignment of sequence used to calculate sequence identity in breakpoint 2 (e.g., CCACT--GCCAGG-AGAACCTCA)
Column 57:  Discordant read pair cluster IDs around breakpoint 1 (e.g., 1)
Column 58:  Discordant read pair cluster IDs around breakpoint 2 (e.g., 1)
Column 59:  Supporting reads count in each cluster around breakpoint 1 (e.g., 4)
Column 60:  Supporting reads count in each cluster around breakpoint 2 (e.g., 3)
 ```
## V. Reference
## VI. Contact
If users have questions, please contact Xiaoming Zhong (xiaomingzhong@uchicago.edu) and Lixing Yang (lixingyang@uchicago.edu).
## VII. Issue
- use tophat 2.1.0 if having following issue:
```
"./SeqAn-1.4.2/seqan/basic/basic_exception.h:236 FAILED!  (Uncaught exception of type St12out_of_range: basic_string::substr: __pos (which is 40) > this->size() (which is 0))"
```
