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
3. Gene annotation contain 6 columns and must have header. X chromosome should be chrX and y chromosome should be chrY.
```
Chr     Start       End         Strand      Symbol      Gene_type
chr7    55019021    55211628    +           EGFR        protein_coding
```
4. SFyNCS can be run as:
```
/path/to/run_SFyNCS.sh -o /path/to/output/direcotry -a /path/to/gene_annotation_file -g /path/to/fasta -s /path/to/directory/to/store/STAR/index/files -t /path/to/directory/to/store/Tophat/index/files
```

#### Output
The output is a tab-delimited file named "fusions.tsv" with the following format:
```
```

