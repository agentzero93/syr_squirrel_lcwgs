# Genomic differentiation between urban and rural eastern gray squirrels in Syracuse, NY #
Repository to house scripts used to analyze low-coverage whole genome sequencing of the eastern gray squirrel (Sciurus carolinensis) in Syracuse, New York.

# Tools used (running list) #
### [Read data QC](https://github.com/agentzero93/syr_squirrel_lcwgs/edit/main/README.md#read-data-qc-only-running-syracuse-samples-as-of-now) ###
1) [FastQC v.0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2) [MultiQC v.1.12](https://multiqc.info/)
3) [Trimmomatic v.0.39](https://github.com/usadellab/Trimmomatic)
### [Alignment](https://github.com/agentzero93/syr_squirrel_lcwgs/edit/main/README.md#alignment-1) ###
1) [BWA v.0.7.17](https://github.com/lh3/bwa)
2) [GATK v.3.8.1](https://gatk.broadinstitute.org/hc/en-us)
3) [Picard v.2.25.6](https://github.com/broadinstitute/picard)
4) [seqkit v.2.3.0](https://github.com/shenwei356/seqkit)
5) [samtools v.1.13](http://www.htslib.org/)
6) [QualiMap v.2.2.2-dev](http://qualimap.conesalab.org/)
7) [MultiQC v.1.12](https://multiqc.info/)
### [SNP generation](https://github.com/agentzero93/syr_squirrel_lcwgs/edit/main/README.md#snp-generation-1) ###
1) [bcftools v.1.16](https://samtools.github.io/bcftools/bcftools.html)
2) [samtools v.1.13](http://www.htslib.org/)
### [Analyses](https://github.com/agentzero93/syr_squirrel_lcwgs/edit/main/README.md#analyses) ###
1)
2)
3)

# Read data QC (only running Syracuse samples as of now) #

Note: 37 samples were originally run at 3x, but chromosome coverage was only at 2x. So, we resequenced these samples at 7x to get to 5x chromosome coverage.

Note: JPV022 doesn't appear to be from Syracuse, will not run in pipeline (46 samples remaining).

1) Check the quality of the raw reads (.fastq.gz) using FastQC and MultiQC.
```bash
fastqc SCCA1009_raw_R1.fastq.gz --outdir ./ --threads 20
fastqc SCCA1009_raw_R2.fastq.gz --outdir ./ --threads 20
...
multiqc ./ --interactive
```
Note: SCCA1009 SCCA1011 SCCA1012 SCCA1017 SCCA1018 SCCA1020 SCCA1026 SCCA1033 SCCA1034 SCCA1040 were run at 10x, the rest were run at 7x.

Note: Raw reads look good, just a bit of adapter contamination.

2) Trim the raw reads of adapter contamination and low-quality bases using Trimmomatic (then reassess quality with FastQC and MultiQC).
```bash
java -jar trimmomatic-0.39.jar PE \
  -threads 20 \
  SCCA1009_raw_R1.fastq.gz SCCA1009_raw_R2.fastq.gz \
  SCCA1009_trimmed_1P.fastq.gz SCCA1009_trimmed_1U.fastq.gz SCCA1009_trimmed_2P.fastq.gz SCCA1009_trimmed_2U.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
  SLIDINGWINDOW:4:25
  
fastqc SCCA1009_trimmed_1P.fastq.gz --outdir ./ --threads 20
fastqc SCCA1009_trimmed_1U.fastq.gz --outdir ./ --threads 20
fastqc SCCA1009_trimmed_2P.fastq.gz --outdir ./ --threads 20
fastqc SCCA1009_trimmed_2P.fastq.gz --outdir ./ --threads 20
...
multiqc ./ --interactive
```
Note: Reads no longer have adapter contamination.

# Alignment #

1) Align the trimmed reads to the reference genome using BWA. e.g., paired-end alignment:
```bash
bwa mem -t 20 egsq_1MBmin SCCA1009_trimmed_1P.fastq.gz SCCA1009_trimmed_2P.fastq.gz | \
  samtools sort --threads 20 -o SCCA1009_paired.bam --output-fmt BAM
```
Note: The reference genome can be downloaded here - https://rapid.ensembl.org/Sciurus_carolinensis_GCA_902686445.2/Info/Index, or via the command-line:
```bash
wget http://ftp.ensembl.org/pub/rapid-release/species/Sciurus_carolinensis/GCA_902686445.2/genome/Sciurus_carolinensis-GCA_902686445.2-unmasked.fa.gz
```
Note: Will not align against the smaller unplaced scaffolds as some had crazy high coverages when I ran through this pipeline the first time (possibly reducing coverage on the chromosomes), so I used seqkit to remove scaffolds less than 1Mb in the genome file.
```bash
seqkit seq --min-len 1000000 egsq_genome.fa.gz > egsq_genome_1MBmin.fa.gz
```
Note: Aligning unpaired reads to the reference genome as well since there is a decent percentage (>10% of total cleaned reads) in the 3x samples. I first concatenated the forward (1U) and reverse (2U) unpaired reads to create one unpaired read file, then aligned. e.g., concatenation and single-end alignment:
```bash
cat SCCA1009_trimmed_1U.fastq.gz SCCA1009_trimmed_2U.fastq.gz > SCCA1009_trimmed_UC.fastq.gz

bwa mem -t 20 egsq_1MBmin SCCA1009_trimmed_UC.fastq.gz | \
  samtools sort --threads 20 -o SCCA1009_unpaired.bam --output-fmt BAM
```
2) Merge the unpaired and paired read alignments using samtools merge and check their quality using samtools, QualiMap, and MultiQC.
```bash
samtools merge --threads 20 \
  -o SCCA1009_merged.bam \
  SCCA1009_paired.bam \
  SCCA1009_unpaired.bam
  
qualimap bamqc -bam SCCA1009_merged.bam \
  -outdir SCCA1009_merged_qualimap \
  -nt 20 \
  --java-mem-size=110G
  
samtools idxstats --threads 20 SCCA1009_merged.bam > SCCA1009_merged_idxstats.txt
samtools flagstat --threads 20 SCCA1009_merged.bam > SCCA1009_merged_flagstat.txt
multiqc ./ --interactive
```
3) Mark duplicate aligned reads using Picard MarkDuplicates.
```bash
java -jar picard.jar MarkDuplicates \
  I=SCCA1009_merged.bam \
  O=SCCA1009_merged_dedup.bam \
  M=SCCA1009_merged_metrics.txt
```
4) Add read group header to the deduplicated bam files using Picard AddOrReplaceReadGroups.
```bash
java -jar picard.jar AddOrReplaceReadGroups \
  I=SCCA1009_merged_dedup.bam \
  O=SCCA1009_merged_dedup_rg.bam \
  RGID=HMJ7JDSX2.2 \
  RGPU=HMJ7JDSX2.2.SCCA1009 \
  RGSM=SCCA1009 \
  RGPL=ILLUMINA \
  RGLB=wgs_SCCA1009
```
5) Perform local realignment around indels using GATK and index bam files (then reassess final bam files with samtools, QualiMap, and MultiQC).
```bash
java -Xmx25g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R egsq_genome_1MBmin.fa \
  -I SCCA1009_merged_dedup_rg.bam \
  -o SCCA1009_merged_dedup_rg.intervals 

java -Xmx25g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R egsq_genome_1MBmin.fa \
  -targetIntervals SCCA1009_merged_dedup_rg.intervals \
  -I SCCA1009_merged_dedup_rg.bam \
  -o SCCA1009_merged_dedup_rg_realigned.bam 

samtools index -@ 20 SCCA1009_merged_dedup_rg.bam

qualimap bamqc -bam SCCA1009_merged_dedup_rg_realigned.bam \
  -outdir SCCA1009_merged_dedup_rg_realigned_qualimap \
  -nt 20 \
  --java-mem-size=110G

samtools idxstats --threads 20 SCCA1009_merged_dedup_rg_realigned.bam > SCCA1009_merged__dedup_rg_realignedidxstats.txt
samtools flagstat --threads 20 SCCA1009_merged_dedup_rg_realigned.bam > SCCA1009_merged__dedup_rg_realignedflagstat.txt
multiqc ./ --interactive
```

# SNP generation #

1) 
```bash

```
2) 
```bash

```
3)
```bash

```
4)
```bash

```
# Analyses #
1)
2)
3)

