# Genomic differentiation between urban and rural eastern gray squirrels in Syracuse, NY #
Repository to house scripts used to analyze low-coverage whole genome sequencing of the eastern gray squirrel (Sciurus carolinensis) in Syracuse, New York.

# Tools used (running list) #
### Read data QC ###
1) FastQC 0.11.9 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2) MultiQC 1.10.1 (https://multiqc.info/)
3) Trimmomatic 0.39 (https://github.com/usadellab/Trimmomatic)
### Alignment ###
1) BWA 0.7.17 (https://github.com/lh3/bwa)
2) seqkit 2.3.0 (https://github.com/shenwei356/seqkit)
3) samtools 1.13 (http://www.htslib.org/)
4) QualiMap (http://qualimap.conesalab.org/)
5) MultiQC
### SNP generation ###
1) GATK 4.2.6.1 (https://gatk.broadinstitute.org/hc/en-us)
2) Picard 2.25.6 (https://github.com/broadinstitute/picard)
3) bcftools 1.16 (https://samtools.github.io/bcftools/bcftools.html)
4) samtools
### Analysis ###
1)
2)
3)

# Data generation (only running Syracuse samples as of now) #

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
```
Note: Reads no longer have adapter contamination.
3) Align the trimmed reads to the reference genome using BWA. e.g., paired-end alignment:
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
seqkit seq --min-len 1000000 genome.fa.gz > genome_1MBmin.fa.gz
```
Note: Aligning unpaired reads to the reference genome as well since there is a decent percentage (>10% of total cleaned reads) in the 3x samples. I first concatenated the forward (1U) and reverse (2U) unpaired reads to create one unpaired read file, then aligned. e.g., concatenation and single-end alignment:
```bash
cat SCCA1009_trimmed_1U.fastq.gz SCCA1009_trimmed_2U.fastq.gz > SCCA1009_trimmed_UC.fastq.gz

bwa mem -t 20 egsq_1MBmin SCCA1009_trimmed_UC.fastq.gz | \
  samtools sort --threads 20 -o SCCA1009_unpaired.bam --output-fmt BAM
```
4) Merge the unpaired and paired read alignments using samtools merge and check their quality using QualiMap and MultiQC.
```bash
samtools merge --threads 20 \
  -o SCCA1009_merged.bam \
  SCCA1009_paired.bam \
  SCCA1009_unpaired.bam
  
qualimap bamqc -bam SCCA1009_merged.bam \
  -outdir SCCA1009_merged_qualimap \
  -nt 20 \
  --java-mem-size=110G
  
multiqc ./ --interactive
```
5) Mark duplicate aligned reads using Picard MarkDuplicates (then reassess bam files with QualiMap and MultiQC).
```bash
java -jar picard.jar MarkDuplicates \
  I=SCCA1009_merged.bam \
  O=SCCA1009_merged_dedup.bam \
  M=SCCA1009_merged_metrics.txt
```
6) Add read group header to the deduplicated bam files using Picard AddOrReplaceReadGroups, then index bam files with samtools index (required by GATK).
```bash
java -jar picard.jar AddOrReplaceReadGroups \
  I=SCCA1009_merged_dedup.bam \
  O=SCCA1009_merged_dedup_rg.bam \
  RGID=HMJ7JDSX2.2 \
  RGPU=HMJ7JDSX2.2.SCCA1009 \
  RGSM=SCCA1009 \
  RGPL=ILLUMINA \
  RGLB=wgs_SCCA1009
  
samtools index -@ 20 SCCA1009_merged_dedup_rg.bam
```
7) First step of the GATK SNP calling pipeline. Run HaplotypeCaller on the deduplicated aligned reads to generate initial variant calls (1 gvcf produced for each sample; 46 here).
```bash
gatk --java-options '-Xmx100g' HaplotypeCaller \
  --input SCCA1009_merged_dedup_rg.bam \
  --output SCCA1009_merged_dedup_rg.gvcf.gz \
  --reference genome_1MBmin.fa.gz \
  -ERC GVCF
```
Note: Only running this pipeline on the 46 Syracuse samples for now. Once the other cities are sequenced, I will run those samples through the pipeline as well.

8) Second step of the GATK SNP calling pipeline. Run GenomicsDBImport to combine the variants called across each sample (1 database produced for each chromosome and unplaced scaffold; 89 here).
```bash

```
9) Third and final step of the GATK SNP calling pipeline. Run GenotypeGVCFs to jointly call variants for every sample across all chromosomes and unplaced scaffolds (1 gvcf produced for each chromosome and unplaced scaffold; 89 here). After the gvcfs are generated, use Picard to combine all 89 files into one large gvcf file.
```bash

```
10) Filter the variants called by GATK using bcftools.
11) Annotate SNPs

# Analyses (running) #
1)
2)
3)

