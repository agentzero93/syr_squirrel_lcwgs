# Genomic differentiation between urban and rural eastern gray squirrels #
Repository to house scripts used to analyze low-coverage whole genome sequencing of the eastern gray squirrel (Sciurus carolinensis)

# Tools used (running list) #
### Read data QC ###
1) FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2) MultiQC (https://multiqc.info/)
3) Trimmomatic (https://github.com/usadellab/Trimmomatic)
### Alignment ###
1) BWA (https://github.com/lh3/bwa)
2) samtools (http://www.htslib.org/)
3) QualiMap (http://qualimap.conesalab.org/)
4) MultiQC
### SNP generation ###
1) GATK (https://gatk.broadinstitute.org/hc/en-us)
2) Picard (https://github.com/broadinstitute/picard)
3) bcftools (https://samtools.github.io/bcftools/bcftools.html)
4) samtools
### Analysis ###
1)
2)
3)

# Data generation (only running Syracuse samples as of now) #
1) Check the quality of the raw reads (.fastq.gz) using FastQC and MultiQC.

Note: SCCA1009 SCCA1011 SCCA1012 SCCA1017 SCCA1018 SCCA1020 SCCA1026 SCCA1033 SCCA1034 SCCA1040 were run at 10x, the rest were run at 3x.

Note: Raw reads look good, just a bit of adapter contamination.

2) Trim the raw reads of adapter contamination and low-quality bases using Trimmomatic (then reassess quality with FastQC and MultiQC).

Note: Reads no longer have adapter contamination.

3) Align the trimmed reads to the reference genome using BWA.

Note: The reference genome can be downloaded here - https://rapid.ensembl.org/Sciurus_carolinensis_GCA_902686445.2/Info/Index

Note: Aliging unpaired reads to the reference genome as well since there is a decent percentage (>10% of total cleaned reads) in the 3x samples. I first concatenated the forward (1U) and reverse (2U) unpaired reads to create one unpaired read file (cat 1U.fastq.gz 2U.fastq.gz > UC.fastq.gz), then aligned.

Note: To parallelize, run across multiple nodes with as many threads as possible on your computer(s). Do not try to run multiple alignments on a single node, this will result in a major slow down (see https://www.biostars.org/p/420062/). 

4) Check initial quality of the alignments using QualiMap and MultiQC.
5) Mark and remove duplicate aligned reads using Picard MarkDuplicates (then verify that duplicates have been removed using QualiMap and MultiQC).

Note: I merged the unpaired and paired read alignments (*_merged.bam) prior to deduplicating the reads using samtools merge.

Note: Read groups got messed up after merging, so had to use Picard AddOrReplaceReadGroups to fix the read groups.

6) First step of the GATK SNP calling pipeline. Run HaplotypeCaller on the deduplicated aligned reads to generate initial variant calls (1 gvcf produced for each sample; 47 here).
7) Second step of the GATK SNP calling pipeline. Run GenomicsDBImport to combine the variants called across each sample (1 database produced for each chromosome; 21 here).
8) Third and final step of the GATK SNP calling pipeline. Run GenotypeGVCFs to jointly call variants for every sample across all chromosomes (1 gvcf produced for each chromosome; 21 here). After the gvcfs are generated, use Picard to combine all 21 into one large gvcf file.
9) Filter the variants called by GATK using bcftools.
10) Annotate SNPs

# Analyses (running) #
1)
2)
3)

