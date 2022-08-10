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
2) Trim the raw reads of adapter contamination and low-quality bases using Trimmomatic (then reassess quality with FastQC and MultiQC).
3) Align the trimmed reads to the reference genome using BWA.
4) Check initial quality of the alignments using QualiMap and MultiQC.
5) Remove duplicate aligned reads using Picard (then verify that duplicates have been removed using QualiMap and MultiQC).
6) First step of the GATK SNP calling pipeline. Run HaplotypeCaller on the deduplicated aligned reads to generate initial variant calls (1 gvcf produced for each sample; 47 here).
7) Second step of the GATK SNP calling pipeline. Run GenomicsDBImport to combine the variants called across each sample (1 database produced for each chromosome; 21 here).
8) Third and final step of the GATK SNP calling pipeline. Run GenotypeGVCFs to jointly call variants for every sample across all chromosomes (1 gvcf produced for each chromosome; 21 here). After the gvcfs are generated, use Picard to combine all 21 into one large gvcf file.
9) Filter the variants called by GATK using bcftools.
10) Annotate SNPs

# Analyses (running) #
1)
2)
3)

