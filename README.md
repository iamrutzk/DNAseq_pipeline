
# DNAseq Analysis

This project analyzes lung adenocarcinoma (ADC) and squamous cell carcinoma (SCC) using NGS data to identify somatic and germline variants. Key steps include quality control (FastQC, fastp), alignment (BWA), variant calling (GATK), and annotation (Ensembl VEP). Findings reveal recurrent TRBV mutations in SCC (immune modulation), a germline PTCH2 variant in ADC (Hedgehog pathway), and subtype-specific *GPRIN2/3* germline variants in SCC. The study highlights molecular heterogeneity, proposes potential biomarkers, and demonstrates population-specific insights for precision oncology. The workflow integrates bioinformatics pipelines with translational research, aiding targeted therapy development.

Table of Contents:

Run these step by step:

1. Raw Data Download  
2. Reference Genome Preparation  
3. Quality Control Analysis  
4. Read Trimming and Filtering  
5. Sequence Alignment  
6. File Format Conversion  
7. Duplicate Removal  
8. Variant Calling  
9. Variant Filtering  
10. Variant Annotation  
11. Data Visualization  
12. Structural Variant Analysis

DNAseq Pipeline for Lung Cancer Variant Analysis:

1. Download Raw Data (FASTQ) from NCBI SRA

# Normal sample
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/006/SRR20761706/SRR20761706_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/006/SRR20761706/SRR20761706_2.fastq.gz

# Squamous Cell Carcinoma (SCC)
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/011/SRR20761711/SRR20761711_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/011/SRR20761711/SRR20761711_2.fastq.gz

# Adenocarcinoma
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/004/SRR20761704/SRR20761704_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/004/SRR20761704/SRR20761704_2.fastq.gz

2. Download Reference Genome (hg38)

wget -c https://hgdownload.soc.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar -xzvf hg38.chromFa.tar.gz
cat chr*.fa > genome.fa

3. Quality Control (FastQC)

# Run FastQC on raw data
fastqc SRR20761706_1.fastq.gz SRR20761706_2.fastq.gz  # Normal
fastqc SRR20761711_1.fastq.gz SRR20761711_2.fastq.gz  # SCC
fastqc SRR20761704_1.fastq.gz SRR20761704_2.fastq.gz  # Adenocarcinoma

4. Trimming (fastp)

# Normal
fastp -i SRR20761706_1.fastq.gz -o trim.SRR20761706_1.fastq.gz -I SRR20761706_2.fastq.gz -O trim.SRR20761706_2.fastq.gz -q 30 -t 2 -T 2

# SCC
fastp -i SRR20761711_1.fastq.gz -o trim.SRR20761711_1.fastq.gz -I SRR20761711_2.fastq.gz -O trim.SRR20761711_2.fastq.gz -q 30 -t 2 -T 2

# Adenocarcinoma
fastp -i SRR20761704_1.fastq.gz -o trim.SRR20761704_1.fastq.gz -I SRR20761704_2.fastq.gz -O trim.SRR20761704_2.fastq.gz -q 30 -t 2 -T 2

# Re-run FastQC on trimmed data
fastqc trim.SRR20761706_1.fastq.gz trim.SRR20761706_2.fastq.gz
fastqc trim.SRR20761711_1.fastq.gz trim.SRR20761711_2.fastq.gz
fastqc trim.SRR20761704_1.fastq.gz trim.SRR20761704_2.fastq.gz

5. Alignment (BWA-MEM)

# Index reference genome
bwa index genome.fa

# Align reads
bwa mem genome.fa trim.SRR20761706_1.fastq.gz trim.SRR20761706_2.fastq.gz > SRR20761706.sam  # Normal
bwa mem genome.fa trim.SRR20761711_1.fastq.gz trim.SRR20761711_2.fastq.gz > SRR20761711.sam  # SCC
bwa mem genome.fa trim.SRR20761704_1.fastq.gz trim.SRR20761704_2.fastq.gz > SRR20761704.sam  # Adenocarcinoma

6. SAM to BAM Conversion & Sorting (samtools)

# Normal
samtools view -u SRR20761706.sam | samtools sort -o sort.SRR20761706.bam

# SCC
samtools view -u SRR20761711.sam | samtools sort -o sort.SRR20761711.bam

# Adenocarcinoma
samtools view -u SRR20761704.sam | samtools sort -o sort.SRR20761704.bam

7. Remove Duplicates (Picard)

# Normal
picard-tools MarkDuplicates I=sort.SRR20761706.bam O=rmdup.sort.SRR20761706.bam M=marked_dup_metrics.txt
picard-tools AddOrReplaceReadGroups I=rmdup.sort.SRR20761706.bam O=reheader.rmdup.sort.SRR20761706.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
samtools index reheader.rmdup.sort.SRR20761706.bam

# SCC
picard-tools MarkDuplicates I=sort.SRR20761711.bam O=rmdup.sort.SRR20761711.bam M=marked_dup_metrics.txt
picard-tools AddOrReplaceReadGroups I=rmdup.sort.SRR20761711.bam O=reheader.rmdup.sort.SRR20761711.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
samtools index reheader.rmdup.sort.SRR20761711.bam

# Adenocarcinoma
picard-tools MarkDuplicates I=sort.SRR20761704.bam O=rmdup.sort.SRR20761704.bam M=marked_dup_metrics.txt
picard-tools AddOrReplaceReadGroups I=rmdup.sort.SRR20761704.bam O=reheader.rmdup.sort.SRR20761704.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
samtools index reheader.rmdup.sort.SRR20761704.bam

8. Variant Calling (GATK)

# Download required resources for Mutect2
wget -c https://www.bcgsc.ca/downloads/morinlab/reference/af-only-gnomad.hg38.vcf.gz
wget -c https://www.bcgsc.ca/downloads/morinlab/reference/af-only-gnomad.hg38.vcf.gz.tbi
wget -c https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.38.vcf.gz
tabix GermlineHetPon.38.vcf.gz

# Germline variants (HaplotypeCaller)
java -jar gatk-package-4.6.0.0-local.jar HaplotypeCaller -R genome.fa -I reheader.rmdup.sort.SRR20761706.bam -O GATK_germline.g.vcf.gz  # Normal
java -jar gatk-package-4.6.0.0-local.jar HaplotypeCaller -R genome.fa -I reheader.rmdup.sort.SRR20761711.bam -O GATK_germline.g.vcf.gz  # SCC
java -jar gatk-package-4.6.0.0-local.jar HaplotypeCaller -R genome.fa -I reheader.rmdup.sort.SRR20761704.bam -O GATK_germline.g.vcf.gz  # Adenocarcinoma

# Somatic variants (Mutect2)
java -jar gatk-package-4.6.0.0-local.jar Mutect2 -R genome.fa -I reheader.rmdup.sort.SRR20761711.bam --germline-resource af-only-gnomad.hg38.vcf.gz --panel-of-normals GermlineHetPon.38.vcf.gz -O somatic_variation.vcf.gz  # SCC
java -jar gatk-package-4.6.0.0-local.jar Mutect2 -R genome.fa -I reheader.rmdup.sort.SRR20761704.bam --germline-resource af-only-gnomad.hg38.vcf.gz --panel-of-normals GermlineHetPon.38.vcf.gz -O somatic_variation.vcf.gz  # Adenocarcinoma

9. Variant Filtering (SnpSift)

# Germline variants
cat GATK_germline.g.vcf | java -jar snpEff/SnpSift.jar filter "((QUAL >= 30) & (MQ >= 30) & (DP >= 10))" > Filtered.GATK_germline.g.vcf

# Somatic variants
cat somatic_variation.vcf | java -jar snpEff/SnpSift.jar filter "(DP >= 10)" > Filtered.somatic_variation.vcf

10. Annotation (Ensembl VEP)

Upload filtered VCF files to Ensembl VEP web tool for annotation
Output will include functional predictions (SIFT, PolyPhen-2) and clinical significance

11. Visualization (IGV)

# Load BAM and VCF files into IGV for manual inspection
samtools index reheader.rmdup.sort.SRR20761706.bam  # Normal
samtools index reheader.rmdup.sort.SRR20761711.bam  # SCC
samtools index reheader.rmdup.sort.SRR20761704.bam  # Adenocarcinoma

12. Structural Variant Analysis (BreakDancer)

BreakDancer-Max -g genome.fa -o output_breakdancer.txt reheader.rmdup.sort.SRR20761711.bam  # Example for SCC

