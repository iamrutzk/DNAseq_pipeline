
# DNAseq Analysis

# Key Features

Cancer Types: Lung adenocarcinoma (ADC), squamous cell carcinoma (SCC)

Variants Detected: Germline and somatic (GATK HaplotypeCaller, Mutect2)

Annotation: Ensembl VEP for biological interpretation (SIFT, PolyPhen)

Findings: TRBV mutations in SCC, germline PTCH2 in ADC, subtype-specific GPRIN2/3 in SCC

Bioinformatics Integration: Workflow enables translational research for targeted therapies and population-specific insights.

# Table of Contents:

Run these step by step:

1. Raw data Download :
    - [Normal Sample](#normal-sample)
    - [Squamous Cell Carcinoma (SCC)](#squamous-cell-carcinoma-scc)
    - [Adenocarcinoma](#adenocarcinoma)
2. Reference Genome Preparation : - [Ref Genome](#ref-genome)
3. Quality Control Analysis : - [Run FastQC on Raw Data](#run-fastqc-on-raw-data)
4. Read Trimming and Filtering :
    - [Normal](#normal)
    - [SCC](#scc)
    - [Adenocarcinoma](#adenocarcinoma-1)
5. Run FASTQC : - [Re-run FastQC on Trimmed Data](#re-run-fastqc-on-trimmed-data)
6. Sequence Alignment
    - [Index Reference Genome](#index-reference-genome)
    - [Align Reads](#align-reads)
7. File Format Conversion :
    - [Normal](#normal-1)
    - [SCC](#scc-1)
    - [Adenocarcinoma](#adenocarcinoma-2) 
8. Duplicate Removal
    - [Normal](#normal-2)
    - [SCC](#scc-2)
    - [Adenocarcinoma](#adenocarcinoma-3)
9. Variant Calling : 
    - [Download Resources for Mutect2](#download-required-resources-for-mutect2)
    - [Germline Variants (HaplotypeCaller)](#germline-variants-haplotypecaller)
    - [Somatic Variants (Mutect2)](#somatic-variants-mutect2) 
10. Variant Filtering and Annotation : 
    - [Germline Variants](#germline-variants)
    - [Somatic Variants](#somatic-variants)  
11. Data Visualization : 
    - [Load BAM and VCF Files into IGV](#load-bam-and-vcf-files-into-igv-for-manual-inspection)
12. Structural Variant Analysis : 
    - [Output Data](#output-data)

# DNAseq Pipeline for Lung Cancer Variant Analysis:

This database is taken for a lung cancer patients, suffering from different type of cancer associated with
lungs cancer. The data was downloaded using NCBI SRA: https://www.ncbi.nlm.nih.gov/sra

1. Download Raw Data (FASTQ) from NCBI SRA
The FASTQ files were downloaded from NCBI SRA database. Both tumor and normal samples were acquired for comparison. The data included squamous cell carcinoma and adenocarcinoma subtypes. The download was performed using wget command with resume capability. Proper file naming conventions were maintained throughout.

# Normal sample
```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/006/SRR20761706/SRR20761706_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/006/SRR20761706/SRR20761706_2.fastq.gz
```
# Squamous Cell Carcinoma (SCC)
```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/011/SRR20761711/SRR20761711_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/011/SRR20761711/SRR20761711_2.fastq.gz
```
# Adenocarcinoma 
```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/004/SRR20761704/SRR20761704_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/004/SRR20761704/SRR20761704_2.fastq.gz
```
2. Download Reference Genome (hg38)
The hg38 human reference genome was obtained from UCSC database. Individual chromosome files were combined into a single genome file. The reference was indexed for efficient alignment. This provided the baseline for all variant calling. The genome preparation ensured compatibility with downstream tools.

# Ref Genome 
```
wget -c https://hgdownload.soc.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar -xzvf hg38.chromFa.tar.gz
```
```
cat chr*.fa > genome.fa
````
3. Quality Control (FastQC)
FastQC was run on all raw sequencing files. Quality metrics were assessed including per-base sequencing quality. Adapter contamination and sequence duplication levels were checked. The reports helped identify problematic samples. This guided the subsequent trimming parameters.

# Run FastQC on raw data
```
fastqc SRR20761706_1.fastq.gz SRR20761706_2.fastq.gz  # Normal
fastqc SRR20761711_1.fastq.gz SRR20761711_2.fastq.gz  # SCC
fastqc SRR20761704_1.fastq.gz SRR20761704_2.fastq.gz  # Adenocarcinoma
```
4. Trimming (fastp)
Low-quality bases and adapters were removed using fastp. Reads below quality threshold of 30 were discarded. Both 5' and 3' ends of reads were evaluated. The trimmed files showed improved quality metrics. FastQC was rerun to verify trimming effectiveness.

# Normal
```
fastp -i SRR20761706_1.fastq.gz -o trim.SRR20761706_1.fastq.gz -I SRR20761706_2.fastq.gz -O trim.SRR20761706_2.fastq.gz -q 30 -t 2 -T 2
```
# SCC
```
fastp -i SRR20761711_1.fastq.gz -o trim.SRR20761711_1.fastq.gz -I SRR20761711_2.fastq.gz -O trim.SRR20761711_2.fastq.gz -q 30 -t 2 -T 2
```
# Adenocarcinoma
```
fastp -i SRR20761704_1.fastq.gz -o trim.SRR20761704_1.fastq.gz -I SRR20761704_2.fastq.gz -O trim.SRR20761704_2.fastq.gz -q 30 -t 2 -T 2
```
# Re-run FastQC on trimmed data
```
fastqc trim.SRR20761706_1.fastq.gz trim.SRR20761706_2.fastq.gz
fastqc trim.SRR20761711_1.fastq.gz trim.SRR20761711_2.fastq.gz
fastqc trim.SRR20761704_1.fastq.gz trim.SRR20761704_2.fastq.gz
```

5. Alignment (BWA-MEM)
Reads were aligned to reference genome using BWA-MEM. The alignment process generated SAM format files. Both tumor and normal samples were processed identically. Proper read group information was incorporated. The alignments were optimized for variant detection.

# Index reference genome
```
bwa index genome.fa
```
# Align reads
```
bwa mem genome.fa trim.SRR20761706_1.fastq.gz trim.SRR20761706_2.fastq.gz > SRR20761706.sam  # Normal
bwa mem genome.fa trim.SRR20761711_1.fastq.gz trim.SRR20761711_2.fastq.gz > SRR20761711.sam  # SCC
bwa mem genome.fa trim.SRR20761704_1.fastq.gz trim.SRR20761704_2.fastq.gz > SRR20761704.sam  # Adenocarcinoma
```
6. SAM to BAM Conversion & Sorting (samtools)
SAM files were converted to sorted BAM format. Duplicate reads were marked using Picard tools. Read groups were added to facilitate analysis. The files were indexed for rapid access. This prepared the data for variant calling.

# Normal
```
samtools view -u SRR20761706.sam | samtools sort -o sort.SRR20761706.bam
```
# SCC
```
samtools view -u SRR20761711.sam | samtools sort -o sort.SRR20761711.bam
```
# Adenocarcinoma
```
samtools view -u SRR20761704.sam | samtools sort -o sort.SRR20761704.bam
```
7. Remove Duplicates (Picard)
PCR duplicates were identified and removed. MarkDuplicates tool from Picard was utilized. Metrics were generated to assess duplication rates. The process improved variant calling accuracy. Final BAM files were ready for analysis.

# Normal
```
picard-tools MarkDuplicates I=sort.SRR20761706.bam O=rmdup.sort.SRR20761706.bam M=marked_dup_metrics.txt
picard-tools AddOrReplaceReadGroups I=rmdup.sort.SRR20761706.bam O=reheader.rmdup.sort.SRR20761706.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
samtools index reheader.rmdup.sort.SRR20761706.bam
```
# SCC
```
picard-tools MarkDuplicates I=sort.SRR20761711.bam O=rmdup.sort.SRR20761711.bam M=marked_dup_metrics.txt
picard-tools AddOrReplaceReadGroups I=rmdup.sort.SRR20761711.bam O=reheader.rmdup.sort.SRR20761711.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
samtools index reheader.rmdup.sort.SRR20761711.bam
```
# Adenocarcinoma
```
picard-tools MarkDuplicates I=sort.SRR20761704.bam O=rmdup.sort.SRR20761704.bam M=marked_dup_metrics.txt
picard-tools AddOrReplaceReadGroups I=rmdup.sort.SRR20761704.bam O=reheader.rmdup.sort.SRR20761704.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
samtools index reheader.rmdup.sort.SRR20761704.bam
```

8. Variant Calling (GATK)
Both germline and somatic variants were identified. GATK HaplotypeCaller detected germline mutations. Mutect2 was used for somatic variant detection. Strict quality filters were applied throughout. The process generated comprehensive VCF files.

# Download required resources for Mutect2
```
wget -c https://www.bcgsc.ca/downloads/morinlab/reference/af-only-gnomad.hg38.vcf.gz
wget -c https://www.bcgsc.ca/downloads/morinlab/reference/af-only-gnomad.hg38.vcf.gz.tbi
wget -c https://www.bcgsc.ca/downloads/morinlab/hmftools-references/amber/GermlineHetPon.38.vcf.gz
tabix GermlineHetPon.38.vcf.gz
```
# Germline variants (HaplotypeCaller)
```
# Normal
java -jar gatk-package-4.6.0.0-local.jar HaplotypeCaller -R genome.fa -I reheader.rmdup.sort.SRR20761706.bam -O GATK_germline.g.vcf.gz
# SCC
java -jar gatk-package-4.6.0.0-local.jar HaplotypeCaller -R genome.fa -I reheader.rmdup.sort.SRR20761711.bam -O GATK_germline.g.vcf.gz
# Adenocarcinoma
java -jar gatk-package-4.6.0.0-local.jar HaplotypeCaller -R genome.fa -I reheader.rmdup.sort.SRR20761704.bam -O GATK_germline.g.vcf.gz  
```
# Somatic variants (Mutect2)
```
# SCC
java -jar gatk-package-4.6.0.0-local.jar Mutect2 -R genome.fa -I reheader.rmdup.sort.SRR20761711.bam --germline-resource af-only-gnomad.hg38.vcf.gz --panel-of-normals GermlineHetPon.38.vcf.gz -O somatic_variation.vcf.gz  
```
```
# Adenocarcinoma
java -jar gatk-package-4.6.0.0-local.jar Mutect2 -R genome.fa -I reheader.rmdup.sort.SRR20761704.bam --germline-resource af-only-gnomad.hg38.vcf.gz --panel-of-normals GermlineHetPon.38.vcf.gz -O somatic_variation.vcf.gz  
```
9. Variant Filtering (SnpSift)
Variants were filtered based on quality metrics. Minimum depth of 10x was required. Mapping quality and base quality thresholds were enforced. SnpSift was used for efficient filtering. Only high-confidence variants were retained.

# Germline variants
```
cat GATK_germline.g.vcf | java -jar snpEff/SnpSift.jar filter "((QUAL >= 30) & (MQ >= 30) & (DP >= 10))" > Filtered.GATK_germline.g.vcf
```
# Somatic variants
```
cat somatic_variation.vcf | java -jar snpEff/SnpSift.jar filter "(DP >= 10)" > Filtered.somatic_variation.vcf
```
10. Annotation (Ensembl VEP)
Variants were annotated using Ensembl VEP. Functional consequences were predicted. Known variants were cross-referenced with databases. Pathogenicity scores (SIFT, PolyPhen) were calculated. The annotation added biological context to variants.
1. Upload filtered VCF files to Ensembl VEP web tool for annotation
2. Output will include functional predictions (SIFT, PolyPhen-2) and clinical significance

11. Visualization (IGV)
IGV was used to visualize alignments and variants. BAM files were properly indexed first. Specific genomic regions could be examined in detail. The visual confirmation supported variant validation. Screenshots were captured for documentation.

# Load BAM and VCF files into IGV for manual inspection
```
# Normal
samtools index reheader.rmdup.sort.SRR20761706.bam 
# SCC
samtools index reheader.rmdup.sort.SRR20761711.bam  
# Adenocarcinoma
samtools index reheader.rmdup.sort.SRR20761704.bam
```

12. Structural Variant Analysis (BreakDancer)
BreakDancer detected larger genomic alterations. Both deletions and amplifications were identified. The structural variants were annotated for significance. This complemented the small variant analysis. The complete variant profile was thus established.
```
# Example for SCC
BreakDancer-Max -g genome.fa -o output_breakdancer.txt reheader.rmdup.sort.SRR20761711.bam
```  
# Output Data
Variant Files: Filtered VCFs (germline/somatic)

QC Reports: FastQC and fastp HTML reports

Annotated Variants: VEP output files (SIFT, PolyPhen, clinical context)

Visual Evidence: IGV screenshot records
