
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

