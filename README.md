________________________________________
**🧬 DNAseq Analysis: Lung Adenocarcinoma (ADC) & Squamous Cell Carcinoma (SCC)**
This project analyzes lung adenocarcinoma (ADC) and squamous cell carcinoma (SCC) using Next Generation Sequencing (NGS) data to identify somatic and germline variants.
Key steps include quality control (FastQC, fastp), alignment (BWA), variant calling (GATK), and annotation (Ensembl VEP).
The analysis highlights:
•	Recurrent TRBV mutations in SCC (immune modulation)
•	Germline PTCH2 variant in ADC (Hedgehog pathway)
•	Subtype-specific GPRIN2/3 germline variants in SCC
Overall, the study demonstrates molecular heterogeneity, identifies potential biomarkers, and provides population-specific insights for precision oncology.
________________________________________
**📚 Table of Contents**
1.	Overview
2.	Workflow Summary
3.	Data Preparation
4.	Quality Control
5.	Trimming
6.	Alignment
7.	Post-processing
8.	Variant Calling
9.	Variant Filtering
10.	Annotation
11.	Visualization
12.	Structural Variant Analysis
13.	Findings
________________________________________
**🧾 Overview**
Category	Tool / Resource	Purpose
Quality Control	FastQC, fastp	Assess and improve read quality
Alignment	BWA-MEM	Align reads to reference genome
File Processing	SAMtools, Picard	Convert, sort, and deduplicate reads
Variant Calling	GATK HaplotypeCaller, Mutect2	Identify germline and somatic variants
Variant Filtering	SnpSift	Retain high-confidence variants
Annotation	Ensembl VEP	Predict biological and clinical significance
Visualization	IGV	Manual inspection of alignments
Structural Variants	BreakDancer	Detect large genomic rearrangements
________________________________________
**⚙️ Workflow Summary**
1. Raw Data Download (SRA)
2. Reference Genome Preparation (hg38)
3. Quality Control (FastQC)
4. Read Trimming and Filtering (fastp)
5. Sequence Alignment (BWA-MEM)
6. File Conversion & Sorting (samtools)
7. Duplicate Removal & Read Grouping (Picard)
8. Variant Calling (GATK)
9. Variant Filtering (SnpSift)
10. Variant Annotation (VEP)
11. Visualization (IGV)
12. Structural Variant Analysis (BreakDancer)
________________________________________
**📥 Raw Data Download**
Download tumor and normal FASTQ files from NCBI SRA:
# Normal sample
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/006/SRR20761706/SRR20761706_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/006/SRR20761706/SRR20761706_2.fastq.gz

# Squamous Cell Carcinoma (SCC)
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/011/SRR20761711/SRR20761711_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/011/SRR20761711/SRR20761711_2.fastq.gz

# Adenocarcinoma
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/004/SRR20761704/SRR20761704_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/004/SRR20761704/SRR20761704_2.fastq.gz
________________________________________
🧬 Reference Genome Preparation
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar -xzvf hg38.chromFa.tar.gz
cat chr*.fa > genome.fa
bwa index genome.fa
________________________________________
**🔍 Quality Control (FastQC)**
# Run FastQC
fastqc SRR20761706_1.fastq.gz SRR20761706_2.fastq.gz
fastqc SRR20761711_1.fastq.gz SRR20761711_2.fastq.gz
fastqc SRR20761704_1.fastq.gz SRR20761704_2.fastq.gz
Metric	Description
Per-base quality	Mean Phred score across read positions
Adapter content	Presence of adapter sequences
GC content	Expected vs observed GC percentage
Sequence duplication	Estimate of duplicated reads
________________________________________
**✂️ Trimming (fastp)**
# Normal
fastp -i SRR20761706_1.fastq.gz -I SRR20761706_2.fastq.gz \
-o trim.SRR20761706_1.fastq.gz -O trim.SRR20761706_2.fastq.gz -q 30

# SCC
fastp -i SRR20761711_1.fastq.gz -I SRR20761711_2.fastq.gz \
-o trim.SRR20761711_1.fastq.gz -O trim.SRR20761711_2.fastq.gz -q 30

# ADC
fastp -i SRR20761704_1.fastq.gz -I SRR20761704_2.fastq.gz \
-o trim.SRR20761704_1.fastq.gz -O trim.SRR20761704_2.fastq.gz -q 30
________________________________________
**🧩 Alignment (BWA-MEM)**
bwa mem genome.fa trim.SRR20761706_1.fastq.gz trim.SRR20761706_2.fastq.gz > SRR20761706.sam  # Normal
bwa mem genome.fa trim.SRR20761711_1.fastq.gz trim.SRR20761711_2.fastq.gz > SRR20761711.sam  # SCC
bwa mem genome.fa trim.SRR20761704_1.fastq.gz trim.SRR20761704_2.fastq.gz > SRR20761704.sam  # ADC
________________________________________
**🧮 SAM to BAM Conversion & Sorting**
samtools view -u SRR20761706.sam | samtools sort -o sort.SRR20761706.bam
samtools view -u SRR20761711.sam | samtools sort -o sort.SRR20761711.bam
samtools view -u SRR20761704.sam | samtools sort -o sort.SRR20761704.bam
________________________________________
**🧹 Remove Duplicates (Picard)**
picard-tools MarkDuplicates I=sort.SRR20761706.bam O=rmdup.sort.SRR20761706.bam M=dup_metrics.txt
picard-tools AddOrReplaceReadGroups I=rmdup.sort.SRR20761706.bam O=reheader.rmdup.sort.SRR20761706.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
samtools index reheader.rmdup.sort.SRR20761706.bam
________________________________________
**🧠 Variant Calling (GATK)**
Germline variants (HaplotypeCaller):
gatk HaplotypeCaller -R genome.fa -I reheader.rmdup.sort.SRR20761706.bam -O GATK_germline.g.vcf.gz
Somatic variants (Mutect2):
gatk Mutect2 -R genome.fa -I reheader.rmdup.sort.SRR20761711.bam \
--germline-resource af-only-gnomad.hg38.vcf.gz \
--panel-of-normals GermlineHetPon.38.vcf.gz -O somatic_variation.vcf.gz
________________________________________
**🧪 Variant Filtering (SnpSift)**
# Germline
cat GATK_germline.g.vcf | java -jar SnpSift.jar filter "((QUAL >= 30) & (MQ >= 30) & (DP >= 10))" > Filtered.GATK_germline.g.vcf

# Somatic
cat somatic_variation.vcf | java -jar SnpSift.jar filter "(DP >= 10)" > Filtered.somatic_variation.vcf
________________________________________
**🧬 Variant Annotation (Ensembl VEP)**
Upload the filtered VCFs to Ensembl VEP for annotation.
Annotation Field	Description
SIFT / PolyPhen	Functional prediction of variant impact
ClinVar	Known clinical associations
COSMIC	Cancer-specific mutation database
gnomAD AF	Population allele frequency
________________________________________
**🧭 Visualization (IGV)**
samtools index reheader.rmdup.sort.SRR20761706.bam
igv.sh
# Load indexed BAM and corresponding VCF files for visual validation
________________________________________
**🧱 Structural Variant Analysis (BreakDancer)**
BreakDancer-Max -g genome.fa -o output_breakdancer.txt reheader.rmdup.sort.SRR20761711.bam
________________________________________
**🔬 Findings**
Subtype	Variant Type	Gene	Pathway / Role
SCC	Somatic	TRBV	Immune modulation
ADC	Germline	PTCH2	Hedgehog signaling
SCC	Germline	GPRIN2/3	Neuronal signaling, subtype-specific
________________________________________
**🧠 Conclusion**
This integrated DNAseq workflow provides:
•	Reliable detection of both somatic and germline variants
•	Identification of subtype-specific molecular markers
•	A foundation for precision oncology through bioinformatics-driven discovery

