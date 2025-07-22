# Complete-Variant-Calling-Analysis-Pipeline
This comprehensive pipeline performs variant calling from paired-end whole genome sequencing data and includes downstream variant analysis and visualization.


**Pipeline Highlights:**

•	🔍 Quality control of raw reads using FastQC

•	🔗 Alignment to reference genome via BWA-MEM

•	🧹 Sorting, duplicate marking, and indexing with Samtools and Picard

•	🎯 Base Quality Score Recalibration (BQSR) using GATK

•	🔎 Variant detection using GATK HaplotypeCaller

•	🚫 Filtering low-quality variants based on multiple metrics

•	📊 Downstream variant summary stats and plots generated in R:

•	Variant quality score distribution

•	Allele frequency histogram

•	Variant type (SNP vs Indel) bar plot

•	💾 Saves all intermediate and final files for reproducibility and further analysis

Example data and setup instructions are provided for easy testing and adaptation.

