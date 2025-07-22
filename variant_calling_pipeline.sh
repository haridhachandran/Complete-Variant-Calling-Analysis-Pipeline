#!/bin/bash

# ================================
# Variant Calling Pipeline
# Usage: bash variant_calling_pipeline.sh sample_name R1.fastq.gz R2.fastq.gz reference.fa
# ================================

# Arguments
SAMPLE=$1
R1=$2
R2=$3
REFERENCE=$4

# Output directories
OUTDIR="results/${SAMPLE}"
mkdir -p "${OUTDIR}"

# Load modules or set paths (adjust if needed)
# module load fastqc bwa samtools gatk bcftools

echo "Step 1: Quality control (FastQC)"
fastqc -o "${OUTDIR}" "$R1" "$R2"

echo "Step 2: Read alignment (BWA-MEM)"
bwa mem -t 4 "${REFERENCE}" "$R1" "$R2" | samtools view -bS - > "${OUTDIR}/${SAMPLE}.bam"

echo "Step 3: Sort and index BAM"
samtools sort -@ 4 -o "${OUTDIR}/${SAMPLE}.sorted.bam" "${OUTDIR}/${SAMPLE}.bam"
samtools index "${OUTDIR}/${SAMPLE}.sorted.bam"

echo "Step 4: Mark duplicates (Picard)"
picard MarkDuplicates I="${OUTDIR}/${SAMPLE}.sorted.bam" O="${OUTDIR}/${SAMPLE}.dedup.bam" \
      M="${OUTDIR}/${SAMPLE}.dedup.metrics" REMOVE_DUPLICATES=true
samtools index "${OUTDIR}/${SAMPLE}.dedup.bam"

echo "Step 5: Base quality score recalibration (BQSR)"
gatk BaseRecalibrator -R "${REFERENCE}" -I "${OUTDIR}/${SAMPLE}.dedup.bam" --known-sites dbsnp.vcf \
      -O "${OUTDIR}/${SAMPLE}.recal.table"
gatk ApplyBQSR -R "${REFERENCE}" -I "${OUTDIR}/${SAMPLE}.dedup.bam" \
      --bqsr-recal-file "${OUTDIR}/${SAMPLE}.recal.table" -O "${OUTDIR}/${SAMPLE}.recal.bam"
samtools index "${OUTDIR}/${SAMPLE}.recal.bam"

echo "Step 6: Variant Calling with GATK HaplotypeCaller"
gatk HaplotypeCaller -R "${REFERENCE}" -I "${OUTDIR}/${SAMPLE}.recal.bam" \
      -O "${OUTDIR}/${SAMPLE}.raw.vcf"

echo "Step 7: Variant filtering"
gatk VariantFiltration -R "${REFERENCE}" -V "${OUTDIR}/${SAMPLE}.raw.vcf" \
      --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "BasicFilter" \
      -O "${OUTDIR}/${SAMPLE}.filtered.vcf"

echo "Step 8: Index filtered VCF"
bcftools index "${OUTDIR}/${SAMPLE}.filtered.vcf"

echo "Pipeline finished! Results in ${OUTDIR}"
