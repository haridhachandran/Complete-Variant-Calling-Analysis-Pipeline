# Load libraries
library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(GenomicRanges)

# Set file paths
vcf_file <- "results/sample1/sample1.filtered.vcf"

# Read VCF
vcf <- readVcf(vcf_file, "hg19")

# Summary statistics
cat("Number of variants:", length(vcf), "\n")

# Extract variant info
variants <- rowRanges(vcf)
info <- info(vcf)

# Create data frame of key variant metrics
var_df <- data.frame(
  CHROM = seqnames(variants),
  POS = start(variants),
  REF = as.character(ref(vcf)),
  ALT = as.character(unlist(alt(vcf))),
  QUAL = qual(vcf),
  FILTER = fixed(vcf)$FILTER,
  DP = info$DP,
  MQ = info$MQ
)

# Basic QC plot: Quality distribution
ggplot(var_df, aes(x = QUAL)) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Variant Quality Score Distribution", x = "Quality Score", y = "Count")

# Allele frequency plot (if available as AF field; else derived from DP/AD)
if ("AF" %in% colnames(info(vcf))) {
  af <- info(vcf)$AF
  af_df <- data.frame(AF = unlist(af))
  ggplot(af_df, aes(x = AF)) +
    geom_histogram(binwidth = 0.05, fill = "tomato", color = "black") +
    theme_minimal() +
    labs(title = "Allele Frequency Distribution", x = "Allele Frequency", y = "Count")
}

# Variant type counts (SNP vs Indel)
variant_types <- ifelse(nchar(ref(vcf)) == 1 & nchar(unlist(alt(vcf))) == 1, "SNP", "Indel")
table(variant_types)

# Plot SNP/Indel counts
barplot(table(variant_types),
        col = c("lightgreen", "salmon"),
        main = "Variant Types",
        ylab = "Count")

# Save plots
ggsave("variant_quality_histogram.png")
ggsave("allele_frequency_histogram.png")

# Export variant summary table
write.csv(var_df, "variant_summary.csv", row.names = FALSE)
