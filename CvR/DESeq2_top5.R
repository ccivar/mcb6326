# DESeq2 analysis — Top 5 Up/Down Regulated Genes
# Outputs: pheatmaps + CSV files

# ── SETUP ─────────────────────────────────────────────────────────────────────
setwd("")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("pheatmap")

library(DESeq2)
library(pheatmap)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
countsdata <- read.delim("Brassica_resistant_counts.txt", row.names = "gene_id")
coldata    <- read.delim("sample_info.txt", row.names = 1)

# Validate alignment between counts matrix and sample info
all(colnames(countsdata) %in% rownames(coldata))
all(colnames(countsdata) == rownames(coldata))

# ── BUILD DESEQ2 DATASET ──────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = countsdata,
  colData   = coldata,
  design    = ~ condition)

# Set control as the reference level
dds$condition <- relevel(dds$condition, ref = "control")

# Pre-filter: keep only rows with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# ── RUN DESEQ2 ────────────────────────────────────────────────────────────────
dds <- DESeq(dds)
res <- results(dds)

# Significant DEGs (padj < 0.05)
res05 <- subset(res, padj < 0.05)

# ── NORMALIZED COUNTS ─────────────────────────────────────────────────────────
ntd             <- normTransform(dds)                          # log2(count+1), for heatmaps
NormalizedCount <- counts(estimateSizeFactors(dds), normalized = TRUE)

# Sample names per condition
control_samples     <- rownames(coldata[coldata$condition == "control",     ])
resistant_samples <- rownames(coldata[coldata$condition == "resistant", ])

# Annotation data frame for heatmap columns
df <- as.data.frame(colData(dds)[, c("cultivar", "condition")])

# ── TOP 5 UPREGULATED ─────────────────────────────────────────────────────────
res_up  <- subset(res, padj < 0.05 & log2FoldChange > 0)
top5_up <- as.data.frame(res_up[order(res_up$padj), ][1:5, ])

top5_up_counts           <- as.data.frame(NormalizedCount[rownames(top5_up), ])
top5_up$mean_control     <- rowMeans(top5_up_counts[, control_samples])
top5_up$mean_resistant <- rowMeans(top5_up_counts[, resistant_samples])
top5_up$regulation       <- "upregulated"
top5_up_final            <- cbind(top5_up, top5_up_counts)

# ── TOP 5 DOWNREGULATED ───────────────────────────────────────────────────────
res_down  <- subset(res, padj < 0.05 & log2FoldChange < 0)
top5_down <- as.data.frame(res_down[order(res_down$padj), ][1:5, ])

top5_down_counts           <- as.data.frame(NormalizedCount[rownames(top5_down), ])
top5_down$mean_control     <- rowMeans(top5_down_counts[, control_samples])
top5_down$mean_resistant <- rowMeans(top5_down_counts[, resistant_samples])
top5_down$regulation       <- "downregulated"
top5_down_final            <- cbind(top5_down, top5_down_counts)

# ── CSV OUTPUT ────────────────────────────────────────────────────────────────
write.csv(top5_up_final,   file = "top5_upregulated_norm_counts.csv",   row.names = TRUE)
write.csv(top5_down_final, file = "top5_downregulated_norm_counts.csv", row.names = TRUE)

top5_combined <- rbind(top5_up_final, top5_down_final)
write.csv(top5_combined,   file = "top5_up_and_down_norm_counts.csv",   row.names = TRUE)

# ── HEATMAPS ──────────────────────────────────────────────────────────────────
# Top 5 upregulated
pheatmap(assay(ntd)[rownames(top5_up_final), ],
         cluster_rows   = TRUE,
         show_rownames  = TRUE,
         cluster_cols   = TRUE,
         annotation_col = df,
         main           = "Control vs resistant — Top 5 Upregulated Genes (padj < 0.05)")

# Top 5 downregulated
pheatmap(assay(ntd)[rownames(top5_down_final), ],
         cluster_rows   = TRUE,
         show_rownames  = TRUE,
         cluster_cols   = TRUE,
         annotation_col = df,
         main           = "Control vs resistant — Top 5 Downregulated Genes (padj < 0.05)")

# Combined (10 genes total)
pheatmap(assay(ntd)[rownames(top5_combined), ],
         cluster_rows   = TRUE,
         show_rownames  = TRUE,
         cluster_cols   = TRUE,
         annotation_col = df,
         main           = "Top 5 Up & Downregulated Genes (padj < 0.05)")