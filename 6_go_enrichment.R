#!/usr/bin/env Rscript

# GO Enrichment Analysis for TES vs GFP DEGs
# Separates upregulated and downregulated genes for pathway analysis

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

# Set up directories
input_file <- "./results/05_deseq2/significant_genes_TES_vs_GFP.txt"
output_dir <- "./results/06_go_enrichment"
upregulated_dir <- file.path(output_dir, "upregulated")
downregulated_dir <- file.path(output_dir, "downregulated")
plots_dir <- file.path(output_dir, "plots")

# Create output directories if they don't exist
dir.create(upregulated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(downregulated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading DEGs data...\n")
degs <- read_tsv(input_file, show_col_types = FALSE)

# Convert ENSG IDs to gene symbols and Entrez IDs
cat("Converting gene IDs...\n")
# Strip version numbers from Ensembl IDs (e.g., ENSG00000060718.22 -> ENSG00000060718)
degs <- degs %>%
  mutate(gene_id_clean = str_remove(gene_id, "\\.\\d+$"))

degs <- degs %>%
  mutate(
    gene_symbol = mapIds(org.Hs.eg.db,
                         keys = gene_id_clean,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first"),
    entrez_id = mapIds(org.Hs.eg.db,
                       keys = gene_id_clean,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  )
# Remove genes with missing mappings
degs_filtered <- degs %>%
  filter(!is.na(gene_symbol), !is.na(entrez_id))

cat("Total genes after ID mapping:", nrow(degs_filtered), "\n")

# Separate upregulated and downregulated genes
# Upregulated: log2FoldChange > 0 (higher in TES)
# Downregulated: log2FoldChange < 0 (lower in TES)
upregulated_genes <- degs_filtered %>%
  filter(log2FoldChange > 0) %>%
  pull(entrez_id)

downregulated_genes <- degs_filtered %>%
  filter(log2FoldChange < 0) %>%
  pull(entrez_id)

cat("Upregulated genes in TES:", length(upregulated_genes), "\n")
cat("Downregulated genes in TES:", length(downregulated_genes), "\n")

# Get background genes (all genes in the analysis)
background_genes <- degs_filtered$entrez_id

# Function to perform GO enrichment analysis
perform_go_enrichment <- function(gene_list, background, direction, ontology = "BP") {
  cat(paste("Performing GO", ontology, "enrichment for", direction, "genes...\n"))
  
  ego <- enrichGO(
    gene = gene_list,
    universe = background,
    OrgDb = org.Hs.eg.db,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  return(ego)
}

# Function to save results and plots
save_results <- function(ego, direction, ontology) {
  # Save results table
  results_file <- file.path(output_dir, tolower(direction), 
                           paste0("GO_", ontology, "_", direction, "_results.csv"))
  write.csv(as.data.frame(ego), results_file, row.names = FALSE)
  
  # Create summary
  summary_file <- file.path(output_dir, tolower(direction),
                           paste0("GO_", ontology, "_", direction, "_summary.txt"))
  sink(summary_file)
  cat("GO", ontology, "Enrichment Summary -", direction, "Genes\n")
  cat("=" , rep("=", 50), "\n", sep = "")
  
  # Access gene universe size from the enrichResult object properly
  ego_df <- as.data.frame(ego)
  cat("Significant terms found:", nrow(ego_df), "\n")
  cat("P-value cutoff: 0.05\n")
  cat("Q-value cutoff: 0.2\n\n")
  
  if (nrow(ego_df) > 0) {
    cat("Top 10 enriched terms:\n")
    print(head(ego_df[order(ego_df$pvalue), ], 10))
  }
  sink()
  
  cat("Results saved to:", results_file, "\n")
  cat("Summary saved to:", summary_file, "\n")
}

# Function to create plots
create_plots <- function(ego, direction, ontology) {
  ego_df <- as.data.frame(ego)
  if (nrow(ego_df) == 0) {
    cat("No significant terms found for plotting\n")
    return()
  }
  
  # Dot plot
  p1 <- dotplot(ego, showCategory = 20, title = paste("GO", ontology, "Enrichment -", direction)) +
    theme_bw() +
    theme(text = element_text(size = 12))
  
  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_dotplot.pdf")),
         plot = p1, width = 12, height = 8)
  
  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_dotplot.png")),
         plot = p1, width = 12, height = 8, dpi = 300)
  
  # Enrichment map
  if (nrow(ego_df) >= 5) {
    p2 <- emapplot(ego, showCategory = 30, title = paste("GO", ontology, "Enrichment Map -", direction)) +
      theme_bw()
    
    ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_emap.pdf")),
           plot = p2, width = 10, height = 8)
    
    ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_emap.png")),
           plot = p2, width = 10, height = 8, dpi = 300)
  }
  
  # Bar plot
  p3 <- barplot(ego, showCategory = 20, title = paste("GO", ontology, "Enrichment -", direction)) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_barplot.pdf")),
         plot = p3, width = 12, height = 8)
  
  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_barplot.png")),
         plot = p3, width = 12, height = 8, dpi = 300)
  cat("Plots saved for", direction, ontology, "\n")
}

# Perform GO enrichment analysis for different ontologies
ontologies <- c("BP", "MF", "CC")  # Biological Process, Molecular Function, Cellular Component

for (direction in c("Upregulated", "Downregulated")) {
  if (direction == "Upregulated") {
    gene_list <- upregulated_genes
  } else {
    gene_list <- downregulated_genes
  }
  
  cat("\n=== Analyzing", direction, "Genes ===\n")
  
  for (ont in ontologies) {
    ego <- perform_go_enrichment(gene_list, background_genes, direction, ont)
    save_results(ego, direction, ont)
    create_plots(ego, direction, ont)
  }
}

# Create comparative summary
cat("\nCreating comparative summary...\n")

summary_data <- data.frame(
  Direction = character(),
  Ontology = character(),
  Significant_Terms = integer(),
  Genes_in_Terms = integer(),
  stringsAsFactors = FALSE
)

for (direction in c("Upregulated", "Downregulated")) {
  for (ont in ontologies) {
    results_file <- file.path(output_dir, tolower(direction),
                             paste0("GO_", ont, "_", direction, "_results.csv"))
    if (file.exists(results_file)) {
      results <- read.csv(results_file)
      if (nrow(results) > 0) {
        summary_data <- rbind(summary_data, data.frame(
          Direction = direction,
          Ontology = ont,
          Significant_Terms = nrow(results),
          Genes_in_Terms = sum(str_count(results$geneID, "/") + 1),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Save comparative summary
summary_file <- file.path(output_dir, "GO_enrichment_comparative_summary.csv")
write.csv(summary_data, summary_file, row.names = FALSE)

# Create comparative plot
if (nrow(summary_data) > 0) {
  p_comparison <- ggplot(summary_data, aes(x = Ontology, y = Significant_Terms, fill = Direction)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "GO Enrichment Comparison: Upregulated vs Downregulated Genes",
         x = "GO Ontology",
         y = "Number of Significant Terms") +
    theme_bw() +
    theme(text = element_text(size = 12),
          legend.position = "top")
  
  ggsave(filename = file.path(plots_dir, "GO_enrichment_comparison.pdf"),
         plot = p_comparison, width = 10, height = 6)
  
  ggsave(filename = file.path(plots_dir, "GO_enrichment_comparison.png"),
         plot = p_comparison, width = 10, height = 6, dpi = 300)
}

# Create gene lists for reference
upregulated_df <- degs_filtered %>%
  filter(log2FoldChange > 0) %>%
  select(gene_id, gene_symbol, entrez_id, log2FoldChange, padj) %>%
  arrange(padj)

downregulated_df <- degs_filtered %>%
  filter(log2FoldChange < 0) %>%
  select(gene_id, gene_symbol, entrez_id, log2FoldChange, padj) %>%
  arrange(padj)

write.csv(upregulated_df, file.path(upregulated_dir, "upregulated_genes_list.csv"), row.names = FALSE)
write.csv(downregulated_df, file.path(downregulated_dir, "downregulated_genes_list.csv"), row.names = FALSE)

cat("\n=== GO Enrichment Analysis Complete ===\n")
cat("Results saved in:", output_dir, "\n")
cat("Upregulated results:", upregulated_dir, "\n")
cat("Downregulated results:", downregulated_dir, "\n")
cat("Plots saved in:", plots_dir, "\n")
cat("Comparative summary:", summary_file, "\n")

# Print final statistics
cat("\nFinal Statistics:\n")
cat("- Total significant genes analyzed:", nrow(degs_filtered), "\n")
cat("- Upregulated genes:", length(upregulated_genes), "\n")
cat("- Downregulated genes:", length(downregulated_genes), "\n")

if (file.exists(summary_file)) {
  summary_stats <- read.csv(summary_file)
  cat("\nSignificant GO terms found:\n")
  print(summary_stats)
}

cat("\nAnalysis completed successfully!\n")