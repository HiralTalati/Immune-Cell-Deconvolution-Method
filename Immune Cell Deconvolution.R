# ============================================================
# 1. PACKAGES
# ============================================================

packages <- c(
  "tidyverse",
  "readr",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "ggpubr",
  "openxlsx"
)

installed <- rownames(installed.packages())

for (p in packages) {
  if (!p %in% installed) {
    install.packages(p)
  }
}

library(tidyverse)
library(readr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(openxlsx)


# ============================================================
# 2. TEN HUB GENES
# ============================================================

hub_genes <- c(
  "NDUFS5",
  "NDUFA1",
  "NDUFB3",
  "COX7A2",
  "ATP5ME",
  "UQCRH",
  "UQCRHL",
  "TOMM7",
  "ATP6V1E1",
  "SLIRP"
)


# ============================================================
# 3. READ IMMUCELLAI OUTPUT
# ============================================================

immucellai <- read.xlsx(
  "C:/Users/tahi3002/OneDrive - NIQ/Desktop/New folder/Bioinformatics/for manuscript/Cell deconvolution/ImmuneCellAbundance_sample.xlsx",
  check.names = FALSE
)

# Rename sample column
colnames(immucellai)[colnames(immucellai) == "sample"] <- "SampleId"

# Check
head(immucellai)
dim(immucellai)
colnames(immucellai)


# ============================================================
# 4. READ HUB-GENE EXPRESSION MATRIX
# ============================================================

expr <- read.xlsx(
  "C:/Users/tahi3002/OneDrive - NIQ/Desktop/New folder/Bioinformatics/for manuscript/Cell deconvolution/BCT_Expr_Matrix_Symbol_hubgenes.xlsx",
  check.names = FALSE
)

head(expr)
dim(expr)
colnames(expr)


# ============================================================
# 5. CONVERT HUB-GENE MATRIX
#    Original format:
#    rows = genes
#    columns = samples
#
#    Required format:
#    rows = samples
#    columns = genes
# ============================================================

# First column contains gene symbols
gene_names <- as.character(expr[[1]])

# Remaining columns contain expression values
expr_values <- expr[, -1, drop = FALSE]

# Convert expression values to numeric
expr_values <- as.data.frame(
  lapply(
    expr_values,
    function(x) as.numeric(as.character(x))
  ),
  check.names = FALSE
)

# Assign gene names as row names
rownames(expr_values) <- gene_names


# Transpose
hub_expr <- as.data.frame(
  t(expr_values),
  check.names = FALSE
)

# Original sample IDs are now row names
hub_expr$SampleId <- rownames(hub_expr)

# Move SampleId to first column
hub_expr <- hub_expr[
  ,
  c("SampleId", setdiff(colnames(hub_expr), "SampleId")),
  drop = FALSE
]


# ============================================================
# 6. CHECK HUB-GENE EXPRESSION MATRIX
# ============================================================

cat("Hub expression dimensions:\n")
print(dim(hub_expr))

cat("\nHub expression columns:\n")
print(colnames(hub_expr))

cat("\nAre all 10 hub genes present?\n")
print(hub_genes %in% colnames(hub_expr))

cat("\nNumber of hub genes present:\n")
print(sum(hub_genes %in% colnames(hub_expr)))

cat("\nAre hub-gene columns numeric?\n")
print(
  sapply(
    hub_expr[, hub_genes, drop = FALSE],
    is.numeric
  )
)


# ============================================================
# 7. CHECK SAMPLE IDs
# ============================================================

cat("\nFirst hub-expression SampleIDs:\n")
print(head(hub_expr$SampleId))

cat("\nFirst ImmuCellAI SampleIDs:\n")
print(head(immucellai$SampleId))

common_samples <- intersect(
  hub_expr$SampleId,
  immucellai$SampleId
)

cat("\nNumber of common samples:\n")
print(length(common_samples))


# ============================================================
# 8. MERGE THE TWO DATASETS
# ============================================================

immu_hub <- inner_join(
  hub_expr,
  immucellai,
  by = "SampleId"
)

cat("\nMerged dataset dimensions:\n")
print(dim(immu_hub))


# ============================================================
# 9. IDENTIFY IMMUCELLAI IMMUNE-CELL VARIABLES
# ============================================================

# Everything in ImmuCellAI except SampleId
immune_cells <- setdiff(
  colnames(immucellai),
  "SampleId"
)

cat("\nImmune-cell columns:\n")
print(immune_cells)

cat("\nNumber of immune-cell columns:\n")
print(length(immune_cells))


# ============================================================
# 10. CRITICAL DATA-TYPE CHECK
# ============================================================

cat("\nAre all hub genes numeric?\n")

hub_numeric_check <- sapply(
  immu_hub[, hub_genes, drop = FALSE],
  is.numeric
)

print(hub_numeric_check)


cat("\nAre all immune-cell variables numeric?\n")

immune_numeric_check <- sapply(
  immu_hub[, immune_cells, drop = FALSE],
  is.numeric
)

print(immune_numeric_check)


# Identify problematic columns, if any
problem_hub <- names(hub_numeric_check)[!hub_numeric_check]

problem_immune <- names(immune_numeric_check)[!immune_numeric_check]

cat("\nNon-numeric hub genes:\n")
print(problem_hub)

cat("\nNon-numeric immune-cell variables:\n")
print(problem_immune)


# ============================================================
# 11. STOP IF NON-NUMERIC VARIABLES ARE FOUND
# ============================================================

if (length(problem_hub) > 0 | length(problem_immune) > 0) {
  
  stop(
    "Non-numeric variables detected. Check problem_hub and problem_immune before continuing."
  )
}


# ============================================================
# 12. SPEARMAN CORRELATION
# ============================================================

results <- data.frame()

for (gene in hub_genes) {
  
  for (cell in immune_cells) {
    
    x <- immu_hub[[gene]]
    y <- immu_hub[[cell]]
    
    # Confirm numeric
    if (!is.numeric(x) || !is.numeric(y)) {
      next
    }
    
    # Complete paired observations
    valid <- complete.cases(x, y)
    
    n_valid <- sum(valid)
    
    # Require at least 10 paired samples
    if (n_valid >= 10) {
      
      test <- suppressWarnings(
        cor.test(
          x[valid],
          y[valid],
          method = "spearman",
          exact = FALSE
        )
      )
      
      results <- rbind(
        results,
        data.frame(
          Gene = gene,
          Immune_Cell = cell,
          Spearman_Rho = as.numeric(test$estimate),
          P_value = test$p.value,
          N = n_valid,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}


# ============================================================
# 13. MULTIPLE-TESTING CORRECTION
# ============================================================

results$FDR <- p.adjust(
  results$P_value,
  method = "BH"
)


# ============================================================
# 14. SORT RESULTS
# ============================================================

results <- results %>%
  arrange(
    FDR,
    desc(abs(Spearman_Rho))
  )


# ============================================================
# 15. INSPECT RESULTS
# ============================================================

cat("\nNumber of correlation tests:\n")
print(nrow(results))

cat("\nTop 20 associations:\n")
print(head(results, 20))


# ============================================================
# 16. SIGNIFICANT ASSOCIATIONS
# ============================================================

significant_results <- results %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)

cat("\nNumber of FDR-significant associations:\n")
print(nrow(significant_results))

print(significant_results)


# ============================================================
# 17. STRONGER SIGNIFICANT ASSOCIATIONS
# ============================================================

strong_results <- results %>%
  filter(
    FDR < 0.05,
    abs(Spearman_Rho) >= 0.30
  ) %>%
  arrange(
    FDR,
    desc(abs(Spearman_Rho))
  )

cat("\nSignificant associations with |rho| >= 0.30:\n")
print(strong_results)


# ============================================================
# 18. SAVE RESULTS
# ============================================================

write.csv(
  results,
  "Hub_gene_ImmuneCellAI_Spearman_all_results.csv",
  row.names = FALSE
)

write.csv(
  significant_results,
  "Hub_gene_ImmuneCellAI_Spearman_FDR_significant.csv",
  row.names = FALSE
)

write.csv(
  strong_results,
  "Hub_gene_ImmuneCellAI_Spearman_strong_significant.csv",
  row.names = FALSE
)

# ============================================================
# HUB GENE × IMMUNE-CELL SPEARMAN CORRELATION HEATMAP

# ============================================================

# ------------------------------------------------------------
# 1. Install/load required packages
# ------------------------------------------------------------

packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer"
)

installed <- rownames(installed.packages())

for (p in packages) {
  if (!p %in% installed) {
    install.packages(p)
  }
}

library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)


# ------------------------------------------------------------
# 2. Read complete correlation results
# ------------------------------------------------------------

results <- read.csv(
  "Hub_gene_ImmuneCellAI_Spearman_all_results.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

head(results)
dim(results)
colnames(results)


# ------------------------------------------------------------
# 3. Make sure numerical columns are numeric
# ------------------------------------------------------------

results$Spearman_Rho <- as.numeric(results$Spearman_Rho)
results$FDR <- as.numeric(results$FDR)


# ------------------------------------------------------------
# 4. Create correlation matrix
#
# IMPORTANT:
# ALL correlations are retained:
# positive, negative and nonsignificant.
# ------------------------------------------------------------

cor_matrix <- results %>%
  select(Gene, Immune_Cell, Spearman_Rho) %>%
  pivot_wider(
    names_from = Immune_Cell,
    values_from = Spearman_Rho
  ) %>%
  as.data.frame()

rownames(cor_matrix) <- cor_matrix$Gene
cor_matrix$Gene <- NULL

cor_matrix <- as.matrix(cor_matrix)

mode(cor_matrix) <- "numeric"


# ------------------------------------------------------------
# 5. Create FDR matrix
# ------------------------------------------------------------

fdr_matrix <- results %>%
  select(Gene, Immune_Cell, FDR) %>%
  pivot_wider(
    names_from = Immune_Cell,
    values_from = FDR
  ) %>%
  as.data.frame()

rownames(fdr_matrix) <- fdr_matrix$Gene
fdr_matrix$Gene <- NULL

fdr_matrix <- as.matrix(fdr_matrix)

mode(fdr_matrix) <- "numeric"


# ------------------------------------------------------------
# 6. Arrange genes and immune cells alphabetically
#
# This reproduces the ordering of the heatmap generated above.
# ------------------------------------------------------------

cor_matrix <- cor_matrix[
  order(rownames(cor_matrix)),
  order(colnames(cor_matrix)),
  drop = FALSE
]

fdr_matrix <- fdr_matrix[
  rownames(cor_matrix),
  colnames(cor_matrix),
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Create significance annotation matrix
#
# *     FDR < 0.05
# **    FDR < 0.01
# ***   FDR < 0.001
#
# Nonsignificant correlations remain in the heatmap.
# ------------------------------------------------------------

sig_matrix <- matrix(
  "",
  nrow = nrow(fdr_matrix),
  ncol = ncol(fdr_matrix),
  dimnames = dimnames(fdr_matrix)
)

sig_matrix[fdr_matrix < 0.05] <- "*"
sig_matrix[fdr_matrix < 0.01] <- "**"
sig_matrix[fdr_matrix < 0.001] <- "***"


# ------------------------------------------------------------
# 8. Define correlation color scale
#
# Fixed range: -0.25 to +0.25
# Center = 0
# ------------------------------------------------------------

heat_colors <- colorRampPalette(
  c(
    "navy",
    "white",
    "firebrick3"
  )
)(100)

heat_breaks <- seq(
  -0.25,
  0.25,
  length.out = 101
)


# ------------------------------------------------------------
# 9. Generate heatmap
# ------------------------------------------------------------

pheatmap(
  cor_matrix,
  
  # Significance stars
  display_numbers = sig_matrix,
  number_color = "black",
  fontsize_number = 8,
  
  # Color scale
  color = heat_colors,
  breaks = heat_breaks,
  
  # IMPORTANT:
  # No clustering — same structure as displayed heatmap
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  # Font sizes
  fontsize_row = 10,
  fontsize_col = 8,
  
  # Rotate immune-cell labels
  angle_col = 45,
  
  # Cell borders
  border_color = "lightgray",
  
  # Legend
  legend = TRUE,
  
  # Legend title
  legend_labels = c(
    "-0.25",
    "-0.20",
    "-0.15",
    "-0.10",
    "-0.05",
    "0",
    "0.05",
    "0.10",
    "0.15",
    "0.20",
    "0.25"
  ),
  
  # Figure title
  main = "Hub Gene–Immune Cell Association Analysis",
  
  # Cell dimensions
  cellwidth = 22,
  cellheight = 30,
  
  # Significance text
  fontsize = 10
)


# ------------------------------------------------------------
# 10. Save high-resolution PNG
# ------------------------------------------------------------

png(
  filename = "Hub_gene_ImmuneCellAI_Spearman_heatmap.png",
  width = 16,
  height = 6.5,
  units = "in",
  res = 600
)

pheatmap(
  cor_matrix,
  
  display_numbers = sig_matrix,
  number_color = "black",
  fontsize_number = 8,
  
  color = heat_colors,
  breaks = heat_breaks,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  fontsize_row = 10,
  fontsize_col = 8,
  
  angle_col = 45,
  
  border_color = "lightgray",
  
  legend = TRUE,
  
  main = "Hub Gene–Immune Cell Association Analysis",
  
  cellwidth = 22,
  cellheight = 30,
  
  fontsize = 10
)

dev.off()


# ------------------------------------------------------------
# 11. Save PDF version as well
# ------------------------------------------------------------

pdf(
  file = "Hub_gene_ImmuneCellAI_Spearman_heatmap.pdf",
  width = 16,
  height = 6.5
)

pheatmap(
  cor_matrix,
  
  display_numbers = sig_matrix,
  number_color = "black",
  fontsize_number = 8,
  
  color = heat_colors,
  breaks = heat_breaks,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  fontsize_row = 10,
  fontsize_col = 8,
  
  angle_col = 45,
  
  border_color = "lightgray",
  
  legend = TRUE,
  
  main = "Hub Gene–Immune Cell Association Analysis",
  
  cellwidth = 22,
  cellheight = 30,
  
  fontsize = 10
)

dev.off()


# ------------------------------------------------------------
# 12. Check final dimensions
# ------------------------------------------------------------

cat(
  "\nFinal heatmap dimensions:",
  nrow(cor_matrix),
  "hub genes ×",
  ncol(cor_matrix),
  "immune-cell populations\n"
)

cat(
  "Total correlations represented:",
  sum(!is.na(cor_matrix)),
  "\n"
)

cat(
  "FDR-significant correlations:",
  sum(fdr_matrix < 0.05, na.rm = TRUE),
  "\n"
)
