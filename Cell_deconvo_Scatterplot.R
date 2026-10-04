# ============================================================
# HUB GENE–IMMUNE CELL ASSOCIATION
# REPRESENTATIVE SPEARMAN SCATTER PLOTS
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

packages <- c(
  "tidyverse",
  "ggplot2",
  "ggpubr"
)

installed <- rownames(installed.packages())

for (p in packages) {
  if (!p %in% installed) {
    install.packages(p)
  }
}

library(tidyverse)
library(ggplot2)
library(ggpubr)


# ------------------------------------------------------------
# 2. Check required objects
# ------------------------------------------------------------

if (!exists("immu_hub")) {
  stop(
    "Object 'immu_hub' does not exist. ",
    "Run the previous data preparation script first."
  )
}

if (!exists("results")) {
  stop(
    "Object 'results' does not exist. ",
    "Run the Spearman correlation script first."
  )
}


# ------------------------------------------------------------
# 3. Basic checks
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DATA CHECK\n")
cat("============================================\n")

cat(
  "immu_hub: ",
  nrow(immu_hub),
  " samples x ",
  ncol(immu_hub),
  " columns\n",
  sep = ""
)

cat(
  "results: ",
  nrow(results),
  " associations\n",
  sep = ""
)

cat("\n")


# ------------------------------------------------------------
# 4. Selected representative associations
# ------------------------------------------------------------

plot_pairs <- data.frame(
  
  Gene = c(
    "ATP6V1E1",
    "COX7A2",
    "SLIRP",
    "NDUFB3",
    "UQCRH",
    "NDUFS5"
  ),
  
  Immune_Cell = c(
    "basophils",
    "cDC1",
    "mast_cell",
    "neutrophils",
    "GC_B",
    "cDC1"
  ),
  
  Panel = c(
    "A",
    "B",
    "C",
    "D",
    "E",
    "F"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Check column names
# ------------------------------------------------------------

required_columns <- unique(
  c(
    plot_pairs$Gene,
    plot_pairs$Immune_Cell
  )
)

missing_columns <- setdiff(
  required_columns,
  colnames(immu_hub)
)


if (length(missing_columns) > 0) {
  
  cat("\nMissing columns:\n")
  print(missing_columns)
  
  stop(
    "\nSome gene/immune-cell columns were not found ",
    "in 'immu_hub'."
  )
}


# ------------------------------------------------------------
# 6. Check FDR column
# ------------------------------------------------------------

if (!"FDR" %in% colnames(results)) {
  
  stop(
    "\nThe 'results' object does not contain an FDR column.\n",
    "Please check your correlation results table."
  )
}


# ------------------------------------------------------------
# 7. Function to create scatter plot
# ------------------------------------------------------------

make_scatter <- function(
    data,
    gene,
    immune_cell,
    panel_label
) {
  
  
  # ----------------------------------------------------------
  # Extract required variables
  # ----------------------------------------------------------
  
  plot_data <- data %>%
    
    select(
      all_of(gene),
      all_of(immune_cell)
    ) %>%
    
    rename(
      Gene_Expression = all_of(gene),
      Immune_Abundance = all_of(immune_cell)
    )
  
  
  # ----------------------------------------------------------
  # Convert to numeric safely
  # ----------------------------------------------------------
  
  plot_data$Gene_Expression <- suppressWarnings(
    as.numeric(
      as.character(
        plot_data$Gene_Expression
      )
    )
  )
  
  
  plot_data$Immune_Abundance <- suppressWarnings(
    as.numeric(
      as.character(
        plot_data$Immune_Abundance
      )
    )
  )
  
  
  # ----------------------------------------------------------
  # Remove missing/non-finite observations
  # ----------------------------------------------------------
  
  plot_data <- plot_data %>%
    
    filter(
      is.finite(Gene_Expression),
      is.finite(Immune_Abundance)
    )
  
  
  # ----------------------------------------------------------
  # Check number of observations
  # ----------------------------------------------------------
  
  n <- nrow(plot_data)
  
  if (n < 10) {
    
    stop(
      paste0(
        "\nInsufficient observations for ",
        gene,
        " vs ",
        immune_cell,
        ". N = ",
        n
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Spearman correlation
  # ----------------------------------------------------------
  
  test <- suppressWarnings(
    cor.test(
      plot_data$Gene_Expression,
      plot_data$Immune_Abundance,
      method = "spearman",
      exact = FALSE
    )
  )
  
  
  rho <- as.numeric(
    test$estimate
  )
  
  p_value <- test$p.value
  
  
  # ----------------------------------------------------------
  # Retrieve FDR from results
  # ----------------------------------------------------------
  
  fdr_match <- results %>%
    
    filter(
      Gene == gene,
      Immune_Cell == immune_cell
    )
  
  
  if (nrow(fdr_match) > 0) {
    
    fdr <- as.numeric(
      fdr_match$FDR[1]
    )
    
  } else {
    
    fdr <- NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Format FDR
  # ----------------------------------------------------------
  
  if (!is.na(fdr)) {
    
    fdr_text <- format.pval(
      fdr,
      digits = 3,
      eps = 0.001
    )
    
  } else {
    
    fdr_text <- "NA"
  }
  
  
  # ----------------------------------------------------------
  # Statistical annotation
  # ----------------------------------------------------------
  
  annotation_text <- paste0(
    
    "Spearman \u03c1 = ",
    sprintf(
      "%.3f",
      rho
    ),
    
    "\nFDR = ",
    fdr_text,
    
    "\nN = ",
    format(
      n,
      big.mark = ","
    )
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    
    plot_data,
    
    aes(
      x = Gene_Expression,
      y = Immune_Abundance
    )
    
  ) +
    
    
    # Sample points
    geom_point(
      
      color = "#3B82A0",
      
      size = 1.8,
      
      alpha = 0.55
    ) +
    
    
    # LOESS trend
    geom_smooth(
      
      method = "loess",
      
      formula = y ~ x,
      
      se = TRUE,
      
      color = "#D95F5F",
      
      fill = "#F2B6B6",
      
      linewidth = 0.9,
      
      alpha = 0.25
    ) +
    
    
    # Statistical annotation
    annotate(
      
      "text",
      
      x = Inf,
      
      y = Inf,
      
      label = annotation_text,
      
      hjust = 1.05,
      
      vjust = 1.15,
      
      size = 3.8,
      
      color = "black"
    ) +
    
    
    # Labels
    labs(
      
      title = paste0(
        panel_label,
        ". ",
        gene,
        " vs ",
        immune_cell
      ),
      
      x = paste0(
        gene,
        " expression"
      ),
      
      y = paste0(
        immune_cell,
        " abundance"
      )
    ) +
    
    
    # Theme
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      
      axis.title = element_text(
        size = 11
      ),
      
      axis.text = element_text(
        size = 10
      ),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    )
  
  
  return(p)
}


# ------------------------------------------------------------
# 8. Generate all six plots
# ------------------------------------------------------------

cat("\n============================================\n")
cat("GENERATING SCATTER PLOTS\n")
cat("============================================\n\n")


scatter_plots <- list()


for (i in seq_len(nrow(plot_pairs))) {
  
  cat(
    "Generating Panel ",
    plot_pairs$Panel[i],
    ": ",
    plot_pairs$Gene[i],
    " vs ",
    plot_pairs$Immune_Cell[i],
    "\n",
    sep = ""
  )
  
  
  scatter_plots[[i]] <- make_scatter(
    
    data = immu_hub,
    
    gene = plot_pairs$Gene[i],
    
    immune_cell = plot_pairs$Immune_Cell[i],
    
    panel_label = plot_pairs$Panel[i]
  )
}


# ------------------------------------------------------------
# 9. Display individual plots
# ------------------------------------------------------------

print(scatter_plots[[1]])

print(scatter_plots[[2]])

print(scatter_plots[[3]])

print(scatter_plots[[4]])

print(scatter_plots[[5]])

print(scatter_plots[[6]])


# ------------------------------------------------------------
# 10. Create combined 2 × 3 figure
# ------------------------------------------------------------

combined_scatter <- ggarrange(
  
  plotlist = scatter_plots,
  
  ncol = 2,
  
  nrow = 3
)


# ------------------------------------------------------------
# 11. Display combined figure
# ------------------------------------------------------------

print(combined_scatter)


# ------------------------------------------------------------
# 12. Save combined PNG
# ------------------------------------------------------------

ggsave(
  
  filename =
    "Hub_gene_ImmuneCellAI_scatterplots.png",
  
  plot =
    combined_scatter,
  
  width = 12,
  
  height = 15,
  
  units = "in",
  
  dpi = 600,
  
  bg = "white"
)


# ------------------------------------------------------------
# 13. Save combined PDF
# ------------------------------------------------------------

ggsave(
  
  filename =
    "Hub_gene_ImmuneCellAI_scatterplots.pdf",
  
  plot =
    combined_scatter,
  
  width = 12,
  
  height = 15,
  
  units = "in",
  
  bg = "white"
)


# ------------------------------------------------------------
# 14. Save individual plots
# ------------------------------------------------------------

for (i in seq_along(scatter_plots)) {
  
  filename <- paste0(
    
    "Scatter_",
    
    plot_pairs$Gene[i],
    
    "_",
    
    plot_pairs$Immune_Cell[i],
    
    ".png"
  )
  
  
  ggsave(
    
    filename = filename,
    
    plot = scatter_plots[[i]],
    
    width = 6,
    
    height = 5,
    
    units = "in",
    
    dpi = 600,
    
    bg = "white"
  )
}


# ------------------------------------------------------------
# 15. Completion message
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SIX SCATTER PLOTS GENERATED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "\nCombined files:\n",
  "1. Hub_gene_ImmuneCellAI_scatterplots.png\n",
  "2. Hub_gene_ImmuneCellAI_scatterplots.pdf\n",
  sep = ""
)

cat(
  "\nIndividual PNG files were also saved.\n"
)