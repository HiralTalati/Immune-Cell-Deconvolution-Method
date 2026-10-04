# ============================================================
# IMMUCELLAI IMMUNE-CELL ABUNDANCE
# CTL vs MCI vs AD
# ALL 1,301 SAMPLES INCLUDED
# Borderline samples recoded as MCI
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

packages <- c(
  "openxlsx",
  "tidyverse",
  "rstatix",
  "ggplot2",
  "ggpubr",
  "pheatmap"
)

installed <- rownames(installed.packages())

for (p in packages) {
  
  if (!p %in% installed) {
    install.packages(p)
  }
}

library(openxlsx)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(pheatmap)


# ------------------------------------------------------------
# 2. File paths
# ------------------------------------------------------------

immucellai_file <- 
  "C:/Users/tahi3002/OneDrive - NIQ/Desktop/New folder/Bioinformatics/for manuscript/Cell deconvolution/ImmuneCellAbundance_sample.xlsx"

metadata_file <- 
  "C:/Users/tahi3002/OneDrive - NIQ/Desktop/New folder/Bioinformatics/for manuscript/Cell deconvolution/metadata.csv"


# ------------------------------------------------------------
# 3. Read ImmuCellAI output
# ------------------------------------------------------------

immucellai <- read.xlsx(
  immucellai_file,
  check.names = FALSE
)

cat(
  "\nImmuCellAI dimensions:",
  nrow(immucellai),
  "samples x",
  ncol(immucellai),
  "columns\n"
)


# ------------------------------------------------------------
# 4. Read metadata
# ------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Metadata dimensions:",
  nrow(metadata),
  "samples x",
  ncol(metadata),
  "columns\n"
)


# ------------------------------------------------------------
# 5. Rename sample column
# ------------------------------------------------------------

if ("sample" %in% colnames(immucellai)) {
  
  colnames(immucellai)[
    colnames(immucellai) == "sample"
  ] <- "SampleID"
}


# ------------------------------------------------------------
# 6. Check sample matching
# ------------------------------------------------------------

common_samples <- intersect(
  immucellai$SampleID,
  metadata$SampleID
)

cat(
  "\nNumber of common samples:",
  length(common_samples),
  "\n"
)

if (length(common_samples) == 0) {
  
  stop(
    "No matching SampleIDs found between ImmuCellAI and metadata."
  )
}


# ------------------------------------------------------------
# 7. Merge datasets
# ------------------------------------------------------------

immune_data <- inner_join(
  immucellai,
  metadata,
  by = "SampleID"
)


cat(
  "\nMerged dataset:",
  nrow(immune_data),
  "samples\n"
)


# ------------------------------------------------------------
# 8. Check original group distribution
# ------------------------------------------------------------

cat(
  "\nOriginal clinical group distribution:\n"
)

print(
  table(
    immune_data$Group,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 9. IMPORTANT:
#    Recode borderline samples as MCI
# ------------------------------------------------------------

immune_data <- immune_data %>%
  
  mutate(
    
    Group = case_when(
      
      Group %in% c(
        "borderline",
        "Borderline",
        "BORDERLINE"
      ) ~ "MCI",
      
      TRUE ~ as.character(Group)
    )
  )


# ------------------------------------------------------------
# 10. Keep CTL, MCI and AD
# ------------------------------------------------------------

immune_data <- immune_data %>%
  
  filter(
    Group %in% c(
      "CTL",
      "MCI",
      "AD"
    )
  )


# ------------------------------------------------------------
# 11. Set factor order
# ------------------------------------------------------------

immune_data$Group <- factor(
  immune_data$Group,
  levels = c(
    "CTL",
    "MCI",
    "AD"
  )
)


# ------------------------------------------------------------
# 12. Final group distribution
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "FINAL CLINICAL GROUP DISTRIBUTION\n",
  "========================================\n"
)

print(
  table(
    immune_data$Group
  )
)

cat(
  "\nTotal samples analyzed:",
  nrow(immune_data),
  "\n"
)


# ------------------------------------------------------------
# 13. Identify immune-cell columns
# ------------------------------------------------------------

metadata_columns <- c(
  "SampleID",
  "Age",
  "Gender",
  "Group",
  "GSEID"
)

immune_cells <- setdiff(
  colnames(immucellai),
  "SampleID"
)

cat(
  "\nNumber of immune-cell populations:",
  length(immune_cells),
  "\n"
)

print(immune_cells)


# ------------------------------------------------------------
# 14. Convert immune-cell columns to numeric
# ------------------------------------------------------------

immune_data[immune_cells] <- lapply(
  
  immune_data[immune_cells],
  
  function(x) {
    
    suppressWarnings(
      as.numeric(
        as.character(x)
      )
    )
  }
)


# ------------------------------------------------------------
# 15. Kruskal-Wallis test
# ------------------------------------------------------------

kw_results <- data.frame()


for (cell in immune_cells) {
  
  temp <- immune_data %>%
    
    select(
      Group,
      all_of(cell)
    ) %>%
    
    rename(
      Abundance = all_of(cell)
    ) %>%
    
    filter(
      !is.na(Abundance),
      is.finite(Abundance)
    )
  
  
  if (
    length(unique(temp$Group)) == 3 &&
    nrow(temp) >= 10
  ) {
    
    test <- kruskal.test(
      Abundance ~ Group,
      data = temp
    )
    
    
    kw_results <- rbind(
      
      kw_results,
      
      data.frame(
        
        Immune_Cell = cell,
        
        Kruskal_Wallis_H =
          as.numeric(
            test$statistic
          ),
        
        P_value =
          test$p.value,
        
        N =
          nrow(temp),
        
        stringsAsFactors = FALSE
      )
    )
  }
}


# ------------------------------------------------------------
# 16. BH-FDR correction
# ------------------------------------------------------------

kw_results$FDR <- p.adjust(
  kw_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 17. Sort by FDR
# ------------------------------------------------------------

kw_results <- kw_results %>%
  
  arrange(
    FDR
  )


# ------------------------------------------------------------
# 18. Add significance label
# ------------------------------------------------------------

kw_results <- kw_results %>%
  
  mutate(
    
    Significance = case_when(
      
      FDR < 0.001 ~ "***",
      
      FDR < 0.01 ~ "**",
      
      FDR < 0.05 ~ "*",
      
      TRUE ~ "NS"
    )
  )


# ------------------------------------------------------------
# 19. Display results
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "KRUSKAL-WALLIS RESULTS\n",
  "========================================\n"
)

print(
  kw_results
)


# ------------------------------------------------------------
# 20. Significant immune cells
# ------------------------------------------------------------

significant_cells <- kw_results %>%
  
  filter(
    FDR < 0.05
  )


cat(
  "\n========================================\n",
  "FDR-SIGNIFICANT IMMUNE CELLS\n",
  "========================================\n"
)

print(
  significant_cells
)

cat(
  "\nNumber of FDR-significant immune cells:",
  nrow(significant_cells),
  "\n"
)


# ------------------------------------------------------------
# 21. Save complete KW results
# ------------------------------------------------------------

write.xlsx(
  
  kw_results,
  
  "ImmuCellAI_CTL_MCI_AD_Kruskal_Wallis_all_results.xlsx",
  
  overwrite = TRUE
)


# ------------------------------------------------------------
# 22. Save significant results
# ------------------------------------------------------------

write.xlsx(
  
  significant_cells,
  
  "ImmuCellAI_CTL_MCI_AD_Kruskal_Wallis_significant.xlsx",
  
  overwrite = TRUE
)


# ------------------------------------------------------------
# 23. Dunn post-hoc analysis
# ------------------------------------------------------------

dunn_results <- data.frame()


if (
  nrow(significant_cells) > 0
) {
  
  for (cell in significant_cells$Immune_Cell) {
    
    temp <- immune_data %>%
      
      select(
        Group,
        all_of(cell)
      ) %>%
      
      rename(
        Abundance = all_of(cell)
      ) %>%
      
      filter(
        !is.na(Abundance),
        is.finite(Abundance)
      )
    
    
    dunn <- temp %>%
      
      dunn_test(
        
        Abundance ~ Group,
        
        p.adjust.method = "BH"
      ) %>%
      
      mutate(
        Immune_Cell = cell
      )
    
    
    dunn_results <- rbind(
      
      dunn_results,
      
      dunn
    )
  }
}


# ------------------------------------------------------------
# 24. Save Dunn results
# ------------------------------------------------------------

write.xlsx(
  
  dunn_results,
  
  "ImmuCellAI_CTL_MCI_AD_Dunn_posthoc_results.xlsx",
  
  overwrite = TRUE
)


# ------------------------------------------------------------
# 25. Group-wise median and IQR
# ------------------------------------------------------------

median_data <- immune_data %>%
  
  select(
    Group,
    all_of(immune_cells)
  ) %>%
  
  pivot_longer(
    
    cols = all_of(immune_cells),
    
    names_to = "Immune_Cell",
    
    values_to = "Abundance"
  ) %>%
  
  group_by(
    Immune_Cell,
    Group
  ) %>%
  
  summarise(
    
    Median =
      median(
        Abundance,
        na.rm = TRUE
      ),
    
    Q1 =
      quantile(
        Abundance,
        0.25,
        na.rm = TRUE
      ),
    
    Q3 =
      quantile(
        Abundance,
        0.75,
        na.rm = TRUE
      ),
    
    IQR =
      IQR(
        Abundance,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 26. Save group-wise summary
# ------------------------------------------------------------

write.xlsx(
  
  median_data,
  
  "ImmuCellAI_CTL_MCI_AD_groupwise_summary.xlsx",
  
  overwrite = TRUE
)


# ------------------------------------------------------------
# 27. Identify top nominally changing cells
# ------------------------------------------------------------

top_nominal <- kw_results %>%
  
  arrange(
    P_value
  ) %>%
  
  slice_head(
    n = 6
  )


cat(
  "\n========================================\n",
  "TOP NOMINAL IMMUNE-CELL DIFFERENCES\n",
  "========================================\n"
)

print(
  top_nominal
)


# ------------------------------------------------------------
# 28. Save top nominal results
# ------------------------------------------------------------

write.xlsx(
  
  top_nominal,
  
  "ImmuCellAI_top_nominal_immune_cell_changes.xlsx",
  
  overwrite = TRUE
)


# ------------------------------------------------------------
# 29. Prepare exploratory visualization
#
# IMPORTANT:
# These are NOT FDR-significant cells.
# They are displayed only as the strongest nominal
# differences and must be interpreted accordingly.
# ------------------------------------------------------------

plot_cells <- top_nominal$Immune_Cell


plot_data <- immune_data %>%
  
  select(
    SampleID,
    Group,
    all_of(plot_cells)
  ) %>%
  
  pivot_longer(
    
    cols = all_of(plot_cells),
    
    names_to = "Immune_Cell",
    
    values_to = "Abundance"
  )


plot_data$Immune_Cell <- factor(
  
  plot_data$Immune_Cell,
  
  levels = plot_cells
)


# ------------------------------------------------------------
# 30. Exploratory violin + boxplot
# ------------------------------------------------------------

immune_violin <- ggplot(
  
  plot_data,
  
  aes(
    x = Group,
    y = Abundance,
    fill = Group
  )
) +
  
  geom_violin(
    
    trim = FALSE,
    
    alpha = 0.65,
    
    scale = "width"
  ) +
  
  geom_boxplot(
    
    width = 0.15,
    
    outlier.shape = NA,
    
    fill = "white",
    
    alpha = 0.85
  ) +
  
  facet_wrap(
    
    ~ Immune_Cell,
    
    scales = "free_y",
    
    ncol = 3
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "CTL" = "#4C78A8",
      
      "MCI" = "#F2A541",
      
      "AD" = "#D95F59"
    )
  ) +
  
  labs(
    
    title =
      "Immune-cell abundance across CTL, MCI and AD",
    
    x = NULL,
    
    y =
      "Estimated immune-cell abundance"
  ) +
  
  theme_classic(
    
    base_size = 12
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 10
      ),
    
    legend.position = "top",
    
    legend.title =
      element_blank()
  )


print(
  immune_violin
)


# ------------------------------------------------------------
# 31. Save exploratory plot
# ------------------------------------------------------------

ggsave(
  
  "ImmuCellAI_top_nominal_immune_cell_changes.png",
  
  plot = immune_violin,
  
  width = 12,
  
  height = 9,
  
  units = "in",
  
  dpi = 600
)


ggsave(
  
  "ImmuCellAI_top_nominal_immune_cell_changes.pdf",
  
  plot = immune_violin,
  
  width = 12,
  
  height = 9,
  
  units = "in"
)


# ------------------------------------------------------------
# 32. Final summary
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "IMMUCELLAI DISEASE-STAGE ANALYSIS COMPLETED\n",
  "============================================\n",
  "\nCTL samples: 487",
  "\nMCI samples: 326",
  "\nAD samples: 488",
  "\nTotal samples: 1301",
  "\nImmune-cell populations tested: ",
  length(immune_cells),
  "\nFDR-significant populations: ",
  nrow(significant_cells),
  "\n============================================\n"
)