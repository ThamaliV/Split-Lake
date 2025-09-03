library(tidyverse)
library(dplyr)
library(corrplot)
library(pheatmap)
library(ggplot2)

# Read the data files
taxa_table <- read.table("taxa_table_small.tsv", 
                         sep = "\t", 
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         row.names = NULL)

asv_table <- read.table("asv_table_small.tsv", 
                        sep = "\t", 
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        check.names = FALSE,
                        row.names = NULL)

location_table <- read.csv("location.csv", header=TRUE, na.strings = c("", "NA"))
location_table <- na.omit(location_table)
location_table <- location_table %>%
  rename(station_ID = Station_ID) %>%
  mutate(station_ID = as.character(station_ID)) %>%
  select(sample_ID, station_ID, location_name)

meta_table <- read.csv("TCN_Sampling_MetadataR.csv", header=TRUE, na.strings = c("", "NA"))
meta_table <- na.omit(meta_table)

# Clean column names
colnames(taxa_table)[1] <- "DNA_Sequence"
colnames(asv_table)[1] <- "Sample_ID"

# Create long format for ASV data
asv_long <- asv_table %>% 
  pivot_longer(
    cols = -Sample_ID,
    names_to = "DNA_Sequence", 
    values_to = "Proportional_Abundance")

# Merge taxa and ASV data
merged_table <- taxa_table %>%
  left_join(asv_long, by = "DNA_Sequence")

# Join with location data
joined_df <- inner_join(location_table, merged_table, by = c("sample_ID" = "Sample_ID"))

# Normalize abundances by station
joined_df1 <- joined_df %>%
  group_by(station_ID) %>%
  mutate(Proportional_Abundance = Proportional_Abundance / sum(Proportional_Abundance)) %>%
  ungroup()

# Aggregate bacterial data by phylum and station
phylum_summary <- joined_df1 %>%
  group_by(station_ID, Phylum) %>%
  summarise(Proportional_Abundance = sum(Proportional_Abundance), .groups = "drop")

# Convert to wide format for correlation analysis
phylum_wide <- phylum_summary %>%
  pivot_wider(names_from = Phylum, values_from = Proportional_Abundance, values_fill = 0)

# Merge with environmental metadata
combined_meta <- left_join(meta_table, phylum_wide, by = "station_ID")

# Get ALL environmental variables (excluding non-numeric columns)
env_vars <- names(meta_table %>% select(where(is.numeric)))

# Get phyla variables (exclude environmental variables and other metadata)
phyla_vars <- setdiff(names(combined_meta %>% select(where(is.numeric))), 
                      c(env_vars, "Longitude", "Latitude"))

# Remove any phyla with zero variance
phyla_vars <- phyla_vars[sapply(combined_meta[phyla_vars], function(x) var(x, na.rm = TRUE)) > 0]

# Calculate Spearman correlations with significance testing
cor_results <- data.frame()

for (phy in phyla_vars) {
  for (env in env_vars) {
    # Create pairwise clean data
    temp_df <- combined_meta[, c(phy, env)] %>% 
      filter(!is.na(.data[[phy]]) & !is.na(.data[[env]]))
    
    # Check if we have enough data points
    if (nrow(temp_df) >= 3 && 
        is.numeric(temp_df[[phy]]) && is.numeric(temp_df[[env]]) &&
        length(unique(temp_df[[phy]])) > 1 && length(unique(temp_df[[env]])) > 1) {
      
      # Calculate correlation
      cor_test <- cor.test(temp_df[[phy]], temp_df[[env]], method = "spearman", exact = FALSE)
      
      cor_results <- rbind(cor_results, data.frame(
        Phylum = phy,
        Variable = env,
        Spearman_r = cor_test$estimate,
        p_value = cor_test$p.value
      ))
    }
  }
}

# Create correlation matrix for heatmap
cor_matrix <- cor_results %>%
  select(Phylum, Variable, Spearman_r) %>%
  pivot_wider(names_from = Variable, values_from = Spearman_r) %>%
  column_to_rownames("Phylum") %>%
  as.matrix()

# Create significance matrix
sig_matrix <- cor_results %>%
  select(Phylum, Variable, p_value) %>%
  pivot_wider(names_from = Variable, values_from = p_value) %>%
  column_to_rownames("Phylum") %>%
  as.matrix()

# Create significance symbols
sig_symbols <- matrix("", nrow = nrow(sig_matrix), ncol = ncol(sig_matrix))
sig_symbols[sig_matrix < 0.001] <- "***"
sig_symbols[sig_matrix < 0.01 & sig_matrix >= 0.001] <- "**"
sig_symbols[sig_matrix < 0.05 & sig_matrix >= 0.01] <- "*"

# Clean up phylum names for better display
clean_phylum_names <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b([a-z])([a-z]+)", "\\U\\1\\L\\2", x, perl = TRUE)
  return(x)
}

rownames(cor_matrix) <- clean_phylum_names(rownames(cor_matrix))

# Clean up variable names for better display
clean_var_names <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b([a-z])([a-z]+)", "\\U\\1\\L\\2", x, perl = TRUE)
  return(x)
}

colnames(cor_matrix) <- clean_var_names(colnames(cor_matrix))


# Using ggplot2 
cor_results_plot <- cor_results %>%
  mutate(
    Phylum_clean = clean_phylum_names(Phylum),
    Variable_clean = clean_var_names(Variable),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Create ggplot version - SQUARE PLOT
p <- ggplot(cor_results_plot, aes(x = Variable_clean, y = Phylum_clean, fill = Spearman_r)) +
  geom_tile() +
  geom_text(aes(label = significance), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#000080",  # Navy Blue
    mid = "#FFFFFF",  # White
    high = "#FF8C00",
    midpoint = 0,
    limits = c(-0.8, 0.8),
    name = "Spearman r"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Spearman Correlation: Bacterial Phyla vs Environmental Variables") +
  coord_fixed(ratio = 1)  # Make it perfectly square

# Save ggplot version - SQUARE PLOT
ggsave("spearman_correlation_heatmap_ggplot_square.pdf", plot = p, width = 12, height = 12, dpi = 300)

# Display the plot in RStudio
print(p)

# Print summary statistics
cat("\nCorrelation Summary:\n")
cat("Number of significant correlations (p < 0.05):", sum(cor_results$p_value < 0.05), "\n")
cat("Number of highly significant correlations (p < 0.01):", sum(cor_results$p_value < 0.01), "\n")
cat("Number of very highly significant correlations (p < 0.001):", sum(cor_results$p_value < 0.001), "\n")

# Show top correlations
cat("\nTop 10 strongest correlations:\n")
cor_results %>%
  arrange(desc(abs(Spearman_r))) %>%
  head(10) %>%
  mutate(
    Phylum_clean = clean_phylum_names(Phylum),
    Variable_clean = clean_var_names(Variable)
  ) %>%
  select(Phylum_clean, Variable_clean, Spearman_r, p_value) %>%
  print()


# Create a simple correlation matrix for selected variables only
selected_vars <- c(env_vars, phyla_vars)
cor_selected <- combined_meta %>% 
  select(all_of(selected_vars)) %>%
  cor(method = "spearman", use = "pairwise.complete.obs")

# Save corrplot version
pdf("spearman_corrplot.pdf", width = 12, height = 10)
corrplot(cor_selected, 
         method = "color", 
         type = "upper", 
         tl.cex = 0.8,
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 0.6,
         col = colorRampPalette(c("#000080", "#FFFFFF", "#FF8C00"))(100))
dev.off()

# Save corrplot version - JPEG with SMALLER numbers
jpeg("spearman_corrplot.jpg", width = 1200, height = 1000, res = 150, quality = 95)
corrplot(cor_selected, 
         method = "color", 
         type = "upper", 
         tl.cex = 0.8,
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 0.4,  # Smaller font size for numbers
         number.font = 1,  # Plain font (not bold) for numbers
         col = colorRampPalette(c("#000080", "#FFFFFF", "#FF8C00"))(100))  # Navy Blue to White to Dark Orange
dev.off()


# Show top correlations
cat("\nTop 10 strongest correlations:\n")
cor_results %>%
  arrange(desc(abs(Spearman_r))) %>%
  head(10) %>%
  mutate(
    Phylum_clean = clean_phylum_names(Phylum),
    Variable_clean = clean_var_names(Variable)
  ) %>%
  select(Phylum_clean, Variable_clean, Spearman_r, p_value) %>%
  print()

# Method 3: Simple corrplot (alternative visualization) - JPEG
# Create a simple correlation matrix for selected variables only (EXCLUDING latitude/longitude)
selected_vars <- c(env_vars, phyla_vars)
cor_selected <- combined_meta %>% 
  select(all_of(selected_vars)) %>%
  cor(method = "spearman", use = "pairwise.complete.obs")



spearman_r <-  cor_test$estimate   # correlation coefficient (rho)
p_value <-  cor_test$p.value       # p-value

print(spearman_r)
print(p_value)

citation("ggplot2")
citation("vegan")
citation("ggfortify")

