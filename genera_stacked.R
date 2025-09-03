library(tidyverse)
library(ggplot2)
library(viridis)
library(RColorBrewer)

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

# Aggregate by genus and station to get total abundance per genus per station
genus_station_summary <- joined_df1 %>%
  group_by(station_ID, Genus) %>%
  summarise(Total_Abundance = sum(Proportional_Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Genus) & Genus != "")

# Get the top 20 most abundant genera across all stations
top_20_genera <- genus_station_summary %>%
  group_by(Genus) %>%
  summarise(Overall_Abundance = sum(Total_Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Overall_Abundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)


# Create the final dataset for plotting
plot_data <- genus_station_summary %>%
  mutate(
    Genus_Category = if_else(Genus %in% top_20_genera, Genus, "Other")
  ) %>%
  group_by(station_ID, Genus_Category) %>%
  summarise(Total_Abundance = sum(Total_Abundance, na.rm = TRUE), .groups = "drop")

# Set the order for stations
station_order <- c("YLC", "AK", "OR", "ALe", "ALw", "NB", "BR", "OFCN", "SCN", "CL", "SLTP", "MB", "KD", "HRB", "GR", "JSI")

# Convert station_ID to factor with specified order
plot_data$station_ID <- factor(plot_data$station_ID, levels = station_order)

# Define colors for the top 20 genera 
genus_colors <- c(
  "#6a3d9a",   # purple
  "#ff7f00",   # orange
  "#e31a1c",   # red
  "#b15928",   # brown
  "#000080",   # navy blue
  "#cab2d6",   # lavender
  "#ff6600",   # bright orange
  "#fb9a99",   # light red
  "#b8860b",   # dark yellow
  "#006400",   # dark green
  "#bc80bd",   # pink purple
  "#33a02c",   # dark green
  "#ffed6f",   # bright yellow
  "#fb8072",   # coral
  "#8dd3c7",   # aqua
  "#fccde5",   # very light pink
  "#008080",   # teal blue
  "#bebada",   # lilac
  "#d9d9d9",   # gray
  "#7f0000"    # maroon
)

# Shuffle the color palette
set.seed(123)  
genus_colors <- sample(genus_colors)

# Assign colors directly to the top 20 genera, with special handling for specific genera

genus_color_mapping <- setNames(genus_colors, top_20_genera)

# Assign specific colors for certain genera
for (genus in top_20_genera) {
  if (grepl("cyanobium", genus, ignore.case = TRUE)) {
    genus_color_mapping[genus] <- "#4B0082"  # dark purple for cyanobium
    cat("Cyanobium found:", genus, "-> Dark Purple\n")
  } else if (grepl("hgcl", genus, ignore.case = TRUE)) {
    genus_color_mapping[genus] <- "#4169E1"  # royal blue for hgcl clade
    cat("HGCL Clade found:", genus, "-> Royal Blue\n")
  } else if (grepl("lacisediminimonas", genus, ignore.case = TRUE)) {
    genus_color_mapping[genus] <- "#90EE90"  # light green for lacisediminimonas
    cat("Lacisediminimonas found:", genus, "-> Light Green\n")
  } else if (grepl("incertae sedis", genus, ignore.case = TRUE)) {
    genus_color_mapping[genus] <- "#808080"  # gray for incertae sedis
    cat("Incertae sedis found:", genus, "-> Gray\n")
  }
}

# Add "Other" as black
genus_color_mapping["Other"] <- "#000000"  # black for Other

color_palette <- genus_color_mapping


# Create the stacked bar plot
stacked_plot <- ggplot(plot_data, aes(x = station_ID, y = Total_Abundance, fill = Genus_Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Top 20 Most Abundant Bacterial Genera by Station",
    subtitle = "Other genera are grouped together",
    x = "Station ID",
    y = "Proportional Abundance",
    fill = "Genus"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.8))

print(stacked_plot)

# Save the plot
ggsave("stacked_bar_plot_top20_genera.pdf", plot = stacked_plot, 
       width = 14, height = 8, dpi = 300)

ggsave("stacked_bar_plot_top20_genera.png", plot = stacked_plot, 
       width = 14, height = 8, dpi = 300)


# Calculate percentage of total abundance represented by top 20 genera
total_abundance <- sum(genus_station_summary$Total_Abundance, na.rm = TRUE)
top20_abundance <- sum(genus_station_summary$Total_Abundance[genus_station_summary$Genus %in% top_20_genera], na.rm = TRUE)
percentage_top20 <- (top20_abundance / total_abundance) * 100

