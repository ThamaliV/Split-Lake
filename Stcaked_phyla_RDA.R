library(tidyverse)
library(dplyr)
library(vegan)
library(gplots)
library(dendextend)
library(Rtsne)
library(vegan)
library(fastcluster)
library(tidyverse)
library(scales)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggsci)
taxa_table <- read.table("E:/Thami-Uni/R/TCNData analysis/taxa_table_small.tsv", 
                         sep = "\t", 
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         row.names = NULL)  # Explicitly set row.names to NULL

asv_table <- read.table("E:/Thami-Uni/R/TCNData analysis/asv_table_small.tsv", 
                        sep = "\t", 
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        check.names = FALSE,
                        row.names = NULL)  # Prevent R from modifying column names
location_table <- read.csv("E:/Thami-Uni/R/TCNData analysis/location.csv", header=TRUE, na.strings = c("", "NA"))
location_table <- na.omit(location_table)
location_table <- location_table %>%
  rename(station_ID = Station_ID) %>%
  mutate(station_ID = as.character(station_ID)) %>%
  select(sample_ID, station_ID, location_name)

meta_table <- read.csv("E:/Thami-Uni/R/TCNData analysis/TCN_sampling_metadataR.csv", header=TRUE, na.strings = c("", "NA"))
meta_table <- na.omit (meta_table)

colnames(taxa_table)[1] <- "DNA_Sequence"
colnames(asv_table)[1] <- "Sample_ID"

asv_long <- asv_table %>% 
  pivot_longer(
    cols = - Sample_ID,
    names_to = "DNA_Sequence", 
    values_to= "Proportional_Abundance")


merged_table <- taxa_table %>%
  left_join(asv_long, by = "DNA_Sequence")

asv_long <- asv_table %>% 
  pivot_longer(
    cols = - Sample_ID,
    names_to = "DNA_Sequence", 
    values_to= "Proportional_Abundance")
joined_df <- inner_join(location_table, merged_table, by = c("sample_ID" = "Sample_ID"))


joined_df1 <- joined_df %>%
  group_by(station_ID) %>%
  mutate(Proportional_Abundance = Proportional_Abundance / sum(Proportional_Abundance)) %>%
  ungroup()
station_order <- c("YLC", "AK", "OR", "ALe","ALw","NB", "BR", "OFCN","SCN", "CL", "SLTP", "MB", "KD", "HRB", "GR", "JSI")
joined_df1$station_ID <- factor(joined_df1$station_ID, levels = station_order)




my_colors <- c(
  "Acidobacteriota" = "#1f78b4",       # blue
  "Actinomycetota" = "#6a3d9a",        # purple
  "Armatimonadota" = "#ff7f00",        # orange
  "Bacillota" = "#e31a1c",             # red
  "Bacteroidota" = "#b15928",          # brown
  "Balneolota" = "#000080",            # navy blue (changed)
  "Bdellovibrionota" = "#cab2d6",      # lavender
  "Campylobacterota" = "#ff6600",      # bright orange (changed)
  "Candidatus Kapabacteria" = "#fb9a99", # light red
  "Chloroflexota" = "#b8860b",         # dark yellow (changed)
  "Cyanobacteriota" = "#006400",       # dark green
  "Gemmatimonadota" = "#bc80bd",       # pink purple
  "Hydrogenedentes" = "#33a02c",       # dark green
  "MBNT15" = "#ffed6f",                # bright yellow
  "Myxococcota" = "#fb8072",           # coral
  "Nitrospirota" = "#8dd3c7",          # aqua
  "Patescibacteria" = "#80b1d3",       # denim blue
  "Planctomycetota" = "#fccde5",       # very light pink
  "Pseudomonadota" = "#008080",        # teal blue (changed)
  "Spirochaetota" = "#bebada",         # lilac
  "Thermodesulfobacteriota" = "#d9d9d9", # gray
  "Thermoproteota" = "#fbfbbf",        # butter yellow
  "Verrucomicrobiota" = "#7f0000"      # maroon
)


jpeg("stacked_phyla2.jpeg", width = 1200, height = 800, res = 150)

stacked <- ggplot(joined_df1, aes(x = station_ID, y = Proportional_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  #scale_fill_d3(palette = "category20")+
 
  scale_fill_manual(values = my_colors) +
  labs(x = "Station", y = "Proportional_Abundance")+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",           # or "bottom" if too wide
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    plot.margin = margin(1, 1, 1, 1, "cm")  # top, right, bottom, left
  )


print(stacked)


dev.off()

joined_df %>%
  group_by(station_ID) %>%
  summarise(total = sum(Proportional_Abundance)) %>%
  print()

unique(joined_df1$station_ID)

......................................................................................
#get the number odf unique phyla, family and genera in data

length(unique(joined_df1$Phylum))

length(unique(joined_df1$Family))

length(unique(joined_df1$Genus))

length(unique(joined_df1$Class))

length(unique(joined_df1$Order))


joined_df1%>%
  group_by(Station_ID) %>%
  slice_max(order_by = Proportional_Abundance)

library(dplyr)
library(tidyverse)
library(corrplot)
library(tidyr)

total_all <- sum(joined_df1$Proportional_Abundance)

percentage_pseudo <- joined_df1 %>%
  filter(Phylum == "Pseudomonadota") %>%
  summarise(Percentage = sum(Proportional_Abundance) / total_all * 100)

print(percentage_pseudo)

meta_table <- read.csv("E:/Thami-Uni/R/TCNData analysis/TCN_Sampling_MetadataR.csv")

merged_station_meta <- left_join(joined_df1, meta_table, by = "station_ID")

colnames(joined_df1)
colnames(meta_table)




#RDA plot

library(ggfortify)

phylum_summary <- joined_df1 %>%
  group_by(station_ID, Phylum) %>%
  summarise(Proportional_Abundance = sum(Proportional_Abundance), .groups = "drop")

phylum_wide <- phylum_summary %>%
  pivot_wider(names_from = Phylum, values_from = Proportional_Abundance, values_fill = 0)

combined_meta <- left_join(meta_table, phylum_wide, by = "station_ID")

# Define environmental variables
env_vars <- c("PH", "Dissolved_Oxygen", "Turbidity_NTU", "Sechchi_depth_m", "Temperature")
env_data <- combined_meta[, env_vars]

# Define phyla columns — all numeric proportional abundance columns
phyla_vars_1 <- setdiff(names(combined_meta), c("Date", "station_Name", "station_ID", "Longitude", "Latitude", env_vars))
phyla_data <- combined_meta[, phyla_vars_1]

# Remove phyla with all zero values (to avoid RDA errors)
phyla_data <- phyla_data[, colSums(phyla_data) > 0]



rda_model <- rda(phyla_data ~ ., data = env_data)
summary(rda_model)

jpeg("rda.jpeg", width = 1600, height = 900, res = 150)
par(mar = c(4, 5, 2, 2))


rda <- plot(rda_model, scaling = 2, type = "n")

# Add only points for samples (sites) - no labels
points(rda_model, display = "sites", pch = 19, col = "blue", scaling = 2)

# Add only points for species (phyla) - no labels
points(rda_model, display = "species", pch = 17, col = "red", scaling = 2)

# Add env variable arrows and labels
env_arrows <- scores(rda_model, display = "bp", scaling = 2)
arrows(0, 0, env_arrows[,1], env_arrows[,2], length = 0.1, col = "darkgreen", lwd = 2)
text(env_arrows[,1], env_arrows[,2], labels = rownames(env_arrows), col = "darkgreen", cex = 0.8, pos = 3)


legend("topright", legend = c("Stations", "Phyla", "Env variables"),
       pch = c(19, 17, NA),
       lty = c(NA, NA, 1),
       col = c("blue", "red", "darkgreen"),
       pt.cex = 1.5,
       bty = "n",
       lwd = 2,
       merge = TRUE,
       cex = 0.7)  # smaller text size)


print(rda)


dev.off()

