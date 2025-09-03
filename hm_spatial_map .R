#heatmap code 

getwd()
library(tidyverse)
library(dplyr)
library(vegan)
library(gplots)
library(dendextend)
library(Rtsne)
library(vegan)
library(fastcluster)
library(ggplot2)
library(tidyverse)
library(scales)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(pheatmap)

getwd()
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



print(head(asv_table))

print("Taxa table columns:")
print(colnames(taxa_table))
print("ASV table columns:")
print(colnames(asv_table))

head(asv_table)

# Ensure the DNA sequence column names match and are of character type
colnames(taxa_table)[1] <- "DNA_Sequence"
colnames(asv_table)[1] <- "Sample_ID"

# Convert DNA_Sequence to character in both tables
taxa_table$DNA_Sequence <- as.character(taxa_table$DNA_Sequence)
asv_table$DNA_Sequence <- as.character(asv_table$DNA_Sequence)

str(asv_table)
asv_long <- asv_table %>% 
            pivot_longer(
            cols = - Sample_ID,
            names_to = "DNA_Sequence", 
            values_to= "Proportional_Abundance")


merged_table <- taxa_table %>%
  left_join(asv_long, by = "DNA_Sequence")

head(merged_table)

sapply(asv_table, class)  # to check if the variables are numeric



..............................................................................................................................................................

#merged_taxa_table <- merge(asv_table, taxa_table, by = "DNA_Sequence")

#merged_taxa_table$DNA_Sequence <- NULL

# View the merged table

#getting genus names on the HM


print(merged_table)

data_melted <- melt(merged_table, 
                    id.vars = c("Sample_ID", "ASV"), 
                    measure.vars = "Proportional_Abundance",
                    variable.name = "Metric", 
                    value.name = "Proportional_Abundance")

data_with_genus <- merge(data_melted, taxa_table, by = "ASV")

genus_abundance <- data_with_genus %>%
  group_by(Sample_ID, Genus) %>%
  summarize(Proportional_Abundance = sum(Proportional_Abundance), .groups = "drop")

# Remove rows with NA in Genus
genus_abundance_clean <- genus_abundance %>%
  filter(!is.na(Genus))

genus_matrix <- genus_abundance_clean %>%
  pivot_wider(names_from = Sample_ID, values_from = Proportional_Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Genus") %>%  # Make genus names rownames
  as.matrix() 

str(genus_matrix)


.............................................................................................................................
#Plotting Split lake data in locationwise

joined_df <- inner_join(location_table, merged_table, by = c("sample_ID" = "Sample_ID"))

#converting joined df as a heat mappable dataframe

heatmap_data <- joined_df %>%
  pivot_wider(names_from = station_ID, values_from = Proportional_Abundance, values_fill = 0)

# Save Genus as merged_table# Save Genus as rownames
heatmap_mat <- as.data.frame(heatmap_data)

rownames(heatmap_mat) <- heatmap_mat$Genus

heatmap_mat <- heatmap_mat %>% select(-Genus)  # drop Genus column before converting to matrix

data_melted_location <- melt(joined_df, 
                    id.vars = c("station_ID", "ASV"), 
                    measure.vars = "Proportional_Abundance",
                    variable.name = "Metric", 
                    value.name = "Proportional_Abundance")

data_with_locations <- merge(data_melted_location, taxa_table, by = "ASV")

genus_abundance <- data_with_locations %>%
  group_by(station_ID, Genus) %>%
  summarize(Proportional_Abundance = sum(Proportional_Abundance), .groups = "drop")

asv_abundance <- data_with_locations %>%
  group_by(station_ID, ASV) %>%
  summarize(Proportional_Abundance = sum(Proportional_Abundance), .groups = "drop")

# Remove rows with NA in Genus
genus_abundance_location <- genus_abundance %>%
  filter(!is.na(Genus))

# Remove rows with NA in Genus
asv_abundance_location <- asv_abundance %>%
  filter(!is.na(ASV))

locationD_matrix <- genus_abundance_location %>%
  pivot_wider(names_from = station_ID, values_from = Proportional_Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Genus") %>%  # Make genus names rownames
  as.matrix()

locationD_matrix_asv <- asv_abundance_location %>%
  pivot_wider(names_from = station_ID, values_from = Proportional_Abundance, values_fill = 0) %>%
  column_to_rownames(var = "ASV") %>%  # Make genus names rownames
  as.matrix()

locationD_matrix_asv = locationD_matrix_asv[rowSums(locationD_matrix_asv > 0) > 2, colSums(locationD_matrix_asv > 0) > 2] #Removing NA values

pdf(file = "hm_SplitLk1.pdf", width = 12, height = 12)

hm_SplitLk <- heatmap.2(t(locationD_matrix_asv)^.25,              # The data (transformed numeric data)
                      trace = "none",                # No trace lines in the heatmap
                      col = viridis(100),            # Use the viridis color palette
                      Rowv = TRUE,                   # Clustering rows (samples)
                      Colv = TRUE,                   # Clustering columns (ASVs/DNA sequences)
                      distfun = function(x) vegdist(x, method="bray"),
                      hclustfun = function(x) hclust(x, method="ward.D"),
                      dendrogram = "both",           # Show dendrograms for rows and columns
                      main = "Split Lake Only",   # Add a title
                      cexRow = 1,                  # Adjust row label size
                      cexCol = 0.4                   # Adjust column label size
                      
)

nclus=6

clusters = cutree(hm_SplitLk$rowDendrogram, k = nclus)

hm_SplitLk <- heatmap.2(t(locationD_matrix_asv)^.25,              # The data (transformed numeric data)
                        trace = "none",                # No trace lines in the heatmap
                        col = viridis(100),            # Use the viridis color palette
                        Rowv = TRUE,                   # Clustering rows (samples)
                        Colv = TRUE,                   # Clustering columns (ASVs/DNA sequences)
                        distfun = function(x) vegdist(x, method="bray"),
                        hclustfun = function(x) hclust(x, method="ward.D"),
                        dendrogram = "both",           # Show dendrograms for rows and columns
                        main = "Split Lake Only",   # Add a title
                        RowSideColors = rainbow(nclus)[as.numeric(as.factor(clusters))],
                        cexRow = 1,                  # Adjust row label size
                        cexCol = 0.4                   # Adjust column label size
                        
)


#par(mar = c(20, 20, 5, 10))

dev.off()


......................................................................................................................



png("hm_SplitLk.png", width = 1800, height = 1400, res = 150)


hm_SplitLk <- heatmap.2(t(locationD_matrix_asv)^.25,              # The data (transformed numeric data)
                        trace = "none",                # No trace lines in the heatmap
                        col = viridis(100),            # Use the viridis color palette
                        Rowv = TRUE,                   # Clustering rows (samples)
                        Colv = TRUE,                   # Clustering columns (ASVs/DNA sequences)
                        distfun = function(x) vegdist(x, method="bray"),
                        hclustfun = function(x) hclust(x, method="ward.D"),
                        dendrogram = "both",           # Show dendrograms for rows and columns
                        main = "Split Lake Only",   # Add a title
                        cexRow = 1,                  # Adjust row label size
                        cexCol = 0.4                   # Adjust column label size
                        
)

nclus=6

clusters = cutree(hm_SplitLk$rowDendrogram, k = nclus)

hm_SplitLk <- heatmap.2(t(locationD_matrix_asv)^.25,              # The data (transformed numeric data)
                        trace = "none",                # No trace lines in the heatmap
                        col = viridis(100),            # Use the viridis color palette
                        Rowv = TRUE,                   # Clustering rows (samples)
                        Colv = TRUE,                   # Clustering columns (ASVs/DNA sequences)
                        distfun = function(x) vegdist(x, method="bray"),
                        hclustfun = function(x) hclust(x, method="ward.D"),
                        dendrogram = "both",           # Show dendrograms for rows and columns
                        main = "Split Lake Only",   # Add a title
                        RowSideColors = rainbow(nclus)[as.numeric(as.factor(clusters))],
                        cexRow = 1,                  # Adjust row label size
                        cexCol = 0.4                   # Adjust column label size
                        
)


#par(mar = c(20, 20, 5, 10))

dev.off()

#mapping the clusters on geological locations

library(ggmap)

clean_meta <- read.csv("E:/Thami-Uni/R/TCNData analysis/clean_meta.csv", header=TRUE)


register_stadiamaps("8911d93b-8768-4156-8a20-72bbb76fe601") #API key
bbox <- make_bbox(lon=clean_meta$Longitude, lat = clean_meta$Latitude, f=0.1)
mapz <- get_stadiamap(bbox, maptype = "stamen_terrain", zoom=12)


clean_meta$clusters = as.factor(clusters[clean_meta$Station_ID])
cluster_colors = setNames(rainbow(nclus),1:nclus)

map_plot <- ggmap(mapz)+
    geom_point(data= clean_meta, aes(x= Longitude, y= Latitude, color=clusters), size=5, alpha=0.6)+
    scale_color_manual(values= cluster_colors)
    
  
#    geom_text(data= clean_meta, aes(x= Longitude, y= Latitude, label= Station_ID), size=5)
map_plot
#dev.off()
#dev.set(which = dev.next())


...........................................................................................


png("map.png", width = 1800, height = 1400, res = 150)

register_stadiamaps("8911d93b-8768-4156-8a20-72bbb76fe601") #API key
bbox <- make_bbox(lon=clean_meta$Longitude, lat = clean_meta$Latitude, f=0.1)
mapz <- get_stadiamap(bbox, maptype = "stamen_terrain", zoom=10)


clean_meta$clusters = as.factor(clusters[clean_meta$Station_ID])
cluster_colors = setNames(rainbow(nclus),1:nclus)

map_plot <- ggmap(mapz)+
  geom_point(data= clean_meta, aes(x= Longitude, y= Latitude, color=clusters), size=5, alpha=0.6)+
  scale_color_manual(values= cluster_colors)

dev.off()
#    geom_text(data= clean_meta, aes(x= Longitude, y= Latitude, label= Station_ID), size=5)
map_plot

#Cyanobacteriales
#Synechococcales
choose_taxa = unique(taxa_table[taxa_table$Phylum=="Cyanobacteriota","Family"])
choose_taxa = choose_taxa[!is.na(choose_taxa)]

filterd_orders <- joined_with_clusters %>% filter(Phylum=="Cyanobacteriota") %>% filter(Family %in% choose_taxa)

grouped_orders <- filterd_orders %>%
    group_by(Family, location_name, Cluster) %>%
    #summarise(sum_value = sum("Proportional_Abundance"))
summarise(sum_value = sum(Proportional_Abundance, na.rm = TRUE))


#location_orders <- grouped_orders%>%
                  #inner_join(meta_table, by="Station_ID")



location_Orders <- grouped_orders %>%
  inner_join(meta_table, by = c("location_name" = "Station_Name"))


cluster_colors2 <- cluster_colors
names(cluster_colors2)= paste0("Cluster_", 1:6)

map_plot2 <- ggmap(mapz)+
  geom_point(data= location_Orders, aes(x= Longitude, y= Latitude, color=Cluster, size= sum_value), alpha=0.6)+
  facet_wrap(~Family, ncol=1)+
  scale_color_manual(values= cluster_colors2)+
  scale_size(range=c(0,20))

ggsave(map_plot2,filename="map_cyano_family.pdf", width=8, height=24)
print(map_plot2)

