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

asv_long <- asv_table %>% 
            pivot_longer(
            cols = - Sample_ID,
            names_to = "DNA_Sequence", 
            values_to= "Proportional_Abundance")


merged_table <- taxa_table %>%
  left_join(asv_long, by = "DNA_Sequence")

head(merged_table)

sapply(asv_table, class)  # to check if the variables are numeric

pdf(file = "heatmap_output.pdf", width = 12, height = 12)

rownames(asv_table) <- asv_table$Sample_ID
num_data <- asv_table[ , sapply(asv_table, is.numeric)]


hm_thami <- heatmap.2(transformed,              # The data (transformed numeric data)
  trace = "none",                # No trace lines in the heatmap
  col = viridis(100),            # Use the viridis color palette
  Rowv = TRUE,                   # Clustering rows (samples)
  Colv = TRUE,                   # Clustering columns (ASVs/DNA sequences)
  dendrogram = "both",           # Show dendrograms for rows and columns
  main = "Your Heatmap Title",   # Add a title
  cexRow = 0.4,                  # Adjust row label size
  cexCol = 0.4                   # Adjust column label size
  
)

par(mar = c(20, 20, 5, 10))

dev.off()

pdf("heatmap_output2.pdf", 
    width = 12,  # width in inches
    height = 12, # height in inches
    pointsize = 12)  # base font size

# Reset graphics state with larger margins
par(mar = c(15, 15, 5, 5))  # bottom, left, top, right margins

# Create the heatmap with adjusted margins
heatmap.2(transformed^.25,  # Using power transformation like in TCN-plots.r
          trace = "none",
          col = viridis,
          distfun = function(x) vegdist(x, method = "bray"),
          hclustfun = function(x) hclust(x, method = "ward.D"),
          margins = c(15, 15),  # increased margins for labels
          cexRow = 0.7,
          cexCol = 0.5,
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          srtCol = 45,  # rotate column labels
          offsetCol = 0.5,  # offset for column labels
          offsetRow = 0.5,  # offset for row labels
          lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0)),  # adjust layout
          lhei = c(1.5, 4, 1),  # adjust height ratios
          lwid = c(1.5, 4, 1))  # adjust width ratios

# Close the PDF device
dev.off()

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

pdf(file = "heatmap_genus_level.pdf", width = 12, height = 12)


hm_thami <- heatmap.2(genus_matrix,              # The data (transformed numeric data)
                      trace = "none",                # No trace lines in the heatmap
                      col = viridis(100),            # Use the viridis color palette
                      Rowv = TRUE,                   # Clustering rows (samples)
                      Colv = TRUE,                   # Clustering columns (ASVs/DNA sequences)
                      dendrogram = "both",           # Show dendrograms for rows and columns
                      main = "Heat map Genus level",   # Add a title
                      cexRow = 0.4,                  # Adjust row label size
                      cexCol = 0.4                   # Adjust column label size
                      
)

#par(mar = c(20, 20, 5, 10))

dev.off()
..................................................................................................

#Getting most abundant 20 genera on HM



genus_means <- rowMeans(genus_matrix)

# Get names of top 30 most abundant genera
top20_genera <- names(sort(genus_means, decreasing = TRUE))[1:20]


genus_top20_matrix <- genus_matrix[top20_genera, ]



pdf(file = "heatmap_top20genus.pdf", width = 12, height = 12)


hm_thami <- heatmap.2(genus_top20_matrix,              # The data (transformed numeric data)
                      trace = "none",                # No trace lines in the heatmap
                      col = viridis(100),            # Use the viridis color palette
                      Rowv = TRUE,                   # Clustering rows (samples)
                      Colv = TRUE,                   # Clustering columns (ASVs/DNA sequences)
                      dendrogram = "both",           # Show dendrograms for rows and columns
                      main = "Top 20 Genus",   # Add a title
                      cexRow = 0.4,                  # Adjust row label size
                      cexCol = 0.4                   # Adjust column label size
                      
)

#par(mar = c(20, 20, 5, 10))

dev.off()


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
#mapping the clusters on geological locations

library(ggmap)


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
................................................................................................................................................

#plotting stacked bar plot for split lake data 

#Get total abundance per Genus

top20_genera <- joined_df %>%
  group_by(Genus) %>%
  summarise(Total_Abundance = sum(Proportional_Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

filtered_df <- joined_df %>%
  filter(Genus %in% top20_genera)

#identify genera to be colored in Grey
genus_level <- unique(filtered_df$Genus)
grey_genere <- genus_level[grepl("incertae sedis", genus_level, ignore.case=TRUE)]
non_grey_genera <- setdiff(genus_level,grey_genere)

#generate viridis colors for non-gret genera
viridis_colors <- viridis::viridis(length(non_grey_genera))
names(viridis_colors) <- non_grey_genera

#Assign grey color to incertae sedis
custom_colors <-  c(setNames(rep("grey", length("grey_genera")),grey_genere),viridis_colors)

#plot
pdf("stacked_20 genus_plot.pdf", width = 10, height = 6)


stacked <- ggplot(filtered_df, aes(x = station_ID, y = Proportional_Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors)+
  labs(x = "Station", y = "Proportional_Abundance", title = "Top 20 Genus Abundance")+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",           # or "bottom" if too wide
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    plot.margin = margin(1, 1, 1, 1, "cm")  # top, right, bottom, left
  )
   

print(stacked)


dev.off()
........................................................................................................
#creating a stacked bar plot with clusters - 1


#fixing the uppercase issue :)
colnames(cluster_info)[colnames(cluster_info) == "Station_ID"] <- "station_ID"


cluster_info <- data.frame(station_ID = names(clusters),
                           Cluster=paste0("Cluster_", clusters))

joined_with_clusters <- joined_df %>%
  inner_join(cluster_info, by = "station_ID")

pdf("stacked_cluster_plot.pdf", width = 10, height = 6)


Stacked_cluster <- ggplot(joined_with_clusters, aes(x = station_ID, y = Proportional_Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Cluster, scales = "free_x", nrow = 1) +
  scale_fill_viridis_d() +
  labs(x = "Station", y = "Proportional Abundance", title = "Cluster Abnundance across Stations") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size=8),
    axis.text.y = element_text(size = 8),
    axix.title= element_text(size=10),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),   # smaller legend keys
    legend.text = element_text(size = 6),  # smaller legend text (adjustments for the plot to make it within the margins)
    legend.title = element_text(size = 8),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.05,'cm'),
    plot.title = element_text(size=11, face="bold"),
    plot.margin = margin(1,1,1,1,"cm")
  )+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))


print(Stacked_cluster)
dev.off()



#creating a stacked bar plot with clusters - 2

asv_clusters <- cutree(hm_SplitLk$colDendrogram, k = 6)

# Create a data frame to map ASVs to clusters
asv_cluster <- data.frame(
  ASV = names(asv_clusters),
  Cluster = paste0("Cluster_", asv_clusters)
)

#melting ASV abundance matrix to long format

library(reshape2)




long_melted <- asv_abundance %>%
  inner_join(asv_cluster, by= "ASV")



cluster_abundance_per_station <- long_melted %>%
  group_by(station_ID, Cluster) %>%
  summarise(Total_Abundance = sum(Proportional_Abundance, na.rm = TRUE)) %>%
  ungroup()

pdf("stacked_cluster_plot2.pdf", width = 10, height = 6)


stacked2 <- ggplot(cluster_abundance_per_station, aes(x = station_ID, y = Total_Abundance, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(x = "Station", y = "Total Abundance", title = "Cluster Abundance per Station") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(stacked2)
dev.off()









...........................................................................................


library(ggmap)


register_stadiamaps("8911d93b-8768-4156-8a20-72bbb76fe601") #API key
bbox <- make_bbox(lon=clean_meta$Longitude, lat = clean_meta$Latitude, f=0.1)
mapz <- get_stadiamap(bbox, maptype = "stamen_terrain", zoom=10)


clean_meta$clusters = as.factor(clusters[clean_meta$Station_ID])
cluster_colors = setNames(rainbow(nclus),1:nclus)

map_plot <- ggmap(mapz)+
  geom_point(data= clean_meta, aes(x= Longitude, y= Latitude, color=clusters), size=5, alpha=0.6)+
  scale_color_manual(values= cluster_colors)


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











#data_melted <- melt(merged_table, id.vars = c("Sample_ID", "ASV"), value.name = "Proportional_Abundance")

# Pivot the data to create a matrix with Sample_IDs as rows and ASVs as columns
data_pivot <- dcast(data_melted, Sample_ID ~ ASV, value.var = "Proportional_Abundance", fun.aggregate = sum)

# Create the heatmap
ggplot(data_melted, aes(x = ASV, y = Sample_ID, fill = Proportional_Abundance)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
  theme_minimal() +
  labs(title = "Heatmap of Proportional Abundance per ASV", x = "ASV", y = "Sample_ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




dist_mat <- dist(t(locationD_matrix_asv))          # transpose if clustering samples (columns)
clust <- hclust(dist_mat)
ordered_samples <- clust$labels[clust$order]



joined_df$station_ID <- factor(joined_df$station_ID, levels = ordered_samples)



print(ordered_samples)
length(ordered_samples)
head(joined_df$station_ID)



pdf("stacked_20 genus_plot_order1.pdf", width = 10, height = 6)

# Generate 20 distinct colors (or use your own color palette)

# Example: define your own color palette for genera
#custom_colors <-  c(setNames(rep("grey", length("grey_genera")),grey_genera),magma)



stacked <- ggplot(filtered_df, aes(x = station_ID, y = Proportional_Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = setNames(rainbow(length(unique(filtered_df$Genus))), unique(filtered_df$Genus)))+
  labs(x = "Station", y = "Proportional_Abundance", title = "Top 20 Genus Abundance")+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",           # or "bottom" if too wide
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    plot.margin = margin(1, 1, 1, 1, "cm")  # top, right, bottom, left
  )


print(stacked)


Joined_df_without_incretae <- joined_df[joined_df$Genus != "Incertae Sedis", ]



dev.off()

...................................................................................................

library(ggthemes)

library(dplyr)

library(vegan)



head(dist_mat)

dim(locationD_matrix_asv)
# e.g., 500 x 16 (ASVs x stations)

head(rownames(locationD_matrix_asv))  # should be ASV IDs
head(colnames(locationD_matrix_asv))  # should be station IDs

t_matrix <- t(locationD_matrix_asv)

dim(t_matrix)
head(rownames(t_matrix))  # these should be station IDs now
head(colnames(t_matrix))  # these should be ASV IDs now

#ordering the X axis according to clustering order

dist_mat <- vegdist(t_matrix^0.25, method = "bray")
clustering <- hclust(dist_mat, method = "ward.D")
head(clustering$labels)

# Get ordered site names
ordered_sites <- clustering$labels[clustering$order]

head(ordered_sites)


Joined_df_without_incretae$station_ID <- factor(Joined_df_without_incretae$station_ID, levels = ordered_sites)



# Sum abundances by Genus per site
genus_summary <- Joined_df_without_incretae %>%
  group_by(station_ID, Genus) %>%
  summarise(Abundance = sum(Proportional_Abundance), .groups = "drop")

# Get top 15 genera across all sites
top_genera <- genus_summary %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  slice_max(order_by = Total, n = 15) %>%
  pull(Genus)

# Filter original data to keep only top genera
filtered_genus_by_site <- genus_summary %>%
  filter(Genus %in% top_genera)

pdf("stacked_15 genus_plot_order.pdf", width = 10, height = 6)


library(ggplot2)
library(ggsci)

stacked_15_order <- ggplot(filtered_genus_by_site, aes(x = station_ID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_d3(palette = "category20")+
  labs(x = "Site", y = "Proportional Abundance", title = "Top 15 Genus-level Composition Across Sites") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(stacked_15_order)


dev.off()
...............................................................................................................

jpeg("stacked_15_genus_plot_order.jpeg", width = 1200, height = 800, res = 150)

stacked_15_order <- ggplot(filtered_genus_by_site, aes(x = station_ID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_d3(palette = "category20")+
  labs(x = "Site", y = "Proportional Abundance", title = "Top 15 Genus-level Composition Across Sites") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(stacked_15_order)


dev.off()
........................................................................................................................
#Stacked bar plot of families/phyla
jpeg("stacked_15_genus_plot_order.jpeg", width = 1200, height = 800, res = 150)

stacked_15_order <- ggplot(filtered_genus_by_site, aes(x = station_ID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_d3(palette = "category20")+
  labs(x = "Site", y = "Proportional Abundance", title = "Top 15 Genus-level Composition Across Sites") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(stacked_15_order)

....................................................................................................................................

#phylogenetic analysis


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

library(DECIPHER)





#pdf("stacked_wo_incretae.pdf", width = 10, height = 6)

#stacked_new <- ggplot(Joined_df_without_incretae, aes(x = station_ID, y = Proportional_Abundance, fill = Family)) +
  #geom_bar(stat = "identity") +
  #scale_fill_tableau(values = setNames(rainbow(length(unique(Joined_df_without_incretae$Family))), unique(Joined_df_without_incretae$Family)))+
  #labs(x = "Station", y = "Proportional_Abundance", title = "Family Abundance")+
  #theme_minimal() +
  #theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = "bottom",           # or "bottom" if too wide
    #legend.key.size = unit(0.4, "cm"),
    #legend.text = element_text(size = 8),
    #plot.margin = margin(1, 1, 1, 1, "cm")  # top, right, bottom, left
  #)


library(ggplot2)
library(ggthemes)




...........................................................................................



pdf(file=file.path(prjdir,"hm.pdf"),width=32,height=32)
hm = heatmap.2(as.matrix(num_data)^.25, trace="none",col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D"))
clusnum=8
oldclus = cutree(hm$rowDendrogram, k=clusnum)
oldorder = unname(rle(oldclus[as.hclust(hm$rowDendrogram)$order])$values)
neworder = (1:clusnum)
names(neworder) = oldorder
clus = unname(neworder[as.character(oldclus)])
names(clus) = names(oldclus)

hm = heatmap.2(as.matrix(num_data)^.25, trace="none", key=F, margins=c(20,30),col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D"),
               RowSideColors = rainbow(max(clus))[clus],
               #  labRow = meta.amun.small$project.date.station.depth.filter, cexRow=0.7,
               labRow = meta.tcn$description, cexRow=0.7,
               labCol = taxnames.small, cexCol=0.5)
dev.off()

heatmap.2(as.matrix(proptab.small)^.5, trace="none", key=F, margins=c(20,30),col=viridis,
          distfun=function(x) vegdist(x, method="bray"),
          hclustfun=function(x) hclust(x, method="ward.D"),
          RowSideColors = rainbow(max(clus))[clus],
          labRow = meta.amun.small$project.date.station.depth.filter, cexRow=0.7,
          labCol = taxnames.small, cexCol=0.5)

heatmap.2(as.matrix(proptab.small), trace="none", key=F, margins=c(20,30),
          distfun=function(x) vegdist(x, method="bray"),
          hclustfun=function(x) hclust(x, method="ward.D"),
          RowSideColors = rainbow(max(clus))[clus],
          labRow = meta.amun.small$project.date.station.depth.filter, cexRow=0.7,
          labCol = taxnames.small, cexCol=0.5)

for(taxon in unique(tax.16S[,"Order"])) {
  proptab.taxon = proptab.small
  rownames(proptab.taxon) = meta.amun.small$rownames
  taxon.pick = grep(taxon,taxnames.small)
  if(all(is.na(taxon.pick))) next
  proptab.taxon = proptab.taxon[!is.na(meta.amun.small$project.date.station.depth),taxon.pick,drop=F]
  proptab.taxon = proptab.taxon[rowSums(proptab.taxon>0)>0,,drop=F]
  if(!all(dim(proptab.taxon)>2)) next
  print(taxon)
  clus.taxon = clus[rownames(proptab.taxon)]
  #pdf(file=file.path(prjdir,"hm.taxon.pdf"), width=24,height=30)
  par(cex.main=3)
  heatmap.2(as.matrix(proptab.taxon)^.25, trace="none", key=F, margins=c(20,30),
            distfun=function(x) vegdist(x, method="bray"),
            hclustfun=function(x) hclust(x, method="ward.D"),
            RowSideColors = rainbow(max(clus))[clus.taxon],
            main=taxon,
            labRow = paste0(meta.amun.small[rownames(proptab.taxon),"date"]," ",meta.amun.small[rownames(proptab.taxon),"project.date.station.depth.filter"]), cexRow=0.7,
            labCol = taxnames.small[colnames(proptab.taxon)], cexCol=0.5)
  
}
dev.off()

# Save the plot as PDF with proper margins
pdf("heatmap_thamali.pdf", 
    width = 12,  # width in inches
    height = 12, # height in inches
    pointsize = 12)  # base font size

# Reset graphics state with larger margins
par(mar = c(15, 15, 5, 5))  # bottom, left, top, right margins

# Create the heatmap with adjusted margins
heatmap.2(heatmap_matrix^.25,  # Using power transformation like in TCN-plots.r
          trace = "none",
          col = viridis,
          distfun = function(x) vegdist(x, method = "bray"),
          hclustfun = function(x) hclust(x, method = "ward.D"),
          margins = c(15, 15),  # increased margins for labels
          cexRow = 0.7,
          cexCol = 0.5,
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          srtCol = 45,  # rotate column labels
          offsetCol = 0.5,  # offset for column labels
          offsetRow = 0.5,  # offset for row labels
          lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0)),  # adjust layout
          lhei = c(1.5, 4, 1),  # adjust height ratios
          lwid = c(1.5, 4, 1))  # adjust width ratios

# Close the PDF device
dev.off()


