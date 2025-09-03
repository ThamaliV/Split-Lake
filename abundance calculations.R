genus_level <- joined_df1 %>%
  group_by(Genus) %>%summarise(Total_Abundance= sum(Proportional_Abundance,na.rm=TRUE)) %>%
  #filter(!is.na(Genus)&Genus !="")
#head(genus_level)

mutate(percent_abundance = (Total_Abundance / sum(Total_Abundance)) * 100) %>%
  arrange(desc(percent_abundance)) %>%
  slice_head(n = 10)
print(genus_level)


phylum_level <- joined_df1 %>%
  group_by(Phylum) %>%summarise(Total_Abundance= sum(Proportional_Abundance,na.rm=TRUE)) %>%
  #filter(!is.na(Genus)&Genus !="")
  #head(genus_level)
  
  mutate(percent_abundance = (Total_Abundance / sum(Total_Abundance)) * 100) %>%
  arrange(desc(percent_abundance)) %>%
  slice_head(n = 10)
print(phylum_level)

family_level <- joined_df1 %>%
  group_by(Family) %>%summarise(Total_Abundance= sum(Proportional_Abundance,na.rm=TRUE)) %>%
  #filter(!is.na(Genus)&Genus !="")
  #head(genus_level)
  
  mutate(percent_abundance = (Total_Abundance / sum(Total_Abundance)) * 100) %>%
  arrange(desc(percent_abundance)) %>%
  slice_head(n = 10)
print(family_level)




cyano_pct <- joined_df1%>%
  summarise(
    total_counts   = sum(count, na.rm = TRUE),
    cyanobium_counts = sum(count[Genus == "Cyanobium PCC-6307"], na.rm = TRUE),
    cyanobium_percent = 100 * (cyanobium_counts / total_counts)
  )

asv_jsi <- joined_df1 %>% filter(station_ID == "JSI")

cyano_pct$cyanobium_percent


cyano_pct <- joined_df1 %>% filter(station_ID == "JSI", Genus=="Cyanobium PCC-6307")%>%
  summarise(Perc_Abundance= sum(Proportional_Abundance,na.rm=TRUE)*100) #%>%
  #filter(!is.na(Genus)&Genus !="")
  #head(genus_level)
  #mutate(percent_abundance = (Proportional_Abundance / sum(Proportional_Abundance)) * 100)
  #arrange(desc(percent_abundance)) %>%
  #slice_head(n = 10)
print(cyano_pct)

abundance <-  joined_df1 %>% 
          group_by(Genus) %>%
  summarise(total= sum(Proportional_Abundance))%>%
    mutate(percent=total/sum(total)*100) %>%
    filter(Genus == "Candidatus Methylopumilus")
print(abundance)

  
  
  
  
