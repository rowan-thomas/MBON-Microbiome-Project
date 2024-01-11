### Core Microbiome Analysis ###
#Acropora
acropora <- subset(Mbon_meta_all, host_taxon_abbreviation=="A. cervicornis")
acropora_polyp_subset <- subset(acropora, env_feature=="Polyp")
rownames(acropora_polyp_subset) <- acropora_polyp_subset$SampleID

#Creating a phyloseq object for acropora polyp samples
acropora_polyp_phy <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                              tax_table(as.data.frame(taxtable_16S) %>%
                                          column_to_rownames("Feature.ID") %>%
                                          as.matrix()), sample_data(acropora_polyp_subset))
acropora_polyp_phy = filter_taxa(acropora_polyp_phy, function(x) sum(x>5) > (0.15*length(x)),TRUE)
acropora_polyp_phy <- prune_taxa(taxa_sums(acropora_polyp_phy) > 0, acropora_polyp_phy)

#Calculate compositional version of the data:
acropora_polyp_rel <- microbiome::transform(acropora_polyp_phy, "compositional")

#Core microbiome at 90% prevalence
acropora_core_taxa_standard_90 <- core_members(acropora_polyp_rel, detection = 0, prevalence = 90/100)
acropora_core_taxa_standard_90 <- data.frame(acropora_core_taxa_standard_90)

#Core microbiome at 85% prevalence
acropora_core_taxa_standard_85 <- core_members(acropora_polyp_rel, detection = 0, prevalence = 85/100)
acropora_core_taxa_standard_85 <- data.frame(acropora_core_taxa_standard_85)

#Joining ASV sequences and taxonomy to core microbiome
acropora_core_taxa_standard_85 = rename(acropora_core_taxa_standard_85, Feature.ID = acropora_core_taxa_standard_85)
acropora_core_taxa_standard_85 <- inner_join(acropora_core_taxa_standard_85, seqs, by = "Feature.ID")
acropora_core_taxa_standard_85 <- inner_join(acropora_core_taxa_standard_85, taxtable_16S, by = "Feature.ID")

acropora_core_taxa_standard_90 = rename(acropora_core_taxa_standard_90, Feature.ID = acropora_core_taxa_standard_90)
acropora_core_taxa_standard_90 <- inner_join(acropora_core_taxa_standard_90, seqs, by = "Feature.ID")
acropora_core_taxa_standard_90 <- inner_join(acropora_core_taxa_standard_90, taxtable_16S, by = "Feature.ID")

#Orbicella 
orbicella <- subset(Mbon_meta_all, host_taxon_abbreviation=="O. faveolata")
orbicella_polyp_subset <- subset(orbicella, env_feature=="Polyp")
rownames(orbicella_polyp_subset) <- orbicella_polyp_subset$SampleID

#Creating a phyloseq object for orbicella polyp samples
orbicella_polyp_phy <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                               tax_table(as.data.frame(taxtable_16S) %>%
                                           column_to_rownames("Feature.ID") %>%
                                           as.matrix()), sample_data(orbicella_polyp_subset))
orbicella_polyp_phy = filter_taxa(orbicella_polyp_phy, function(x) sum(x>5) > (0.15*length(x)),TRUE)

#Calculate compositional version of the data:
orbicella_polyp_rel <- microbiome::transform(orbicella_polyp_phy, "compositional")

#Core microbiota analysis
orbicella_core_taxa_standard_90 <- core_members(orbicella_polyp_rel, detection = 0, prevalence = 90/100)
orbicella_core_taxa_standard_90 <- data.frame(orbicella_core_taxa_standard_90)

orbicella_core_taxa_standard_85 <- core_members(orbicella_polyp_rel, detection = 0, prevalence = 85/100)
orbicella_core_taxa_standard_85 <- data.frame(orbicella_core_taxa_standard_85)

#Joining ASV sequences and taxonomy to core microbiome
orbicella_core_taxa_standard_85 = rename(orbicella_core_taxa_standard_85, Feature.ID = orbicella_core_taxa_standard_85)
orbicella_core_taxa_standard_85 <- inner_join(orbicella_core_taxa_standard_85, seqs, by = "Feature.ID")
orbicella_core_taxa_standard_85 <- inner_join(orbicella_core_taxa_standard_85, taxtable_16S, by = "Feature.ID")

orbicella_core_taxa_standard_90 = rename(orbicella_core_taxa_standard_90, Feature.ID = orbicella_core_taxa_standard_90)
orbicella_core_taxa_standard_90 <- inner_join(orbicella_core_taxa_standard_90, seqs, by = "Feature.ID")
orbicella_core_taxa_standard_90 <- inner_join(orbicella_core_taxa_standard_90, taxtable_16S, by = "Feature.ID")

#Siderastrea
siderastrea <- subset(Mbon_meta_all, host_taxon_abbreviation=="S. siderea")
siderastrea_polyp_subset <- subset(siderastrea, env_feature=="Polyp")
rownames(siderastrea_polyp_subset) <- siderastrea_polyp_subset$SampleID

#Creating a phyloseq object for orbicella polyp samples
siderastrea_polyp_phy <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                                 tax_table(as.data.frame(taxtable_16S) %>%
                                             column_to_rownames("Feature.ID") %>%
                                             as.matrix()), sample_data(siderastrea_polyp_subset))
siderastrea_polyp_phy = filter_taxa(siderastrea_polyp_phy, function(x) sum(x>5) > (0.15*length(x)),TRUE)

#Calculate compositional version of the data:
siderastrea_polyp_rel <- microbiome::transform(siderastrea_polyp_phy, "compositional")

#Core microbiota analysis
siderastrea_core_taxa_standard_90 <- core_members(siderastrea_polyp_rel, detection = 0, prevalence = 90/100)
siderastrea_core_taxa_standard_90 <- data.frame(siderastrea_core_taxa_standard_90)

siderastrea_core_taxa_standard_85 <- core_members(siderastrea_polyp_rel, detection = 0, prevalence = 85/100)
siderastrea_core_taxa_standard_85 <- data.frame(siderastrea_core_taxa_standard_85)

#Joining ASV sequences and taxonomy to core microbiome
siderastrea_core_taxa_standard_85 = rename(siderastrea_core_taxa_standard_85, Feature.ID = siderastrea_core_taxa_standard_85)
siderastrea_core_taxa_standard_85 <- inner_join(siderastrea_core_taxa_standard_85, seqs, by = "Feature.ID")
siderastrea_core_taxa_standard_85 <- inner_join(siderastrea_core_taxa_standard_85, taxtable_16S, by = "Feature.ID")

siderastrea_core_taxa_standard_90 = rename(siderastrea_core_taxa_standard_90, Feature.ID = siderastrea_core_taxa_standard_90)
siderastrea_core_taxa_standard_90 <- inner_join(siderastrea_core_taxa_standard_90, seqs, by = "Feature.ID")
siderastrea_core_taxa_standard_90 <- inner_join(siderastrea_core_taxa_standard_90, taxtable_16S, by = "Feature.ID")

#Shared core microbiome OF/AC/SS
siderastrea_core_taxa_standard_85 <- select(siderastrea_core_taxa_standard_85, -c("Confidence", "Name", "ASVID"))
orbicella_core_taxa_standard_85 <- select(orbicella_core_taxa_standard_85, -c("Confidence", "Name", "ASVID"))
acropora_core_taxa_standard_85 <- select(acropora_core_taxa_standard_85, -c("Confidence", "Name", "ASVID"))

host_shared_core_microbiome_85 <- inner_join(siderastrea_core_taxa_standard_85, orbicella_core_taxa_standard_85, by=c("Feature.ID", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
host_shared_core_microbiome_85 <- inner_join(host_shared_core_microbiome_85, acropora_core_taxa_standard_85, by=c("Feature.ID", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

host_shared_core_microbiome_90 <- inner_join(siderastrea_core_taxa_standard_90, orbicella_core_taxa_standard_90, by=c("Feature.ID", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Confidence", "Name", "ASVID"))
host_shared_core_microbiome_90 <- inner_join(host_shared_core_microbiome_90, acropora_core_taxa_standard_90, by=c("Feature.ID", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Confidence", "Name", "ASVID"))
host_shared_core_microbiome_90 <- left_join(host_shared_core_microbiome_90, seqs, by="Feature.ID")

list_venn <- list("S. siderea" = siderastrea_core_taxa_standard_85$ASVID,
                  "O. faveolata" = orbicella_core_taxa_standard_85$ASVID,
                  "A. cervicornis" = acropora_core_taxa_standard_85$ASVID)

list_venn90 <- list("S. siderea" = siderastrea_core_taxa_standard_90$ASVID,
                  "O. faveolata" = orbicella_core_taxa_standard_90$ASVID,
                  "A. cervicornis" = acropora_core_taxa_standard_90$ASVID)

install.packages("ggvenn")

library(ggvenn)
library(RColorBrewer)

levels(Mbon_meta_all$host_taxon_abbreviation)
brewer.pal(n=3, "Paired")

host_core90 <- ggvenn(list_venn90,
       show_elements = FALSE,
       fill_color = c("#B2DF8A", "#1F78B4", "#A6CEE3"),
       fill_alpha = 0.3,
       stroke_linetype = "blank",
       set_name_size = 4,
       show_percentage = FALSE,
       text_size = 4)

host_core90

ggsave("~/Documents/University of Miami/MBON Project/figures/host_core_microbiome_venn.tif",
       device = tiff,
       width = 4, 
       height = 4, 
       units = "in",
       dpi=300)

#Sediment - Water - Polyp
#Sediment
sediment <- subset(Mbon_meta_all, env_feature=="Sediment")

sediment_phy <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                              tax_table(as.data.frame(taxtable_16S) %>%
                                          column_to_rownames("Feature.ID") %>%
                                          as.matrix()), sample_data(sediment))
sediment_phy= filter_taxa(sediment_phy, function(x) sum(x>5) > (0.15*length(x)),TRUE)
sediment_phy <- prune_taxa(taxa_sums(sediment_phy) > 0, sediment_phy)

#Calculate compositional version of the data:
sediment_rel <- microbiome::transform(sediment_phy, "compositional")

#Core microbiome at 90% prevalence
sediment_core_taxa_standard_90 <- core_members(sediment_rel, detection = 0, prevalence = 90/100)
sediment_core_taxa_standard_90 <- data.frame(sediment_core_taxa_standard_90)

#Joining ASV sequences and taxonomy to core microbiome
sediment_core_taxa_standard_90 = rename(sediment_core_taxa_standard_90, Feature.ID = sediment_core_taxa_standard_90)
sediment_core_taxa_standard_90 <- inner_join(sediment_core_taxa_standard_90, seqs, by = "Feature.ID")
sediment_core_taxa_standard_90 <- inner_join(sediment_core_taxa_standard_90, taxtable_16S, by = "Feature.ID", "Sequence")

#Water
water <- subset(Mbon_meta_all, env_feature=="Water")

water_phy <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                        tax_table(as.data.frame(taxtable_16S) %>%
                                    column_to_rownames("Feature.ID") %>%
                                    as.matrix()), sample_data(water))
water_phy= filter_taxa(water_phy, function(x) sum(x>5) > (0.15*length(x)),TRUE)
water_phy <- prune_taxa(taxa_sums(water_phy) > 0, water_phy)

#Calculate compositional version of the data:
water_rel <- microbiome::transform(water_phy, "compositional")

#Core microbiome at 90% prevalence
water_core_taxa_standard_90 <- core_members(water_rel, detection = 0, prevalence = 90/100)
water_core_taxa_standard_90 <- data.frame(water_core_taxa_standard_90)

#Joining ASV sequences and taxonomy to core microbiome
water_core_taxa_standard_90 = rename(water_core_taxa_standard_90, Feature.ID = water_core_taxa_standard_90)
water_core_taxa_standard_90 <- inner_join(water_core_taxa_standard_90, seqs, by = "Feature.ID")
water_core_taxa_standard_90 <- inner_join(water_core_taxa_standard_90, taxtable_16S, by = "Feature.ID", "Sequence")

#Polyp
polyp <- subset(Mbon_meta_all, env_feature=="Polyp")

polyp_phy <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                     tax_table(as.data.frame(taxtable_16S) %>%
                                 column_to_rownames("Feature.ID") %>%
                                 as.matrix()), sample_data(polyp))
polyp_phy= filter_taxa(polyp_phy, function(x) sum(x>5) > (0.15*length(x)),TRUE)
polyp_phy <- prune_taxa(taxa_sums(polyp_phy) > 0, polyp_phy)

#Calculate compositional version of the data:
polyp_rel <- microbiome::transform(polyp_phy, "compositional")

#Core microbiome at 90% prevalence
polyp_core_taxa_standard_90 <- core_members(polyp_rel, detection = 0, prevalence = 90/100)
polyp_core_taxa_standard_90 <- data.frame(polyp_core_taxa_standard_90)

#Joining ASV sequences and taxonomy to core microbiome
polyp_core_taxa_standard_90 = rename(polyp_core_taxa_standard_90, Feature.ID = polyp_core_taxa_standard_90)
polyp_core_taxa_standard_90 <- inner_join(polyp_core_taxa_standard_90, seqs, by = "Feature.ID")
polyp_core_taxa_standard_90 <- inner_join(polyp_core_taxa_standard_90, taxtable_16S, by = "Feature.ID", "Sequence")

SAMPLE_shared_core_microbiome <- inner_join(polyp_core_taxa_standard_90, water_core_taxa_standard_90, by=c("Feature.ID"))
SAMPLE_shared_core_microbiome <- inner_join(SAMPLE_shared_core_microbiome, sediment_core_taxa_standard_90, by=c("Feature.ID"))

water_polyp_shared <- inner_join(polyp_core_taxa_standard_90, water_core_taxa_standard_90, by=c("Feature.ID"))

list_venn_env <- list("Sediment" = sediment_core_taxa_standard_90$ASVID,
                  "Water" = water_core_taxa_standard_90$ASVID,
                  "Polyp" = polyp_core_taxa_standard_90$ASVID)

ggvenn(list_venn_env,
       show_elements = FALSE,
       fill_color =c("#E69F00", "#56B4E9", "#AE0F1B"),
       fill_alpha = 0.3,
       stroke_linetype = "blank",
       set_name_size = 4,
       show_percentage = FALSE,
       text_size = 4)

ggsave("~/Documents/University of Miami/MBON Project/figures/env_core_microbiome_venn.tif",
       device = tiff,
       width = 4, 
       height = 4, 
       units = "in",
       dpi=300)
