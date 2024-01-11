### Loading Necessary Libraries ###

library("pairwiseAdonis")
library("devtools")
library("qiime2R")
library("tidyverse")
library("vegan")
library("ANCOMBC")
library("phyloseq")
library("dplyr")
library("tibble")
library("microbiome") 
library("ggplot2")
library("DT")
library("microseq")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

### Importing the Files ###
#Metadata file
Mbon_meta_all <- read.csv("~/Documents/University of Miami/MBON Project/Metadata/Mbon_meta_all.csv")
Mbon_meta_all <- Mbon_meta_all %>% rename(SampleID = X.SampleID)

ASVtable <- read_qza("~/Documents/University of Miami/MBON Project/qiime/table/table-merged-no-MitoChl-filter.qza")
ASVtable_data <- ASVtable$data

#Taxonomy data
ASVtaxa <- read_qza("~/Documents/University of Miami/MBON Project/qiime/taxonomy/taxonomy_MBON.qza")
taxtable_16S <- ASVtaxa$data %>% as_tibble() %>% 	separate(Taxon, sep=";", c("Kingdom", "Phylum", 	"Class", "Order", "Family", "Genus", "Species"))

taxtable_16S$Kingdom<-gsub("d__","",as.character(taxtable_16S$Kingdom))
taxtable_16S$Phylum<-gsub("p__","",as.character(taxtable_16S$Phylum))
taxtable_16S$Class<-gsub("c__","",as.character(taxtable_16S$Class))
taxtable_16S$Order<-gsub("o__","",as.character(taxtable_16S$Order))
taxtable_16S$Family<-gsub("f__","",as.character(taxtable_16S$Family))
taxtable_16S$Genus<-gsub("g__","",as.character(taxtable_16S$Genus))
taxtable_16S$Species<-gsub("s__","",as.character(taxtable_16S$Species))

taxtable_16S <- taxtable_16S %>%
  mutate(Name = str_c(Genus, Species))

taxtable_16S$ASVID <- 1:nrow(taxtable_16S)
taxtable_16S$ASVID <- paste("ASV", taxtable_16S$ASVID)

#Sequence reads
seqs <- readFasta("~/Documents/University of Miami/MBON Project/qiime/rep-seqs/sequences.fasta.html")
colnames(seqs) <- c("Feature.ID","Sequence")

### Cleaning up the Metadata ###

#Sample Site
Mbon_meta_all %>% group_by(sample_site_id) %>% summarize(SampleID=n())
Mbon_meta_all <- Mbon_meta_all %>% mutate(sample_site_id = recode(sample_site_id, 
                                                                  Mr = "MR", 
                                                                  Sr = "SR", 
                                                                  Ws = "WS",
                                                                  Tw = "TN",
                                                                  CH = "CR",
                                                                  TR = "TN"))

#Recode month numbers as names
Mbon_meta_all <- Mbon_meta_all %>% mutate(month = recode(month, 
                                                         "2" = "February", 
                                                         "4" = "April", 
                                                         "7" = "July",
                                                         "11" = "November",
                                                         "12" = "December"))

#Create new column in metadata with key designation
Mbon_meta_all <- Mbon_meta_all %>% add_column(Key = NA)
Mbon_meta_all$Key <- recode(Mbon_meta_all$sample_site_id, 
                            ER = "Upper", 
                            MR = "Upper",
                            SR = "Middle", 
                            TN = "Middle",
                            CR = "Middle", 
                            LK = "Lower", 
                            WS = "Lower")

Mbon_meta_all %>% group_by(Key) %>% summarize(SampleID=n())
Mbon_meta_all$Key <- as.factor(Mbon_meta_all$Key)

#Recode reef sites as names
Mbon_meta_all <- Mbon_meta_all %>% mutate(sample_site_id = recode(sample_site_id, 
                                                         "CR" = "Cheeca Rocks", 
                                                         "LK" = "Looe Key", 
                                                         "MR" = "Molasses Reef",
                                                         "SR" = "Sombrero Reef",
                                                         "TN" = "Tennessee Reef",
                                                         "WS" = "Western Sambo",
                                                         "ER" = "Emerald Reef"))

#Recode species as names
Mbon_meta_all <- Mbon_meta_all %>% mutate(host_taxon_abbreviation = recode(host_taxon_abbreviation, 
                                                                  "AC" = "A. cervicornis", 
                                                                  "OF" = "O. faveolata", 
                                                                  "SS" = "S. siderea"))

#Summarize by species
Mbon_meta_all %>% group_by(host_taxon_abbreviation) %>% summarize(SampleID=n())
Mbon_meta_all <- Mbon_meta_all %>% mutate(host_taxon_abbreviation = recode(host_taxon_abbreviation, AG = "AC"))
Mbon_meta_all %>% group_by(host_taxon_abbreviation) %>% summarize(SampleID=n())
Mbon_meta_all$host_taxon_abbreviation <- as.factor(Mbon_meta_all$host_taxon_abbreviation)

#Summarize by environment feature
Mbon_meta_all %>% group_by(env_feature) %>% summarize(SampleID = n())
Mbon_meta_all <- Mbon_meta_all %>% mutate(env_feature = recode(env_feature, polyp = "Polyp", sediment = "Sediment", water = "Water"))
Mbon_meta_all %>% group_by(env_feature) %>% summarize(SampleID = n())
Mbon_meta_all$env_feature <- as.factor(Mbon_meta_all$env_feature)

#Recoding host taxon name column
Mbon_meta_all$host_taxon_name <- recode(Mbon_meta_all$host_taxon_abbreviation, "A. cervicornis" = "Acropora cervicornis", "O. faveolata" = "Orbicella faveolata", "S. siderea" = "Siderastrea siderea")

#Recode na months as NA and remove
Mbon_meta_all$month[Mbon_meta_all$month=="na"]=NA #renames invalid values as NA
Mbon_meta_all <- na.omit(Mbon_meta_all)

Mbon_meta_all$month <- as.factor(Mbon_meta_all$month)
Mbon_meta_all$month <- factor(Mbon_meta_all$month, levels = c("February", "April", "July", "November", "December"))
levels(Mbon_meta_all$month)

### Calculating Alpha  Diversity ###
library(ggpubr)
library(rcompanion)
theme_set(theme_bw())

rownames(Mbon_meta_all) <- Mbon_meta_all$SampleID 

#Creating a phyloseq object
physeq_all <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                      tax_table(as.data.frame(taxtable_16S) %>%
                                  column_to_rownames("Feature.ID") %>%
                                  as.matrix()), sample_data(Mbon_meta_all))

#Rarefy the data
ps_rare <- phyloseq::rarefy_even_depth(physeq_all, rngseed = 123, replace = TRUE, sample.size=1000)           

#Estimating richness
richness_rarefied <- estimate_richness(ps_rare)
richness_rarefied <- tibble::rownames_to_column(richness_rarefied, "SampleID")

#Subsetting Shannon diversity metrics
Shannon_alpha_diversity_rarefied <- data.frame(richness_rarefied$SampleID, richness_rarefied$Shannon)
Shannon_alpha_diversity_rarefied = Shannon_alpha_diversity_rarefied %>%
  rename(SampleID = richness_rarefied.SampleID,
         Shannon = richness_rarefied.Shannon)

#Join shannon values to metadata
Shannon_alpha_diversity_rarefied$SampleID = gsub(".", "-", Shannon_alpha_diversity_rarefied$SampleID, fixed = TRUE)
Mbon_meta_all_with_alpha <- inner_join(Shannon_alpha_diversity_rarefied, Mbon_meta_all, by = "SampleID")

#Checking normalcy
plotNormalHistogram(Shannon_alpha_diversity_rarefied$Shannon, main="Shannon index", xlab="") #not normal

#Plotting alpha diversity by reef site
alpha_reef <- plot_richness(ps_rare, x="sample_site_id", measures="Shannon", color = "sample_site_id")+
  geom_boxplot(alpha=0.5)+ 
  theme_bw()+
  labs(color = "Reef Site")+
  theme(
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size=15),
        axis.title.y = element_text(vjust = +5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),)+
  coord_cartesian(ylim = c(0, 7.5))+
  ylab("Shannon Index Value")+
  scale_color_brewer(palette = "Paired")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

alpha_reef

#Plotting alpha diversity by host species
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "not significant"))
comparisons_host <- list( c("A. cervicornis", "O. faveolata"), c("A. cervicornis", "S. siderea"), c("S. siderea", "O. faveolata"))

alpha_hosts <- plot_richness(ps_rare, x="host_taxon_abbreviation", measures="Shannon", color = "host_taxon_abbreviation")+
  geom_boxplot(alpha=0.5)+ 
  theme_bw()+
  labs(color = "Host Species")+
  theme(
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size=15),
        axis.title.y = element_text(vjust = +5),
        axis.title.x = element_blank(),
        plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),)+
  coord_cartesian(ylim = c(0, 7.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_host, label = "p.signif", symnum.args = symnum.args, size=3)+
  ylab("Shannon Index Value")+
  scale_color_brewer(palette = "Paired")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

alpha_hosts

ggsave("~/Documents/University of Miami/MBON Project/figures/alpha_hosts.png",
       width = 6, height = 4, units = "in", dpi=300)

#Calculate mean alpha diversity values for each host species
aggregate(Mbon_meta_all_with_alpha$Shannon, list(Mbon_meta_all_with_alpha$host_taxon_abbreviation), FUN=mean) 

#Plotting alpha diversity by sample type
comparisons_env <- list( c("Polyp", "Sediment"), c("Polyp", "Water"), c("Sediment", "Water"))
alpha_environment <- plot_richness(ps_rare, x="env_feature", measures="Shannon", color = "env_feature")+
  geom_boxplot(alpha=0.5)+ 
  theme_bw()+
  labs(color = "Sample Type")+
  theme(
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size=15),
        axis.title.y = element_text(vjust = +5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),)+
  coord_cartesian(ylim = c(0, 7.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_env, label = "p.signif", symnum.args = symnum.args)+
  ylab("Shannon Index Value")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())+
  scale_color_manual(values=c("#AE0F1B", "#E69F00", "#56B4E9"))

alpha_environment

#Calculating mean alpha diversity values for sample types
aggregate(Mbon_meta_all_with_alpha$Shannon, list(Mbon_meta_all_with_alpha$env_feature), FUN=mean) 

#Calculating wilcoxon test stats for alpha diversity
library(tidyverse)
library(rstatix)
library(ggpubr)\
library(writexl)

stat.test.env <- Mbon_meta_all_with_alpha %>% 
  wilcox_test(Shannon ~ env_feature)
stat.test.env = select(stat.test.env, -".y.")
write_xlsx(stat.test.env,"~/Documents/University of Miami/MBON Project/figures/stat.test.env.xlsx")

stat.test.reef <- Mbon_meta_all_with_alpha %>% 
  wilcox_test(Shannon ~ sample_site_id)
stat.test.reef = select(stat.test.reef, -".y.")
write_xlsx(stat.test.reef,"~/Documents/University of Miami/MBON Project/figures/stat.test.reef.xlsx")

stat.test.host <- Mbon_meta_all_with_alpha %>% 
  wilcox_test(Shannon ~ host_taxon_abbreviation)
stat.test.host = select(stat.test.host, -".y.")
view(stat.test.host)
write_xlsx(stat.test.host,"~/Documents/University of Miami/MBON Project/figures/stat.test.host.xlsx")

### Beta Diversity Analysis ###
#Reef site
#Creating NMDS plots to show beta diversity between reef sites
rownames(Mbon_meta_all) <- Mbon_meta_all$SampleID 
physeq_16S <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                      tax_table(as.data.frame(taxtable_16S) %>%
                                  column_to_rownames("Feature.ID") %>%
                                  as.matrix()), sample_data(Mbon_meta_all))
physeq_16S = filter_taxa(physeq_16S, function(x) sum(x>5) > (0.15*length(x)),TRUE)

#NMDS Plots
RA_all = microbiome::transform(physeq_16S, "compositional")
all_ord = ordinate(RA_all, "NMDS", "bray")

reef_NMDS <- plot_ordination(physeq_16S, all_ord, color="sample_site_id", shape = "env_feature") +
  geom_point(size = 2, )+
  labs(color = "Reef Site")+
  labs(shape = "Sample Type")+
  scale_color_brewer(palette = "Paired")+
  scale_shape_manual(values = c(1, 7, 8)) + 
  theme_classic()+
  theme(axis.title.x=element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

reef_NMDS

reefs_combined <- ggarrange(alpha_reef, reef_NMDS,
                            labels = c("A", "B"),
                            ncol = 2, nrow = 1,
                            common.legend = TRUE,
                            legend = "bottom")
reefs_combined

ggsave("~/Documents/University of Miami/MBON Project/figures/reefs_combined.tif",
       device = tiff,
       width = 7.5, 
       height = 4, 
       units = "in", 
       dpi=300)

reef_pairwise <- pairwise.adonis(t(otu_table(RA_all)), sample_data(RA_all)$sample_site_id, sim.method 	= "bray", p.adjust.m = "bonferroni") 
view(reef_pairwise)
write_xlsx(reef_pairwise,"~/Documents/University of Miami/MBON Project/figures/reef_pairwise.xlsx")

key_pairwise <- pairwise.adonis(t(otu_table(RA_all)), sample_data(RA_all)$Key, sim.method 	= "bray", p.adjust.m = "bonferroni") 
view(key_pairwise)
write_xlsx(key_pairwise,"~/Documents/University of Miami/MBON Project/figures/key_pairwise.xlsx")

#Host Species
#Subset only samples taken from polyps
polyp <- subset(Mbon_meta_all, env_feature=="Polyp")

#Creating NMDS plots to show beta diversity between species
polyp <- subset(Mbon_meta_all, env_feature=="Polyp")
polyp_physeq <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                           tax_table(as.data.frame(taxtable_16S) %>%
                                       column_to_rownames("Feature.ID") %>%
                                       as.matrix()), sample_data(polyp))
polyp_physeq_filtered = filter_taxa(polyp_physeq, function(x) sum(x>5) > (0.15*length(x)),TRUE)

polyp_RA = microbiome::transform(polyp_physeq_filtered, "compositional")
polyp_ord = ordinate(polyp_RA, "NMDS", "bray")

host_NMDS <- plot_ordination(polyp_physeq_filtered, polyp_ord, color="host_taxon_abbreviation") +
  geom_point(size = 2, )+
  labs(color = "Host Species")+
  scale_color_brewer(palette = "Paired")+
  scale_shape_manual(values = c(1, 7, 8)) + 
  theme_classic()+
  theme(axis.title.x=element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

host_NMDS

hosts_combined <- ggarrange(alpha_hosts, host_NMDS,
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1,
                       common.legend = TRUE,
                       legend = "bottom")
hosts_combined

ggsave("~/Documents/University of Miami/MBON Project/figures/hosts_combined.tif",
       device = tiff,
       width = 7.5, 
       height = 4, 
       units = "in", 
       dpi=300)


#Calculating pairwise comparisons for host species
polyp_host_pairwise <- pairwise.adonis(t(otu_table(polyp_RA)), sample_data(polyp_RA)$host_taxon_abbreviation, sim.method 	= "bray", p.adjust.m = "bonferroni") 
view(polyp_host_pairwise)
write_xlsx(polyp_host_pairwise,"~/Documents/University of Miami/MBON Project/figures/host_pairwise.xlsx")


#Beta Diversity Between Sample Types
rownames(Mbon_meta_all) <- Mbon_meta_all$SampleID 
physeq_16S <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                      tax_table(as.data.frame(taxtable_16S) %>%
                                  column_to_rownames("Feature.ID") %>%
                                  as.matrix()), sample_data(Mbon_meta_all))
physeq_16S = filter_taxa(physeq_16S, function(x) sum(x>5) > (0.15*length(x)),TRUE)

#NMDS Plots
RA_all = microbiome::transform(physeq_16S, "compositional")
all_ord = ordinate(RA_all, "NMDS", "bray")

sample_NMDS <- plot_ordination(ps_rare, all_ord, color="env_feature", shape = "host_taxon_abbreviation") +
  geom_point(size = 2.5)+
  labs(color = "Sample Type", shape = "Host Species")+
  scale_color_manual(values=c("#AE0F1B", "#E69F00", "#56B4E9"))+
  scale_shape_manual(values = c(1, 7, 8)) + 
  theme_classic()+
  theme(axis.title.x=element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

sample_NMDS

ggsave("NMDS_sample_type_presentation2.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

all_sampletype_pairwise <- pairwise.adonis(t(otu_table(RA_all)), sample_data(RA_all)$env_feature, sim.method 	= "bray", p.adjust.m = "bonferroni") 
view(all_sampletype_pairwise)
write_xlsx(all_sampletype_pairwise,"~/Documents/University of Miami/MBON Project/figures/sampletype_pairwise.xlsx")


samples_combined <- ggarrange(alpha_environment, sample_NMDS,
                            labels = c("A", "B"),
                            ncol = 2, nrow = 1,
                            common.legend = TRUE,
                            legend = "bottom")


samples_combined

ggsave("~/Documents/University of Miami/MBON Project/figures/samples_combined.tif",
       device = tiff,
       width = 7.5, 
       height = 4, 
       units = "in", 
       dpi=300)

#Testing homogeneity
library(vegan)
ps_clr = microbiome::transform(polyp_physeq, "clr")

dis <- vegdist(otu_table(t(ps_clr)), method ="euclidean")

#homogeneity of group dispersions (variances).
mod <- betadisper(dis, sample_data(ps_clr)$host_taxon_abbreviation)
mod
## Permutation test for F
permutest(mod)
TukeyHSD(mod)
