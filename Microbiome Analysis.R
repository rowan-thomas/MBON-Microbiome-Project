### Loading Necessary Libraries ###
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
library("devtools")
library("pairwiseAdonis")
library("qiime2R")
library("tidyverse")
library("vegan")
library("ANCOMBC")
library("phyloseq")
library("dplyr")
library("tibble")
library("microbiome") 
library("ggplot2")
library("dplyr")
library("DT")
library("microseq")

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

Mbon_meta_all %>% group_by(host_taxon_abbreviation) %>% summarize(SampleID=n())

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
Mbon_meta_all$Key <- recode(Mbon_meta_all$sample_site_id, CR = "Upper", ER = "Upper", TR = "Upper", MR = "Upper",
                            LK = "Lower", WS = "Lower",
                            SR = "Middle", TN = "Lower",
                            CH = "NA", Tw = "NA")
Mbon_meta_all %>% group_by(Key) %>% summarize(SampleID=n())
Mbon_meta_all$Key <- as.factor(Mbon_meta_all$Key)

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
Mbon_meta_all$host_taxon_name <- recode(Mbon_meta_all$host_taxon_abbreviation, AC = "Acropora cervicornis", OF = "Orbicella faveolata", SS = "Siderastrea siderea")




### Alpha  Diversity ###
library(ggpubr)
library(rcompanion)
library(tibble)
theme_set(theme_bw())

rownames(Mbon_meta_all) <- Mbon_meta_all$SampleID 

#Creating a phyloseq object
physeq_all <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                      tax_table(as.data.frame(taxtable_16S) %>%
                                  column_to_rownames("Feature.ID") %>%
                                  as.matrix()), sample_data(Mbon_meta_all))

#Filtering
physeq_all = filter_taxa(physeq_all, function(x) sum(x>5) > (0.15*length(x)),TRUE)

#Estimating richness
richness_all <-estimate_richness(physeq_all)
richness_all <- tibble::rownames_to_column(richness_all, "SampleID")

#Subsetting Shannon diversity metrics
Shannon_alpha_diversity_all <- data.frame(richness_all$SampleID, richness_all$Shannon)
Shannon_alpha_diversity_all = Shannon_alpha_diversity_all %>%
  rename(SampleID = richness_all.SampleID,
         Shannon = richness_all.Shannon)

#Join shannon to metadata
Shannon_alpha_diversity_all$SampleID = gsub(".", "-", Shannon_alpha_diversity_all$SampleID, fixed = TRUE)

Mbon_meta_all_with_alpha <- inner_join(Shannon_alpha_diversity_all, Mbon_meta_all, by = "SampleID")

#Checking normalcy
plotNormalHistogram(Shannon_alpha_diversity_all$Shannon, main="Shannon index", xlab="") #not normal

#Plotting alpha diversity by host species
par(mfrow = c(2,1))

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "not significant"))
comparisons_host <- list( c("AC", "OF"), c("AC", "SS"), c("SS", "OF"))
alpha_hosts <- plot_richness(physeq_all, x="host_taxon_abbreviation", measures="Shannon", color = "host_taxon_abbreviation")+
  geom_boxplot(alpha=0.5)+ 
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=45,hjust=1,vjust=1, size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size=15),
        axis.title.y = element_text(vjust = +5),
        axis.title.x = element_text(vjust = -2.5),
        plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),)+
  coord_cartesian(ylim = c(0, 6.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_host, label = "p.signif", symnum.args = symnum.args)+
  xlab("Host Species")+
  scale_color_brewer(palette = "Paired")+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text.x= element_text(size=15))

library(tidyverse)
library(rstatix)
library(ggpubr)

stat.test <- Mbon_meta_all_with_alpha %>% 
  wilcox_test(Shannon ~ env_feature)

view(stat.test)
stat.test = select(stat.test, -".y.")

gridExtra::grid.table(stat.test)


view(Mbon_meta_all_with_alpha)

aggregate(Mbon_meta_all_with_alpha$Shannon, list(Mbon_meta_all_with_alpha$host_taxon_abbreviation), FUN=mean) 

alpha_hosts
ggsave("alpha_diversity_host_species.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

#Plotting alpha diversity by sample type
comparisons_env <- list( c("Polyp", "Sediment"), c("Polyp", "Water"), c("Sediment", "Water"))
alpha_environment <- plot_richness(physeq_all, x="env_feature", measures="Shannon", color = "env_feature")+
  geom_boxplot(alpha=0.5)+ 
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=45,hjust=1,vjust=1, size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size=15),
        axis.title.y = element_text(vjust = +5),
        axis.title.x = element_text(vjust = -2.5),
        plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),)+
  coord_cartesian(ylim = c(0, 6.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_env, label = "p.signif", symnum.args = symnum.args)+
  xlab("Environment Type")+
  theme(strip.background =element_rect(fill="white"))+
  scale_color_manual(values=c("#AE0F1B", "#E69F00", "#56B4E9"))

alpha_environment
ggsave("alpha_diversity_environment_type.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

ggsave("alpha_diversity_combined.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

aggregate(Mbon_meta_all_with_alpha$Shannon, list(Mbon_meta_all_with_alpha$env_feature), FUN=mean) 


### Beta Diversity Analysis ###

#Part 1: Beta Diversity in Polyps Comparing Host Species
#Subset only samples taken from polyps
polyp <- subset(Mbon_meta_all, env_feature=="Polyp")

#Create the polyp phyloseq object
rownames(polyp) <- polyp$SampleID
polyp_physeq <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                        tax_table(as.data.frame(taxtable_16S) %>%
                                    column_to_rownames("Feature.ID") %>%
                                    as.matrix()), sample_data(polyp))

#Filtering
polyp_ps_filter = filter_taxa(polyp_physeq, function(x) sum(x>10) > (0.30*length(x)),TRUE)

#CLR transformation
polyp_clr = microbiome::transform(polyp_ps_filter, "clr")
polyp_ra <- microbiome::transform(polyp_ps_filter, 'compositional') %>%
  tax_glom("Order")
ps_ra_melt =  psmelt(polyp_ra) 
dim(ps_ra_melt)

ps_ra_all = transform_sample_counts(polyp_ps_filter, function(x) x / sum(x)) %>%
  tax_glom("Family")%>% psmelt()

library(RColorBrewer)
nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(22, "Paired"))(nb.cols)

RA_plot_Order = ps_ra_all %>%
  subset(Abundance > 0.04) %>%
  ggplot(  
    aes(x = sample_site_id, y=Abundance, fill=Order)) + 
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  facet_wrap(~ host_taxon_name, scales="free") +
  theme(axis.title.x=element_text(size = 12, vjust = -1),
        axis.title.y = element_text(size = 12, vjust = +2),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size =10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'))+
  theme(strip.text.x = element_text(size = 10))+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Relative Abundance") +
  xlab("Reef")

RA_plot_Order

### ACROPORA POLYP RELATIVE ABUNDANCE ###
polyp <- subset(Mbon_meta_all, env_feature=="Polyp")
AC_polyp <- subset(Mbon_meta_all, host_taxon_abbreviation=="AC")

#Create the polyp phyloseq object
rownames(AC_polyp) <- AC_polyp$SampleID
AC_polyp_physeq <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                        tax_table(as.data.frame(taxtable_16S) %>%
                                    column_to_rownames("Feature.ID") %>%
                                    as.matrix()), sample_data(AC_polyp))

#Filtering
AC_polyp_ps_filter = filter_taxa(AC_polyp_physeq, function(x) sum(x>10) > (0.30*length(x)),TRUE)

#CLR transformation
AC_polyp_clr = microbiome::transform(AC_polyp_ps_filter, "clr")
AC_polyp_ra <- microbiome::transform(AC_polyp_ps_filter, 'compositional') %>%
  tax_glom("Order")
AC_ps_ra_melt =  psmelt(AC_polyp_ra) 

AC_ps_ra_all = transform_sample_counts(AC_polyp_ps_filter, function(x) x / sum(x)) %>%
  tax_glom("Family")%>% psmelt()

library(RColorBrewer)
nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(22, "Paired"))(nb.cols)

RA_plot_Order_AC = AC_ps_ra_all %>%
  subset(Abundance > 0.04) %>%
  ggplot(  
    aes(x = month, y=Abundance, fill=Order)) + 
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  facet_wrap(~ host_taxon_name, scales="free") +
  theme(axis.title.x=element_text(size = 12, vjust = -1),
        axis.title.y = element_text(size = 12, vjust = +2),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, size =10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'))+
  theme(strip.text.x = element_text(size = 10))+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Relative Abundance") +
  xlab("Month")

RA_plot_Order_AC

ggsave("AC_polyp_RA.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

### ORBICELLA POLYP RELATIVE ABUNDANCE ###
polyp <- subset(Mbon_meta_all, env_feature=="Polyp")
OF_polyp <- subset(Mbon_meta_all, host_taxon_abbreviation=="OF")
OF_polyp$month[OF_polyp$month=="na"]=NA #renames invalid values as NA
OF_polyp <- na.omit(OF_polyp)

OF_polyp$month <- as.factor(OF_polyp$month)
OF_polyp$month <- factor(OF_polyp$month, levels = c("February", "April", "July", "November", "December"))
levels(OF_polyp$month)

#Create the polyp phyloseq object
rownames(OF_polyp) <- OF_polyp$SampleID
OF_polyp_physeq <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                           tax_table(as.data.frame(taxtable_16S) %>%
                                       column_to_rownames("Feature.ID") %>%
                                       as.matrix()), sample_data(OF_polyp))

#Filtering
OF_polyp_ps_filter = filter_taxa(OF_polyp_physeq, function(x) sum(x>10) > (0.30*length(x)),TRUE)

#CLR transformation
OF_polyp_clr = microbiome::transform(OF_polyp_ps_filter, "clr")
OF_polyp_ra <- microbiome::transform(OF_polyp_ps_filter, 'compositional') %>%
  tax_glom("Order")
OF_ps_ra_melt =  psmelt(OF_polyp_ra) 

OF_ps_ra_all = transform_sample_counts(OF_polyp_ps_filter, function(x) x / sum(x)) %>%
  tax_glom("Order")%>% psmelt()

library(RColorBrewer)
nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(22, "Paired"))(nb.cols)

RA_plot_Order_OF = OF_ps_ra_all %>%
  subset(Abundance > 0.04) %>%
  ggplot(  
    aes(x = month, y=Abundance, fill=Order)) + 
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=c("#234d20",
                             "#36802d", 
                             "#77ab59",
                             "#e8f4ea",
                             "#f0f7da", 
                             "#005073", 
                             "#107dac", 
                             "#189ad3",
                             "#1ebbd7",
                             "#71c7ec",
                             "#ffbaba", 
                             "#ff7b7b",
                             "#ff0000",
                             "#a70000",
                             "#efbbff", 
                             "#d896ff", 
                             "#be29ec", 
                             "#ffc2cd",
                             "#fff9ae",
                             "#ffd7b5",
                             "#ffb38a",
                             "#ff9248",
                             "#ff6700",
                             "#a3c1ad")) +
  theme_classic() +
  facet_wrap(~ host_taxon_name, scales="free") +
  theme(axis.title.x=element_text(size = 12, vjust = -1),
        axis.title.y = element_text(size = 12, vjust = +2),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(angle = 45, size =10, hjust = 1, vjust =1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'))+
  theme(strip.text.x = element_text(size = 10))+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Relative Abundance") +
  xlab("Month")

RA_plot_Order_OF

ggsave("OF_polyp_RA.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

### SIDERASTREA POLYP RELATIVE ABUNDANCE ###
polyp <- subset(Mbon_meta_all, env_feature=="Polyp")
SS_polyp <- subset(Mbon_meta_all, host_taxon_abbreviation=="SS")
SS_polyp$month[SS_polyp$month=="na"]=NA #renames invalid values as NA
SS_polyp <- na.omit(SS_polyp)

SS_polyp$month <- as.factor(SS_polyp$month)
SS_polyp$month <- factor(SS_polyp$month, levels = c("February", "April", "July", "November", "December"))
levels(SS_polyp$month)

#Create the polyp phyloseq object
rownames(SS_polyp) <- SS_polyp$SampleID
SS_polyp_physeq <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                           tax_table(as.data.frame(taxtable_16S) %>%
                                       column_to_rownames("Feature.ID") %>%
                                       as.matrix()), sample_data(SS_polyp))

#Filtering
SS_polyp_ps_filter = filter_taxa(SS_polyp_physeq, function(x) sum(x>10) > (0.30*length(x)),TRUE)

#CLR transformation
SS_polyp_clr = microbiome::transform(SS_polyp_ps_filter, "clr")
SS_polyp_ra <- microbiome::transform(SS_polyp_ps_filter, 'compositional') %>%
  tax_glom("Order")
SS_ps_ra_melt =  psmelt(SS_polyp_ra) 

SS_ps_ra_all = transform_sample_counts(SS_polyp_ps_filter, function(x) x / sum(x)) %>%
  tax_glom("Order")%>% psmelt()

library(RColorBrewer)
nb.cols <- 24

RA_plot_Order_SS = SS_ps_ra_all %>%
  subset(Abundance > 0.04) %>%
  ggplot(  
    aes(x = month, y=Abundance, fill=Order)) + 
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=c("#234d20",
                             "#36802d", 
                             "#77ab59",
                             "#e8f4ea",
                             "#f0f7da", 
                             "#005073", 
                             "#107dac", 
                             "#189ad3",
                             "#1ebbd7",
                             "#71c7ec",
                             "#ffbaba", 
                             "#ff7b7b",
                             "#ff0000",
                             "#a70000",
                             "#efbbff", 
                             "#d896ff", 
                             "#be29ec", 
                             "#ffc2cd",
                             "#fff9ae",
                             "#ffd7b5",
                             "#ffb38a",
                             "#ff9248",
                             "#ff6700",
                             "#a3c1ad")) +
  theme_classic() +
  facet_wrap(~ host_taxon_name, scales="free") +
  theme(axis.title.x=element_text(size = 12, vjust = -1),
        axis.title.y = element_text(size = 12, vjust = +2),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(angle = 45, size =10, hjust = 1, vjust =1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'))+
  theme(strip.text.x = element_text(size = 10))+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Relative Abundance") +
  xlab("Month")

RA_plot_Order_SS

ggsave("SS_polyp_RA.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

#Creating NMDS plots to show beta diversity between species
polyp_RA = microbiome::transform(polyp_physeq, "compositional")
polyp_ord = ordinate(polyp_RA, "NMDS", "bray")

NMDS_colors <- brewer.pal(n=8, "Paired")
plot_ordination(polyp_physeq, polyp_ord, color="host_taxon_abbreviation") +
  geom_point(size = 2, )+
  labs(color = "Host Species")+
  scale_color_manual(values=NMDS_colors)+
  scale_shape_manual(values = c(1, 7, 8)) + 
  theme_classic()+
  theme(axis.title.x=element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

ggsave("NMDS_hosts-pres.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

polyp_host_pairwise <- pairwise.adonis(t(otu_table(polyp_ra)), sample_data(polyp_ra)$host_taxon_abbreviation, sim.method 	= "bray", p.adjust.m = "bonferroni") 
view(polyp_host_pairwise)
gridExtra::grid.table(polyp_host_pairwise)


#Testing reef significance
polyp_reef_pairwise <- pairwise.adonis(t(otu_table(polyp_ra)), sample_data(polyp_ra)$sample_site_id, sim.method 	= "bray", p.adjust.m = "bonferroni") 
polyp_key_pairwise <- pairwise.adonis(t(otu_table(polyp_ra)), sample_data(polyp_ra)$Key, sim.method 	= "bray", p.adjust.m = "bonferroni") 

gridExtra::grid.table(polyp_reef_pairwise)
gridExtra::grid.table(polyp_key_pairwise)


#Part 2: Beta Diversity Between Sample Types
rownames(Mbon_meta_all) <- Mbon_meta_all$SampleID 
physeq_16S <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                      tax_table(as.data.frame(taxtable_16S) %>%
                                  column_to_rownames("Feature.ID") %>%
                                  as.matrix()), sample_data(Mbon_meta_all))
physeq_16S = filter_taxa(physeq_16S, function(x) sum(x>10) > (0.30*length(x)),TRUE)

#NMDS Plots
RA_all = microbiome::transform(physeq_16S, "compositional")
all_ord = ordinate(RA_all, "NMDS", "bray")

plot_ordination(physeq_16S, all_ord, color="env_feature", shape = "host_taxon_abbreviation") +
  geom_point(size = 2.5)+
  labs(color = "Sample Type", shape = "Host Species")+
  scale_color_manual(values=c("#AE0F1B", "#E69F00", "#56B4E9"))+
  scale_shape_manual(values = c(1, 7, 8)) + 
  theme_classic()+
  theme(axis.title.x=element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(color=env_feature, fill = env_feature), show.legend=FALSE)+
  scale_fill_manual(values = c("#AE0F1B", "#E69F00", "#56B4E9"))

?stat_ellipse
  
ggsave("NMDS_sample_type_presentation2.jpeg", plot = last_plot(), device = "jpeg", 
       width = 6,
       height = 5,
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

all_clr = microbiome::transform(physeq_16S, "clr")
all_sampletype_pairwise <- pairwise.adonis(t(otu_table(RA_all)), sample_data(RA_all)$env_feature, sim.method 	= "bray", p.adjust.m = "bonferroni") 

gridExtra::grid.table(all_sampletype_pairwise)

#Part 4: Compute beta dissimilarity indices by sample

dis_cer <- vegdist(t(otu_table(RA_all)), method ="bray")
dim(dis_cer)
dim(sample_data(RA_all))

mod_cer <- betadisper(dis_cer, sample_data(RA_all)$host_taxon_abbreviation)
mod_dist_cer= as.data.frame(mod_cer$distances)
mod_dist_cer <- tibble::rownames_to_column(mod_dist_cer, "SampleID")

#Join to metadata file
Mbon_meta_all_with_beta <- inner_join(mod_dist_cer, Mbon_meta_all, by="SampleID")
Mbon_meta_all_with_beta <- Mbon_meta_all_with_beta %>% 
  rename("beta_distances" = "mod_cer$distances")


