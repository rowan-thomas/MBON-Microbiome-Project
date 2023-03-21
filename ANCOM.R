### ANCOM Analysis ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")
library(ANCOMBC)

#ANCOM for bacteria driving differences between hosts (polyp only)
rownames(polyp) <- polyp$SampleID
polyp_physeq <-phyloseq(otu_table(ASVtable_data, taxa_are_rows= T),
                        tax_table(as.data.frame(taxtable_16S) %>%
                                    column_to_rownames("Feature.ID") %>%
                                    as.matrix()), sample_data(polyp))
polyp_physeq = filter_taxa(polyp_physeq, function(x) sum(x > 5) > (0.15*length(x)), TRUE)
polyp_ancombc = ancombc(phyloseq = polyp_physeq, 
                        formula = "host_taxon_abbreviation",
                        p_adj_method = "holm",
                        group = "host_taxon_abbreviation",
                        lib_cut = 1000,
                        struc_zero = TRUE, 
                        neg_lb = TRUE,
                        conserve = FALSE, 
                        alpha = 0.001, 
                        global = TRUE)
citation("ANCOMBC")

levels(Mbon_meta_all$host_taxon_abbreviation)
#AC - OF - SS

#Primary: individual comparisons, negatives mean significant in AC, positives mean significant in compared
polyp_res = polyp_ancombc$res
col_name = c("Feature.ID", "OF", "SS")

polyp_lfc = polyp_res$lfc
polyp_lfc = rownames_to_column(polyp_lfc, var="Feature.ID")
colnames(polyp_lfc) = col_name

#ASVs significant in acropora compared to orbicella
acropora_orbicella_lfc = select(polyp_lfc, -'SS')

#ASVs significant in acropora compared to siderastrea
acropora_siderastrea_lfc = select(polyp_lfc, -'OF')

### LFC Graphs ###
#Acropora vs. Orbicella
acropora_orbicella_lfc = left_join(acropora_orbicella_lfc, taxtable_16S, by = "Feature.ID")
acropora_orbicella_lfc = left_join(acropora_orbicella_lfc, seqs, by = "Feature.ID")

acropora_orbicella_lfc = acropora_orbicella_lfc %>% mutate(Host =
                     case_when(OF <= 0 ~ "Acropora", 
                               OF > 0 ~ "Orbicella"))

view(acropora_orbicella_lfc)

acropora_significant_compared_OF <- subset(acropora_orbicella_lfc, OF < 0)
AC_sig_OF_sd <- sd(acropora_significant_compared_OF$OF)
AC_sig_OF_mean <- mean(acropora_significant_compared_OF$OF)
AC_sig_OF_mean - AC_sig_OF_sd # -1.394006 cut off for outliers significant in acropora
2*1.394006

orbicella_significant_compared_AC <- subset(acropora_orbicella_lfc, OF > 0)
OF_sig_AC_sd <- sd(orbicella_significant_compared_AC$OF)
OF_sig_AC_mean <- mean(orbicella_significant_compared_AC$OF)
OF_sig_AC_mean + OF_sig_AC_sd # 1.238354 cut off for outliers significant in acropora
2*1.238354


ancom_lfc_AC_OF <- ggplot(acropora_orbicella_lfc, aes(x=OF, y = Order)) + 
  geom_point(aes(color=Host), size=1)+
  labs(x = "Log-Fold Changes", y = "ASV Order")+
  scale_color_manual(values = c("Acropora" = "#A6CEE3", "Orbicella" = "#1F78B4"))+
  geom_vline(xintercept = -2.788012, linetype="dotted", 
                  color = "black", size=1)+
  geom_vline(xintercept = 2.476708, linetype="dotted", 
             color = "black", size=1)+
  theme_bw() +
  theme(axis.text = element_text(size = 8))

ancom_lfc_AC_OF

ggsave("lfc_AC_OF.jpeg", plot = last_plot(), device = "jpeg", 
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

#Subset acropora_orbicella_lfc to get most significant ASVs
acropora_orbicella_lfc_most_sig <- subset(acropora_orbicella_lfc, OF <= -2.788012 | OF >= 2.476708)
2*1.238354
acropora_orbicella_lfc_most_sig <- left_join(acropora_orbicella_lfc_most_sig, seqs, by = 'Feature.ID')

#Acropora vs. Siderastrea
acropora_siderastrea_lfc = select(polyp_lfc, -'OF')
acropora_siderastrea_lfc = left_join(acropora_siderastrea_lfc, taxtable_16S, by = "Feature.ID")
acropora_siderastrea_lfc = acropora_siderastrea_lfc %>% mutate(Host =
                                   case_when(SS <= 0 ~ "Acropora", 
                                             SS > 0 ~ "Siderastrea"))
acropora_significant_compared_SS <- subset(acropora_siderastrea_lfc, SS < 0)
AC_sig_SS_sd <- sd(acropora_significant_compared_SS$SS)
AC_sig_SS_mean <- mean(acropora_significant_compared_SS$SS)
AC_sig_SS_mean - AC_sig_SS_sd # -1.50802 cut off for outliers significant in acropora
1.50802*2

acropora_significant_compared_SS_cutoff <- subset(acropora_significant_compared_SS, SS < -1.50802)

siderastrea_significant_compared_AC <- subset(acropora_siderastrea_lfc, SS > 0)
SS_sig_AC_sd <- sd(siderastrea_significant_compared_AC$SS)
SS_sig_AC_mean <- mean(siderastrea_significant_compared_AC$SS)
SS_sig_AC_mean + SS_sig_AC_sd # 1.387094 cut off for outliers significant in acropora
1.387094*2

siderastrea_significant_compared_AC_cutoff <- subset(siderastrea_significant_compared_AC, SS > 1.387094)

combined_cutoff_AC_SS <- rbind(acropora_significant_compared_SS_cutoff, siderastrea_significant_compared_AC_cutoff)

ancom_lfc_AC_SS <- ggplot(acropora_siderastrea_lfc, aes(x=SS, y = Order)) + 
  geom_point(aes(color=Host), size=1)+
  labs(x = "Log-Fold Changes", y = "ASV Order")+
  scale_color_manual(values = c("Siderastrea" = "#B2DF8A", "Acropora" = "#A6CEE3"))+
  geom_vline(xintercept = -3.01604, linetype="dotted", 
             color = "black", size=1)+
  geom_vline(xintercept = 2.774188, linetype="dotted", 
             color = "black", size=1)+
  theme_bw()+
  theme(axis.text = element_text(size = 8))

ancom_lfc_AC_SS

ggsave("lfc_AC_SS.jpeg", plot = last_plot(), device = "jpeg", 
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

ancom_lfc_AC_SS_cutoff <- ggplot(combined_cutoff_AC_SS, aes(x=SS, y = Order)) + 
  geom_point(aes(color=Host), size=3)+
  labs(x = "Log-Fold Changes", y = "ASV Order")+
  scale_color_manual(values = c("Significant in Siderastrea" = "#B2DF8A", "Significant in Acropora" = "#A6CEE3"))+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

ancom_lfc_AC_SS_cutoff

ggsave("lfc_AC_SS_cutoff.jpeg", plot = last_plot(), device = "jpeg", 
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

ggarrange(ancom_lfc_AC_OF, ancom_lfc_AC_SS, ncol=1, nrow=2)

ggsave("lfc_both.jpeg", plot = last_plot(), device = "jpeg", 
       path = "~/Documents/University of Miami/MBON Project/figures",
       dpi = 300)

#Subset acropora_siderastrea_lfc to get most significant ASVs
acropora_siderastrea_lfc_most_sig <- subset(acropora_siderastrea_lfc, SS <= -3.01604 | SS >= 2.774188)
view(acropora_siderastrea_lfc_most_sig)



