### ANCOM Analysis ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")
library(ANCOMBC)
library(phyloseq)
library(ggrepel)
sessionInfo()

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

levels(Mbon_meta_all$host_taxon_abbreviation)
#AC - OF - SS

#Primary: individual comparisons, negatives mean significant in AC, positives mean significant in compared
polyp_res = polyp_ancombc$res
col_name = c("Feature.ID", "OF", "SS")

polyp_lfc = polyp_res$lfc
polyp_lfc = as.data.frame(polyp_lfc)
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
                     case_when(OF <= 0 ~ "A. cervicornis", 
                               OF > 0 ~ "O. faveolata"))

view(acropora_orbicella_lfc)
acropora_orbicella_lfc %>% summarize(OF=n())
acropora_siderastrea_lfc %>% summarize(SS=n())

acropora_significant_compared_OF <- subset(acropora_orbicella_lfc, OF < 0)
AC_sig_OF_sd <- sd(acropora_significant_compared_OF$OF)
AC_sig_OF_mean <- mean(acropora_significant_compared_OF$OF)
AC_sig_OF_mean - AC_sig_OF_sd 
2*-1.394006

acropora_significant_compared_OF_cutoff <- subset(acropora_significant_compared_OF, OF < -2.788012)
view(acropora_significant_compared_OF)

orbicella_significant_compared_AC <- subset(acropora_orbicella_lfc, OF > 0)
OF_sig_AC_sd <- sd(orbicella_significant_compared_AC$OF)
OF_sig_AC_mean <- mean(orbicella_significant_compared_AC$OF)
OF_sig_AC_mean + OF_sig_AC_sd # 1.238354 cut off for outliers significant in acropora
2*1.238354

orbicella_significant_compared_AC_cutoff <- subset(orbicella_significant_compared_AC, OF > 2.476708)
view(orbicella_significant_compared_AC_cutoff)

ANCOM_AC_OF_cutoff <- rbind(acropora_significant_compared_OF_cutoff, orbicella_significant_compared_AC_cutoff)
ANCOM_AC_OF_cutoff <- ANCOM_AC_OF_cutoff %>% drop_na(Order)

ancom_lfc_AC_OF <- ggplot(ANCOM_AC_OF_cutoff, aes(x=OF, y = Order)) + 
  geom_point(aes(color=Host), size=1)+
  geom_vline(xintercept = -1.394006, linetype="dotted", 
             color = "black", linewidth=1)+
  geom_vline(xintercept = 1.238354, linetype="dotted", 
             color = "black", linewidth=1)+
  xlim(-9, 7)+
  geom_label_repel(
    label=ANCOM_AC_OF_cutoff$ASVID, 
    nudge_y = 0.2, nudge_x = 1.2,
    label.padding = 0.2,
    size=2
  )+
  labs(x = "Log-Fold Changes", y = "ASV Order")+
  scale_color_manual(values = c("A. cervicornis" = "#A6CEE3", "O. faveolata" = "#1F78B4"))+
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = "bottom")

ancom_lfc_AC_OF

#Acropora vs. Siderastrea
acropora_siderastrea_lfc = left_join(acropora_siderastrea_lfc, taxtable_16S, by = "Feature.ID")
acropora_siderastrea_lfc = left_join(acropora_siderastrea_lfc, seqs, by = "Feature.ID")
acropora_siderastrea_lfc = acropora_siderastrea_lfc %>% mutate(Host =
                                   case_when(SS <= 0 ~ "A. cervicornis", 
                                             SS > 0 ~ "S. siderea"))

view(acropora_siderastrea_lfc)

acropora_significant_compared_SS <- subset(acropora_siderastrea_lfc, SS < 0)
AC_sig_SS_sd <- sd(acropora_significant_compared_SS$SS)
AC_sig_SS_mean <- mean(acropora_significant_compared_SS$SS)
AC_sig_SS_mean - AC_sig_SS_sd # -1.50802 cut off for outliers significant in acropora
2*-1.50802

acropora_significant_compared_SS <- subset(acropora_siderastrea_lfc, SS < -3.01604)
view(acropora_significant_compared_SS)

siderastrea_significant_compared_AC <- subset(acropora_siderastrea_lfc, SS > 0)
SS_sig_AC_sd <- sd(siderastrea_significant_compared_AC$SS)
SS_sig_AC_mean <- mean(siderastrea_significant_compared_AC$SS)
SS_sig_AC_mean + SS_sig_AC_sd # 1.387094 cut off for outliers significant in acropora
2*1.387094

siderastrea_significant_compared_AC <- subset(acropora_siderastrea_lfc, SS > 2.774188)
view(siderastrea_significant_compared_AC)

ANCOM_AC_SS_cutoff <- rbind(acropora_significant_compared_SS, siderastrea_significant_compared_AC)
ANCOM_AC_SS_cutoff <- ANCOM_AC_SS_cutoff %>% drop_na(Order)

ancom_lfc_AC_SS <- ggplot(ANCOM_AC_SS_cutoff, aes(x=SS, y = Order)) + 
  geom_vline(xintercept = -1.50802, linetype="dotted", 
             color = "black", linewidth=1)+
  geom_vline(xintercept = 1.387094, linetype="dotted", 
             color = "black", linewidth=1)+
  xlim(-9, 9)+
  geom_point(aes(color=Host), size=1)+
  geom_label_repel(
    label=ANCOM_AC_SS_cutoff$ASVID,
    nudge_y=0.1, nudge_x = 1,
    label.padding = 0.2,
    size=2
  )+
  labs(x = "Log-Fold Changes", y = "ASV Order")+
  scale_color_manual(values = c("S. siderea" = "#B2DF8A", "A. cervicornis" = "#A6CEE3"))+
  theme_bw()+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = "bottom")


ancom_lfc_AC_SS

ancom_both <- ggarrange(ancom_lfc_AC_OF, ancom_lfc_AC_SS,
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1)
ancom_both


ggsave("~/Documents/University of Miami/MBON Project/figures/ancom_both.tif", 
       device = tiff, 
       width = 7.5,
       height = 3.5,
       units = "in",
       dpi = 300)


