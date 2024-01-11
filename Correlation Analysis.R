library(readr)
library(rerddap)
library(lubridate)
library(dplyr)
library(flexdashboard)
library(reshape2)
library(leaflet)
library(ggplot2)
library(vegan)
library(xts)
library(dygraphs)
library(plotly)
library(mapdata)
library(RColorBrewer)
library(tibble)
library(tidyr)
palette(brewer.pal(8, "Set2"))

#Read in lat/long information for sample sites and clean up dataframe
MBON_sampling_dates <- read.csv("~/Documents/University of Miami/MBON Project/metadata/MBON_months_coordinates_from Stephanie intern sept2022_CDS updated_9-22-2022.csv")
MBON_sampling_dates <- subset(MBON_sampling_dates, select = c("year", "month", "latitude", "longitude", "sample_site_id", "full.sample.date"))     
MBON_sampling_dates <- MBON_sampling_dates %>%
  rename(full_sample_date = full.sample.date)

## TEMPERATURE CORRELATION ##
## set dataset source (monthly SST)
SSTsource = info("jplMURSST41mday")
## Get sst for each site and clean up the data
MR <- griddap(SSTsource, 
               time=c("2016-11-14", "2018-02-16"),
               longitude = c(-80.3692, -80.3692),
               latitude = c(25.0083, 25.0083),
               fields = "sst",
               fmt = "csv")

MR <- MR %>%
  dplyr::mutate(year = lubridate::year(time), 
                month = lubridate::month(time), 
                day = lubridate::day(time))

MR <- MR %>% add_column(sample_site_id = NA)
MR$sample_site_id = "MR"

WS <- griddap(SSTsource, 
              time=c("2016-11-07", "2018-02-02"),
              longitude = c(-81.7080, -81.7080),
              latitude = c(24.4733, 24.4733),
              fields = "sst",
              fmt = "csv")

WS <- WS %>%
  dplyr::mutate(year = lubridate::year(time), 
                month = lubridate::month(time), 
                day = lubridate::day(time))

WS <- WS %>% add_column(sample_site_id = NA)
WS$sample_site_id = "WS"

TN <- griddap(SSTsource, 
              time=c("2016-11-14", "2018-02-01"),
              longitude = c(-80.7500, -80.7500),
              latitude = c(24.7667, 24.7667),
              fields = "sst",
              fmt = "csv")

TN <- TN %>%
  dplyr::mutate(year = lubridate::year(time), 
                month = lubridate::month(time), 
                day = lubridate::day(time))

TN <- TN %>% add_column(sample_site_id = NA)
TN$sample_site_id = "TN"


SR <- griddap(SSTsource, 
              time=c("2016-11-14", "2018-02-16"),
              longitude = c(-81.1106, -81.1106),
              latitude = c(24.6263, 24.6263),
              fields = "sst",
              fmt = "csv")

SR <- SR %>%
  dplyr::mutate(year = lubridate::year(time), 
                month = lubridate::month(time), 
                day = lubridate::day(time))

SR <- SR %>% add_column(sample_site_id = NA)
SR$sample_site_id = "SR"


LK <- griddap(SSTsource, 
              time=c("2016-11-14", "2018-02-02"),
              longitude = c(-81.4069, -81.4069),
              latitude = c(24.6618, 24.6618),
              fields = "sst",
              fmt = "csv")

LK <- LK %>%
  dplyr::mutate(year = lubridate::year(time), 
                month = lubridate::month(time), 
                day = lubridate::day(time))

LK <- LK %>% add_column(sample_site_id = NA)
LK$sample_site_id = "LK"

CR <- griddap(SSTsource, 
              time=c("2016-11-14", "2018-02-01"),
              longitude = c(-80.6100, -80.6100),
              latitude = c(24.9024, 24.9024),
              fields = "sst",
              fmt = "csv")

CR <- CR %>%
  dplyr::mutate(year = lubridate::year(time), 
                month = lubridate::month(time), 
                day = lubridate::day(time))

CR <- CR %>% add_column(sample_site_id = NA)
CR$sample_site_id = "CR"

#Combine all SST measurements into one
blended_sst <- rbind(MR, WS, CR, LK, SR, TN)

blended_sst <- blended_sst %>% mutate(month = recode(month, 
                                                     "1" = "January",
                                                     "2" = "February",
                                                     "3" = "March",
                                                     "4" = "April", 
                                                     "5" = "May",
                                                     "6" = "June",
                                                     "7" = "July",
                                                     "8" = "August",
                                                     "9" = "September",
                                                     "10" = "October",
                                                     "11" = "November",
                                                     "12" = "December"))

blended_sst <- blended_sst %>% mutate(sample_site_id = recode(sample_site_id, 
                                                                  "CR" = "Cheeca Rocks", 
                                                                  "LK" = "Looe Key", 
                                                                  "MR" = "Molasses Reef",
                                                                  "SR" = "Sombrero Reef",
                                                                  "TN" = "Tennessee Reef",
                                                                  "WS" = "Western Sambo",
                                                                  "ER" = "Emerald Reef"))

#Change month to a factor (instead of character)
blended_sst$month <- as.factor(blended_sst$month)
class(blended_sst$month)

blended_sst$year <- as.factor(blended_sst$year)
class(blended_sst$year)

#Repeat with metadata file (with alpha diversity values already added)
Mbon_meta_all_with_alpha$month <- as.factor(Mbon_meta_all_with_alpha$month)
class(Mbon_meta_all_with_alpha$month)

class(Mbon_meta_all_with_alpha$year)
Mbon_meta_all_with_alpha$year <- as.factor(Mbon_meta_all_with_alpha$year)

#Merge SST points to metadata
Mbon_meta_all_with_alpha <- left_join(Mbon_meta_all_with_alpha, blended_sst, by = c("sample_site_id", "year", "month"))
Mbon_meta_all_with_alpha = select(Mbon_meta_all_with_alpha, -c("latitude", "longitude", "time", "day.y"))

### PLOT ALPHA DIVERSITY AGAINST SST ###
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(ggpubr)
library(tidyr)

#Subset to give only polyps and plot sst to alpha
polyp_with_alpha <- subset(Mbon_meta_all_with_alpha, env_feature=="Polyp")
polyp_temp_cor <- polyp_with_alpha %>%
  drop_na(sst) %>%
  ggplot(aes(x = sst, y= Shannon, shape = host_taxon_abbreviation)) +
  geom_point(aes(color = month), alpha = 0.7, size = 1) +
  geom_smooth(method="lm", se=FALSE) + 
  theme_bw()+
  facet_wrap(~ env_feature, scales="free")+
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6)) +
  labs(x="Temperature (°C)")+
  ylab("Shannon Diversity Index")+
  scale_color_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"))+
  labs(color = "Month", shape = "Host Species")+
  ylim(0,5.9)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x = "right",
               label.y = "top",
               vstep = 0.05, size =2)
polyp_temp_cor

#Subset to give only water and plot sst to alpha
water_with_alpha <- subset(Mbon_meta_all_with_alpha, env_feature=="Water")
water_with_alpha$month <- factor(water_with_alpha$month, levels = c("February", "April", "July", "November", "December"))
levels(water_with_alpha$month)

#Correlate water sst with alpha diversity
water_temp_cor <- water_with_alpha %>%
  drop_na(sst) %>%
  ggplot(aes(x = sst, y= Shannon)) +
  geom_point(aes(color = month), alpha = 0.7, size = 1) +
  geom_smooth(method="lm", se=FALSE) + 
  theme_bw()+
  facet_wrap(~ env_feature, scales="free")+
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6)) +
  labs(x="Temperature")+
  scale_color_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"))+
  ylab("Shannon Diversity Index")+
  labs(color = "Month")+
  ylim(0,5.9)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x = "right",
               label.y = "top",
               vstep = 0.05, size = 2)
water_temp_cor

#Subset to give only sediment
sed_with_alpha <- subset(Mbon_meta_all_with_alpha, env_feature=="Sediment")
sed_with_alpha$month <- factor(sed_with_alpha$month, levels = c("February", "April", "July", "November", "December"))
levels(sed_with_alpha$month)

#Correlate sed sst with alpha diversity
sed_temp_cor <- sed_with_alpha %>%
  drop_na(sst) %>%
  ggplot(aes(x = sst, y= Shannon)) +
  geom_point(aes(color = month), alpha = 0.7, size = 1) +
  geom_smooth(method="lm", se=FALSE) + 
  theme_bw()+
  facet_wrap(~ env_feature, scales="free")+
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6))+
  labs(x="Temperature")+
  scale_color_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"))+
  ylab("Shannon Diversity Index")+
  labs(color = "Month")+
  ylim(0,5.9)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x = "right",
               label.y = "top",
               vstep = 0.05, size = 2)

sed_temp_cor

### CHL Analysis ###
library(readxl)
chl_monthly <- read.csv("~/Documents/University of Miami/MBON Project/Chlorophyll\ Data/chl_mir_monthly.csv")
time_vector <- read_excel("Downloads/time_vector.xlsx", 
                      col_names = FALSE)


#Convert first row to column names in time vector
install.packages("janitor")
library(janitor)
time_vector <- time_vector %>%
  row_to_names(row_number = 1)

#Round to 6 decimal places
class(time_vector$dec_yr)
class(chl_monthly$Dec_yr)

options(digits=10)
time_vector$dec_yr <- as.numeric(time_vector$dec_yr)
time_vector$dec_yr <- round(time_vector$dec_yr, digits = 6)

chl_monthly$Dec_yr <- as.numeric(chl_monthly$Dec_yr)

CR_chl <- chl_monthly %>% select(Dec_yr, CR, CR.1)

CR_chl <- CR_chl %>%
  rename("mean" = "CR",
         "sd" = "CR.1")
CR_chl = CR_chl[-1,]

ER_chl <- chl_monthly %>% select(Dec_yr, ER, ER.1)
ER_chl <- ER_chl %>%
  rename(mean = ER,
         sd = ER.1)
ER_chl = ER_chl[-1,]

LK_chl <- chl_monthly %>% select(Dec_yr, LK, LK.1)
LK_chl <- LK_chl %>%
  rename(mean = LK,
         sd = LK.1)
LK_chl = LK_chl[-1,]

MR_chl <- chl_monthly %>% select(Dec_yr, MR, MR.1)
MR_chl <- MR_chl %>%
  rename(mean = MR,
         sd = MR.1)
MR_chl = MR_chl[-1,]

SR_chl <- chl_monthly %>% select(Dec_yr, SR, SR.1)
SR_chl <- SR_chl %>%
  rename(mean = SR,
         sd = SR.1)
SR_chl = SR_chl[-1,]

TN_chl <- chl_monthly %>% select(Dec_yr, TN, TN.1)
TN_chl <- TN_chl %>%
  rename(mean = TN,
         sd = TN.1)
TN_chl = TN_chl[-1,]

WS_chl <- chl_monthly %>% select(Dec_yr, WS, WS.1)
WS_chl <- WS_chl %>%
  rename(mean = WS,
         sd = WS.1)
WS_chl = WS_chl[-1,]

#Rename first column in time_vector 
time_vector <- time_vector %>%
  rename(Dec_yr = dec_yr)

CR_chl <- left_join(CR_chl, time_vector, by = "Dec_yr")
ER_chl <- left_join(ER_chl, time_vector, by = "Dec_yr")
LK_chl <- left_join(LK_chl, time_vector, by = "Dec_yr")
MR_chl <- left_join(MR_chl, time_vector, by = "Dec_yr")
SR_chl <- left_join(SR_chl, time_vector, by = "Dec_yr")
TN_chl <- left_join(TN_chl, time_vector, by = "Dec_yr")
WS_chl <- left_join(WS_chl, time_vector, by = "Dec_yr")

CR_chl <- CR_chl %>% add_column(sample_site_id = NA)
CR_chl$sample_site_id = "CR"

ER_chl <- ER_chl %>% add_column(sample_site_id = NA)
ER_chl$sample_site_id = "ER"

LK_chl <- LK_chl %>% add_column(sample_site_id = NA)
LK_chl$sample_site_id = "LK"

MR_chl <- MR_chl %>% add_column(sample_site_id = NA)
MR_chl$sample_site_id = "MR"

SR_chl <- SR_chl %>% add_column(sample_site_id = NA)
SR_chl$sample_site_id = "SR"

TN_chl <- TN_chl %>% add_column(sample_site_id = NA)
TN_chl$sample_site_id = "TN"

WS_chl <- WS_chl %>% add_column(sample_site_id = NA)
WS_chl$sample_site_id = "WS"

blended_chl <- rbind(MR_chl, WS_chl, CR_chl, LK_chl, SR_chl, TN_chl, ER_chl)

#Rename columns to match with Mbon metadata file
blended_chl <- blended_chl %>%
  rename(year = yr,
         month = mo,
         chl_mean = mean,
         chl_sd = sd)

class(blended_chl$year)
blended_chl$year <- as.factor(blended_chl$year)

class(blended_chl$month)
blended_chl$month <- as.factor(blended_chl$month)

class(Mbon_meta_all_with_alpha$year)
Mbon_meta_all_with_alpha$year <- as.factor(Mbon_meta_all_with_alpha$year)

blended_chl <- blended_chl %>% mutate(month = recode(month,
                                                     "1" = "January",
                                                     "2" = "February", 
                                                     "3" = "March",
                                                     "4" = "April",
                                                     "5" = "May",
                                                     "6" = "June",
                                                     "7" = "July",
                                                     "8" = "August",
                                                     "9" = "September",
                                                     "10" = "October",
                                                     "11" = "November",
                                                     "12" = "December"))

blended_chl <- blended_chl %>% mutate(sample_site_id = recode(sample_site_id, 
                                                              "CR" = "Cheeca Rocks", 
                                                              "LK" = "Looe Key", 
                                                              "MR" = "Molasses Reef",
                                                              "SR" = "Sombrero Reef",
                                                              "TN" = "Tennessee Reef",
                                                              "WS" = "Western Sambo",
                                                              "ER" = "Emerald Reef"))
Mbon_meta_all_with_alpha <- left_join(Mbon_meta_all_with_alpha, blended_chl, by = c("sample_site_id", "year", "month"))
Mbon_meta_all_with_alpha = select(Mbon_meta_all_with_alpha, -c("Dec_yr"))

#Subset to give only polyps and plot chl to alpha
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(ggpubr)
library(scales)

polyp_with_alpha <- subset(Mbon_meta_all_with_alpha, env_feature=="Polyp")
polyp_chl_cor <- polyp_with_alpha %>%
  drop_na(chl_mean) %>%
  ggplot(aes(x = as.numeric(chl_mean), y= Shannon, shape = host_taxon_abbreviation)) +
  geom_point(aes(color = month), alpha = 0.7, size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  geom_smooth(method="lm", se=FALSE) + 
  theme_bw()+
  facet_wrap(~ env_feature, scales="free")+
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        axis.title.y=element_blank()) +
  labs(x="Chlorophyll (mg/L)")+
  scale_color_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"))+
  ylab("Shannon Diversity Index")+
  labs(color = "Month", shape = "Host Species")+
  ylim(0,5.9)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x = "right",
               label.y = "top",
               vstep = 0.05,
               size = 2)

polyp_chl_cor

class(Mbon_meta_all_with_alpha$Shannon)
Mbon_meta_all_with_alpha$chl_mean <- as.numeric(Mbon_meta_all_with_alpha$chl_mean)

sed_with_alpha <- subset(Mbon_meta_all_with_alpha, env_feature=="Sediment")
sed_chl_cor <- sed_with_alpha %>%
  drop_na(chl_mean) %>%
  ggplot(aes(x = as.numeric(chl_mean), y= Shannon)) +
  geom_point(aes(color = month), alpha = 0.7, size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  geom_smooth(method="lm", se=FALSE) + 
  theme_bw()+
  facet_wrap(~ env_feature, scales="free")+
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        axis.title.y=element_blank()) +
  labs(x="Chlorophyll (mg/L)")+
  scale_color_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"))+
  ylab("Shannon Diversity Index")+
  labs(color = "Month")+
  ylim(0,5.9)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x = "right",
               label.y = "top",
               vstep = 0.05,
               size = 2)

sed_chl_cor

water_with_alpha <- subset(Mbon_meta_all_with_alpha, env_feature=="Water")

water_chl_cor <- water_with_alpha %>%
  drop_na(chl_mean) %>%
  ggplot(aes(x = as.numeric(chl_mean), y= Shannon)) +
  geom_point(aes(color = month), alpha = 0.7, size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  geom_smooth(method="lm", se=FALSE) + 
  theme_bw()+
  facet_wrap(~ env_feature, scales="free")+
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        axis.title.y=element_blank()) +
  labs(x="Chlorophyll (mg/L)")+
  scale_color_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"))+
  ylab("Shannon Diversity Index")+
  labs(color = "Month")+
  ylim(0,5.9)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x = "right",
               label.y = "top",
               vstep = 0.05,
               size = 2)
water_chl_cor

#Display all three plots together:

cor_all <- ggarrange(polyp_temp_cor, polyp_chl_cor, water_temp_cor, water_chl_cor,  sed_temp_cor,  sed_chl_cor,
                     labels = c("A", "B", "C", "D", "E", "F"),
                     hjust = -0.2,
                     ncol = 2, nrow = 3)
cor_all

ggsave("~/Documents/University of Miami/MBON Project/figures/correlation_combined.tif",
       device = tiff,
       width = 7.5, 
       height = 8.75, 
       units = "in", 
       dpi=300)

#Making CCA plots to show effect of environmental factors

library("phyloseq")
library("ggplot2")
library("ape")
library("vegan")
library("goeveg")
library("dplyr")
library("plyr")
library(scales)
library(grid)
library(reshape2)
library(grid)
library("cowplot")
library(gridExtra)
library(ggpubr)

#Vegan package CCA
Mbon_meta_all_with_alpha <- tibble::column_to_rownames(Mbon_meta_all_with_alpha, "SampleID")

environment <- select(Mbon_meta_all_with_alpha, sst, chl_mean) #create dataframe with just sst values
view(environment)
environment <- environment %>% drop_na(sst) #drop NA values
environment <- environment %>% drop_na(chl_mean) #drop NA values

view(ASVtable_data)
count.data.t <- t(ASVtable_data)

dim(environment)
dim(count.data.t)

uneeded<-which(!rownames(count.data.t) %in% rownames(environment))    
count.data.t2<-count.data.t[-uneeded,]

dim(count.data.t2)
sst2<-data.matrix(environment, rownames.force = NA)
view(sst2)

cca <- cca(count.data.t2, sst2)
plot(cca)



