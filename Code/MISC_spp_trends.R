###########################
#This script merges the map figure with the trend figure for year~dissimialrity for each domain
#######################
##VERSIONS##
#######################
#R 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
#macOS Big Sur 11.7

#######################
##PACKAGES##
#######################

library(data.table)
library(ggplot2)
library(cowplot)

#######################
##DATA##
#######################
#Pull in cleaned EBS data from 1_PullEBS_data.R

EBS.domains.r.full <- fread(file.path("Data","EBS.domains.r.full.csv"))

######################
##ANALYSIS##
#######################

#how many species in each domain?
EBS.domains.r.full[,domain_richness := uniqueN(accepted_name),.(domain)]
species_domains <- unique(EBS.domains.r.full[,.(domain, domain_richness)])

#most abundant species (by #s and by biomass)
EBS.domains.r.full[,summed_abundance := sum(num_cpue, na.rm = T), accepted_name][,summed_biomass := sum(wgt_cpue, na.rm = T), accepted_name]
summed_spp_abundance <- unique(EBS.domains.r.full[,.(accepted_name, summed_abundance, summed_biomass)])

#top 5 species
#Limanda aspera
#Gadus chalcogrammus
#Hippoglossoides elassodon
#Pleuronectes quadrituberculatus
#Atheresthes stomias

setorder(summed_spp_abundance, -summed_biomass) #set order of summed_spp_abundance by biomass

top_spp <- summed_spp_abundance[1:5, accepted_name]

#data table
top_spp.dt <- EBS.domains.r.full[accepted_name %in% top_spp, ][,wgt_cpue_summed := sum(wgt_cpue, na.rm = T), .(accepted_name, year, domain)]
#unique values
top_spp.dt <- unique(top_spp.dt[,.(year, domain, accepted_name, wgt_cpue_summed)])

#top species

top_spp_overtime <- ggplot(data = top_spp.dt) +
  geom_line(aes(x = year, y = wgt_cpue_summed/1000, color = domain), size = 2) +
  facet_grid(accepted_name~domain, switch = "y") +
  scale_colour_manual(values = c("#AA4499","#44AA99","#999933")) +
  labs(x = "Year", y = "Summed CPUE (1000s of kg/km^2)") +
  theme_bw() +
  theme(legend.position = "null", text = element_text(size = 17))

ggsave(top_spp_overtime, path = file.path("Figures","Supplement"),
       filename = "top_spp_overtime.jpg", height = 16, width = 12, unit = "in")

top_spp_overtime_freey <- ggplot(data = top_spp.dt) +
  geom_line(aes(x = year, y = wgt_cpue_summed/1000, color = domain), size = 2) +
  facet_grid(accepted_name~domain, scales = "free", switch = "y") +
  scale_colour_manual(values = c("#AA4499","#44AA99","#999933")) +
  labs(x = "Year", y = "Summed CPUE (1000s of kg/km^2)") +
  theme_bw() +
  theme(legend.position = "null", text = element_text(size = 17))

ggsave(top_spp_overtime_freey, path = file.path("Figures","Supplement"),
       filename = "top_spp_overtime_freey.jpg", height = 16, width = 12, unit = "in")

#pollock only (Gadus chalcogrammus)
pollock <- EBS.domains.r.full[accepted_name == "Gadus chalcogrammus",][,wgt_cpue_summed := sum(wgt_cpue, na.rm = T),.(year, domain)]

pollock_byyear <- unique(pollock[,.(year, domain, wgt_cpue_summed)])

#pollock biomass overtime
pollock_overtime <- ggplot(data = pollock_byyear) +
  geom_smooth(aes(x = year, y = wgt_cpue_summed/1000)) +
  geom_point(aes(x = year, y = wgt_cpue_summed/1000, color = domain)) +
  facet_wrap(~domain) +
  scale_colour_manual(values = c("#AA4499","#44AA99","#999933")) +
  labs(x = "Year", y = "Summed CPUE (1000s of kg/km^2)", title = "Pollock") +
  theme_classic() +
  theme(legend.position = "null")

ggsave(pollock_overtime, path = file.path("Figures","Supplement"),
       filename = "pollock_overtime.jpg", height = 6, width = 8, unit = "in")

#all standardized, and in one plot
top_spp_overtime[,wgt_cpue_summed_scaled := scale(),.(domain, )]
