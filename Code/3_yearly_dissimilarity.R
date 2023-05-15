###########################
#This script calculates annual average pairwise dissimilarity
#Specifically, it calculates Bray Curtis Balanced Variation (metric of dissimilarity)
#######################
##VERSIONS##
#######################
#R 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
#macOS Big Sur 11.7

#######################
##PACKAGES##
#######################

library(data.table)
library(geosphere)
library(betapart)
library(ggplot2)
library(cowplot)
library(RobustLinearReg)
library(mgcv)
#######################
##DATA##
#######################

#Pull in cleaned EBS data from 1_PullEBS_data.R

EBS.domains.r.full <- fread(file.path("Data","EBS.domains.r.full.csv"))

#######################
##ANALYSIS
#######################

#check that there is one species record for each haul id
stopifnot(
  nrow(unique(EBS.domains.r.full[,.(accepted_name, haul_id)])) == nrow(EBS.domains.r.full)
          )

#calculate avg # haul ids per year domain and across all domains
EBS.domains.r.full[, hauls_per_year := uniqueN(haul_id), .(year)][, avg_hauls_per_year := median(hauls_per_year)]
EBS.domains.r.full[, hauls_per_year_bydomain := uniqueN(haul_id), .(year, domain)][, avg_hauls_per_year_bydomain := median(hauls_per_year_bydomain), domain]


#list years
EBS.domains.r.full[, year:= as.numeric(year)] #make numeric
setorder(EBS.domains.r.full,  year)
years <- unique(EBS.domains.r.full[, year])

#haul id keys
haul_ids <- unique(EBS.domains.r.full[, haul_id])
haul_ids_key <- data.table(haul_id = haul_ids,  key_ID = seq(1, length(haul_ids),  by = 1))


#convert haul_ids to numeric key_ids
EBS.domains.r.full <- EBS.domains.r.full[haul_ids_key,  on = "haul_id"]

###################
###Calculate dissimilarity for full region
###################

#set up data table
distances_dissimilarities_allyears_fullEBS <- data.table("haul_id1" = integer(), 
                                                 "haul_id2" = integer(), 
                                                 "distance" = numeric(), 
                                                 "bray_curtis_dissimilarity_balanced" = numeric(), 
                                                 year = integer(),
                                                 "jaccard_dissimilarity_turnover" = numeric(), 
                                                 "jaccard_dissimilarity_nestedness" =  numeric(), 
                                                 "bray_curtis_dissimilarity_gradient" = as.numeric(), 
                                                 "jaccard_dissimilarity_total" = numeric(), 
                                                 "bray_curtis_dissimilarity_total" = numeric())

for (j in 1:length(years)) {
  EBS.singleyear <- EBS.domains.r.full[year == years[j], ]
  
  #distances among cells
  setorder(EBS.singleyear,  key_ID)
  
  latitude_longitude_haul_id <- unique(EBS.singleyear[, .(latitude, longitude, key_ID)])
  distances <- distm(latitude_longitude_haul_id[, .(longitude, latitude)])
  key_IDs <- latitude_longitude_haul_id[, key_ID]
  
  colnames(distances) <- rownames(distances) <- key_IDs
  
  #wide to long
  haul_id_distances.l <- reshape2::melt(distances, varnames = (c("haul_id1",  "haul_id2")),  value.name = "distance")
  
  #make into data table
  haul_id_distances.l <- data.table(haul_id_distances.l)
    
    #if some rows have wgt_cpue missing,  get rid of these rows
    EBS.singleyear <- EBS.singleyear[complete.cases(EBS.singleyear[, wgt_cpue]), ]
    
    #longitude to wide data for community matrix,  column names are cell then species
    EBS.singleyear_wide <- dcast(EBS.singleyear,  key_ID + year ~ accepted_name,  value.var = "wgt_cpue",  fun.aggregate = sum) 
    
    
    ncols <- ncol(EBS.singleyear_wide)
    communitymatrix <- EBS.singleyear_wide[, 3:ncols] #community matrix
    communitymatrix.occurence <- communitymatrix
    communitymatrix.occurence[communitymatrix.occurence > 0] <- 1
    
    #list of haul_id keys
    key_IDs_subset <- EBS.singleyear_wide[,key_ID]
    
    dissimilarities_abundance <- beta.pair.abund(communitymatrix,  index.family = "bray") #BC dissimilarity 
    dissimilarities_occurrence <- beta.pair(communitymatrix.occurence,  index.family = "jaccard") #jaccard dissimilarity
    
    #make into matrix
    dissimilarities_abundance_balanced.m <- as.matrix(dissimilarities_abundance$beta.bray.bal,  labels=TRUE) #bal = balanced
    dissimilarities_abundance_gradient.m <- as.matrix(dissimilarities_abundance$beta.bray.gra,  labels=TRUE) #gra = gradient
    dissimilarities_abundance_total.m <- as.matrix(dissimilarities_abundance$beta.bray,  labels=TRUE) #total
    
    dissimilarities_occurrence_turnover.m <- as.matrix(dissimilarities_occurrence$beta.jtu,  labels=TRUE) #jtu = turnover
    dissimilarities_occurrence_nestedness.m <- as.matrix(dissimilarities_occurrence$beta.jne,  labels=TRUE) #jne = nestedness
    dissimilarities_occurrence_total.m <- as.matrix(dissimilarities_occurrence$beta.jac,  labels=TRUE) #total
    
    colnames(dissimilarities_abundance_balanced.m) <- rownames(dissimilarities_abundance_balanced.m) <- key_IDs_subset
    colnames(dissimilarities_abundance_gradient.m) <- rownames(dissimilarities_abundance_gradient.m) <- key_IDs_subset
    colnames(dissimilarities_abundance_total.m) <- rownames(dissimilarities_abundance_total.m) <- key_IDs_subset
    colnames(dissimilarities_occurrence_turnover.m) <- rownames(dissimilarities_occurrence_turnover.m) <- key_IDs_subset
    colnames(dissimilarities_occurrence_nestedness.m) <- rownames(dissimilarities_occurrence_nestedness.m) <- key_IDs_subset
    colnames(dissimilarities_occurrence_total.m) <- rownames(dissimilarities_occurrence_total.m) <- key_IDs_subset
    
    #reshape dissimilarities
    dissimilarities_abundance_balanced.l <- reshape2::melt(dissimilarities_abundance_balanced.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_balanced")
    dissimilarities_abundance_gradient.l <- reshape2::melt(dissimilarities_abundance_gradient.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_gradient")
    dissimilarities_abundance_total.l <- reshape2::melt(dissimilarities_abundance_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_total")
    
    dissimilarities_occurrence_turnover.l <- reshape2::melt(dissimilarities_occurrence_turnover.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_turnover")
    dissimilarities_occurrence_nestedness.l <- reshape2::melt(dissimilarities_occurrence_nestedness.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_nestedness")
    dissimilarities_occurrence_total.l <- reshape2::melt(dissimilarities_occurrence_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_total")
    
    #and then to data table format
    dissimilarities_abundance_balanced.l <- data.table(dissimilarities_abundance_balanced.l)
    dissimilarities_abundance_gradient.l <- data.table(dissimilarities_abundance_gradient.l)
    dissimilarities_abundance_total.l <- data.table(dissimilarities_abundance_total.l)
    dissimilarities_occurrence_turnover.l <- data.table(dissimilarities_occurrence_turnover.l)
    dissimilarities_occurrence_nestedness.l <- data.table(dissimilarities_occurrence_nestedness.l)
    dissimilarities_occurrence_total.l <- data.table(dissimilarities_occurrence_total.l)
    
    #add year for these values
    dissimilarities_abundance_balanced.l[,  "year" := years[j]]
    dissimilarities_occurrence_turnover.l[,  "year" := years[j]]
    
    
    dissimilarities_abundance_gradient.l[,  "year" := years[j]]
    dissimilarities_occurrence_nestedness.l[,  "year" := years[j]]
    
    dissimilarities_abundance_total.l[,  "year" := years[j]]
    dissimilarities_occurrence_total.l[,  "year" := years[j]]
    
    #merge distance with dissimilarity for this year with both metrics of dissimilarity
    dissimilarities_full <- haul_id_distances.l[dissimilarities_abundance_balanced.l,  on = c("haul_id1",  "haul_id2")]
    dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_turnover.l,  on = c("haul_id1",  "haul_id2",  "year")]
    
    dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_total.l,  on = c("haul_id1",  "haul_id2",  "year")]
    
    dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_nestedness.l,  on = c("haul_id1",  "haul_id2",  "year")]
    
    dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_gradient.l,  on = c("haul_id1",  "haul_id2",  "year")]
    
    dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_total.l,  on = c("haul_id1",  "haul_id2",  "year")]
    
    #add to data table
    distances_dissimilarities_allyears_fullEBS <- rbind(distances_dissimilarities_allyears_fullEBS,  dissimilarities_full)

    
    print(paste0(j, "/", length(years)))
    
}

##Take average for each year
distances_dissimilarities_allyears_fullEBS[, bray_curtis_dissimilarity_balanced_mean := 
                                      mean(bray_curtis_dissimilarity_balanced, na.rm = T),.(year)][,
                                     bray_curtis_dissimilarity_gradient_mean := 
                                      mean(bray_curtis_dissimilarity_gradient, na.rm = T),.(year)][,
                                     bray_curtis_dissimilarity_total_mean := 
                                      mean(bray_curtis_dissimilarity_total, na.rm = T),.(year)][,
                                     jaccard_dissimilarity_turnover_mean :=
                                      mean(jaccard_dissimilarity_turnover, na.rm = T),.(year)][,
                                     jaccard_dissimilarity_nestedness_mean := 
                                      mean(jaccard_dissimilarity_nestedness, na.rm = T),.(year)][,
                                     jaccard_dissimilarity_total_mean := 
                                      mean(jaccard_dissimilarity_total, na.rm = T),.(year)]

EBS.distances_dissimilarities_allyears_fullEBS.r <- unique(distances_dissimilarities_allyears_fullEBS[,.(year,bray_curtis_dissimilarity_balanced_mean, bray_curtis_dissimilarity_gradient_mean, bray_curtis_dissimilarity_total_mean, jaccard_dissimilarity_turnover_mean, jaccard_dissimilarity_nestedness_mean, jaccard_dissimilarity_total_mean)])

EBS.distances_dissimilarities_allyears_fullEBS.r[,domain := "Full"]


###################
###Calculate dissimilarity for each domain independently
###################
domains <- levels(as.factor(EBS.domains.r.full[,domain]))

#set up data table
distances_dissimilarities_allyears_by_domain <- data.table("haul_id1" = integer(), 
                                                         "haul_id2" = integer(), 
                                                         "distance" = numeric(), 
                                                         "bray_curtis_dissimilarity_balanced" = numeric(), 
                                                         year = integer(),
                                                         "jaccard_dissimilarity_turnover" = numeric(), 
                                                         "jaccard_dissimilarity_nestedness" =  numeric(), 
                                                         "bray_curtis_dissimilarity_gradient" = as.numeric(), 
                                                         "jaccard_dissimilarity_total" = numeric(), 
                                                         "bray_curtis_dissimilarity_total" = numeric(),
                                                         "domain" = as.character())

for (i in 1:length(domains)) {
  EBS.single.domain <- EBS.domains.r.full[domain == domains[i], ]
    for (j in 1:length(years)) {
      EBS.singleyear <- EBS.single.domain[year == years[j], ]
      
      #distances among cells
      setorder(EBS.singleyear,  key_ID)
      
      latitude_longitude_haul_id <- unique(EBS.singleyear[, .(latitude, longitude, key_ID)])
      distances <- distm(latitude_longitude_haul_id[, .(longitude, latitude)])
      key_IDs <- latitude_longitude_haul_id[, key_ID]
      
      colnames(distances) <- rownames(distances) <- key_IDs
      
      #wide to long
      haul_id_distances.l <- reshape2::melt(distances, varnames = (c("haul_id1",  "haul_id2")),  value.name = "distance")
      
      #make into data table
      haul_id_distances.l <- data.table(haul_id_distances.l)
      
      #if some rows have wgt_cpue missing,  get rid of these rows
      EBS.singleyear <- EBS.singleyear[complete.cases(EBS.singleyear[, wgt_cpue]), ]
      
      #longitude to wide data for community matrix,  column names are cell then species
      EBS.singleyear_wide <- dcast(EBS.singleyear,  key_ID + year ~ accepted_name,  value.var = "wgt_cpue",  fun.aggregate = sum) 
      
      
      ncols <- ncol(EBS.singleyear_wide)
      communitymatrix <- EBS.singleyear_wide[, 3:ncols] #community matrix
      communitymatrix.occurence <- communitymatrix
      communitymatrix.occurence[communitymatrix.occurence > 0] <- 1
      
      #list of haul_id keys
      key_IDs_subset <- EBS.singleyear_wide[,key_ID]
      
      dissimilarities_abundance <- beta.pair.abund(communitymatrix,  index.family = "bray") #BC dissimilarity 
      dissimilarities_occurrence <- beta.pair(communitymatrix.occurence,  index.family = "jaccard") #jaccard dissimilarity
      
      #make into matrix
      dissimilarities_abundance_balanced.m <- as.matrix(dissimilarities_abundance$beta.bray.bal,  labels=TRUE) #bal = balanced
      dissimilarities_abundance_gradient.m <- as.matrix(dissimilarities_abundance$beta.bray.gra,  labels=TRUE) #gra = gradient
      dissimilarities_abundance_total.m <- as.matrix(dissimilarities_abundance$beta.bray,  labels=TRUE) #total
      
      dissimilarities_occurrence_turnover.m <- as.matrix(dissimilarities_occurrence$beta.jtu,  labels=TRUE) #jtu = turnover
      dissimilarities_occurrence_nestedness.m <- as.matrix(dissimilarities_occurrence$beta.jne,  labels=TRUE) #jne = nestedness
      dissimilarities_occurrence_total.m <- as.matrix(dissimilarities_occurrence$beta.jac,  labels=TRUE) #total
      
      colnames(dissimilarities_abundance_balanced.m) <- rownames(dissimilarities_abundance_balanced.m) <- key_IDs_subset
      colnames(dissimilarities_abundance_gradient.m) <- rownames(dissimilarities_abundance_gradient.m) <- key_IDs_subset
      colnames(dissimilarities_abundance_total.m) <- rownames(dissimilarities_abundance_total.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_turnover.m) <- rownames(dissimilarities_occurrence_turnover.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_nestedness.m) <- rownames(dissimilarities_occurrence_nestedness.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_total.m) <- rownames(dissimilarities_occurrence_total.m) <- key_IDs_subset
      
      #reshape dissimilarities
      dissimilarities_abundance_balanced.l <- reshape2::melt(dissimilarities_abundance_balanced.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_balanced")
      dissimilarities_abundance_gradient.l <- reshape2::melt(dissimilarities_abundance_gradient.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_gradient")
      dissimilarities_abundance_total.l <- reshape2::melt(dissimilarities_abundance_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_total")
      
      dissimilarities_occurrence_turnover.l <- reshape2::melt(dissimilarities_occurrence_turnover.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_turnover")
      dissimilarities_occurrence_nestedness.l <- reshape2::melt(dissimilarities_occurrence_nestedness.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_nestedness")
      dissimilarities_occurrence_total.l <- reshape2::melt(dissimilarities_occurrence_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_total")
      
      #and then to data table format
      dissimilarities_abundance_balanced.l <- data.table(dissimilarities_abundance_balanced.l)
      dissimilarities_abundance_gradient.l <- data.table(dissimilarities_abundance_gradient.l)
      dissimilarities_abundance_total.l <- data.table(dissimilarities_abundance_total.l)
      dissimilarities_occurrence_turnover.l <- data.table(dissimilarities_occurrence_turnover.l)
      dissimilarities_occurrence_nestedness.l <- data.table(dissimilarities_occurrence_nestedness.l)
      dissimilarities_occurrence_total.l <- data.table(dissimilarities_occurrence_total.l)
      
      #add year for these values
      dissimilarities_abundance_balanced.l[,  "year" := years[j]]
      dissimilarities_occurrence_turnover.l[,  "year" := years[j]]
      
      
      dissimilarities_abundance_gradient.l[,  "year" := years[j]]
      dissimilarities_occurrence_nestedness.l[,  "year" := years[j]]
      
      dissimilarities_abundance_total.l[,  "year" := years[j]]
      dissimilarities_occurrence_total.l[,  "year" := years[j]]
      
      #merge distance with dissimilarity for this year with both metrics of dissimilarity
      dissimilarities_full <- haul_id_distances.l[dissimilarities_abundance_balanced.l,  on = c("haul_id1",  "haul_id2")]
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_turnover.l,  on = c("haul_id1",  "haul_id2",  "year")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_total.l,  on = c("haul_id1",  "haul_id2",  "year")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_nestedness.l,  on = c("haul_id1",  "haul_id2",  "year")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_gradient.l,  on = c("haul_id1",  "haul_id2",  "year")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_total.l,  on = c("haul_id1",  "haul_id2",  "year")]
      
      dissimilarities_full[,"domain" := domains[i]]
      
      #add to data table
      distances_dissimilarities_allyears_by_domain <- rbind(distances_dissimilarities_allyears_by_domain,  dissimilarities_full)
      
      
      print(paste0(j, "/", length(years)))
    }  
}

##Take average for each year and domain
distances_dissimilarities_allyears_by_domain[, bray_curtis_dissimilarity_balanced_mean := 
                                             mean(bray_curtis_dissimilarity_balanced, na.rm = T),.(year, domain)][,
                                                                                                          bray_curtis_dissimilarity_gradient_mean := 
                                                                                                            mean(bray_curtis_dissimilarity_gradient, na.rm = T),.(year, domain)][,
                                                                                                                                                                         bray_curtis_dissimilarity_total_mean := 
                                                                                                                                                                           mean(bray_curtis_dissimilarity_total, na.rm = T),.(year, domain)][,
                                                                                                                                                                                                                                     jaccard_dissimilarity_turnover_mean :=
                                                                                                                                                                                                                                       mean(jaccard_dissimilarity_turnover, na.rm = T),.(year, domain)][,
                                                                                                                                                                                                                                                                                                jaccard_dissimilarity_nestedness_mean := 
                                                                                                                                                                                                                                                                                                  mean(jaccard_dissimilarity_nestedness, na.rm = T),.(year, domain)][,
                                                                                                                                                                                                                                                                                                                                                             jaccard_dissimilarity_total_mean := 
                                                                                                                                                                                                                                                                                                                                                               mean(jaccard_dissimilarity_total, na.rm = T),.(year, domain)]

EBS.distances_dissimilarities_allyears_by_domain.r <- unique(distances_dissimilarities_allyears_by_domain[,.(year,bray_curtis_dissimilarity_balanced_mean, bray_curtis_dissimilarity_gradient_mean, bray_curtis_dissimilarity_total_mean, jaccard_dissimilarity_turnover_mean, jaccard_dissimilarity_nestedness_mean, jaccard_dissimilarity_total_mean, domain)])

EBS.distances_dissimilarities_allyears <- rbind(EBS.distances_dissimilarities_allyears_fullEBS.r, EBS.distances_dissimilarities_allyears_by_domain.r)

#GAM for plot

bray_curtis_balanced_gam_full_mod <- gam(bray_curtis_dissimilarity_balanced_mean~s(year, bs = "cr"), data = EBS.distances_dissimilarities_allyears[domain == "Full",])
bray_curtis_balanced_gam_inner_mod <- gam(bray_curtis_dissimilarity_balanced_mean~s(year, bs = "cr"), data = EBS.distances_dissimilarities_allyears[domain == "Inner",])
bray_curtis_balanced_gam_middle_mod <- gam(bray_curtis_dissimilarity_balanced_mean~s(year, bs = "cr"), data = EBS.distances_dissimilarities_allyears[domain == "Middle",])
bray_curtis_balanced_gam_outer_mod <- gam(bray_curtis_dissimilarity_balanced_mean~s(year, bs = "cr"), data = EBS.distances_dissimilarities_allyears[domain == "Outer",])

#linear models for comparison
bray_curtis_balanced_lm_full_mod <- lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",])
bray_curtis_balanced_lm_inner_mod <- lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",])
bray_curtis_balanced_lm_middle_mod <- lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",])
bray_curtis_balanced_lm_outer_mod <- lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",])

#GAM always lower
MuMIn::AICc(bray_curtis_balanced_gam_full_mod)- MuMIn::AICc(bray_curtis_balanced_lm_full_mod)
MuMIn::AICc(bray_curtis_balanced_gam_inner_mod)- MuMIn::AICc(bray_curtis_balanced_lm_inner_mod)
MuMIn::AICc(bray_curtis_balanced_gam_middle_mod)- MuMIn::AICc(bray_curtis_balanced_lm_middle_mod)
MuMIn::AICc(bray_curtis_balanced_gam_outer_mod)- MuMIn::AICc(bray_curtis_balanced_lm_outer_mod)

#add predicted values
EBS.distances_dissimilarities_allyears[domain == "Full", b_divers_gam_pred := predict(bray_curtis_balanced_gam_full_mod)]
EBS.distances_dissimilarities_allyears[domain == "Inner", b_divers_gam_pred := predict(bray_curtis_balanced_gam_inner_mod)]
EBS.distances_dissimilarities_allyears[domain == "Outer", b_divers_gam_pred := predict(bray_curtis_balanced_gam_outer_mod)]
EBS.distances_dissimilarities_allyears[domain == "Middle", b_divers_gam_pred := predict(bray_curtis_balanced_gam_middle_mod)]


###Visualize change over time (Bray Curtis Balanced)
ggplot(EBS.distances_dissimilarities_allyears) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean, color = domain)) +
  labs(color = "Domain", x = "Year", y = "β diversity") +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  theme_classic()
  
year_beta_bydomain <- ggplot(EBS.distances_dissimilarities_allyears) +
    geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean, color = domain)) +
    labs(color = "Domain", x = "Year", y = "β diversity") +
    geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean), method = "lm", se = F, color = "red") +
    geom_line(aes(x = year, y = b_divers_gam_pred), color = "darkgrey", linewidth = 1) +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
    facet_wrap(~domain) +
    theme_classic() +
  theme(legend.position = "null")

saveRDS(year_beta_bydomain, file.path("Figures","year_beta_bydomain.Rds"))



#slopes
slope_full <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_inner <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_middle <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_outer <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)

slope_full
slope_inner
slope_middle
slope_outer

#slopes
slope_se_full <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,2],2)
slope_se_inner <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,2],2)
slope_se_middle <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,2],2)
slope_se_outer <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,2],2)

slope_se_full
slope_se_inner
slope_se_middle
slope_se_outer

#R^2 values summary(model)
R2_full <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_inner <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_middle <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_outer <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)

R2_full
R2_inner
R2_middle
R2_outer

#p_values summary(model)
p_value_full <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_inner <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_middle <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_outer <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)

p_value_full
p_value_inner
p_value_middle
p_value_outer


#place model values on plot
year_beta_bydomain_annotate <- ggdraw(xlim = c(0,10), ylim = c(0,8)) +
  draw_plot(year_beta_bydomain, x = 0, y = 0, width = 10, height = 8) +
  geom_text(aes(x = 2.5, y = 0.8, label = paste0("Slope = ",slope_middle, "  R^2 = ",R2_middle, "  p = ",p_value_middle))) +
  geom_text(aes(x = 6.8, y = 0.8, label = paste0("Slope = ",slope_outer, "  R^2 = ",R2_outer, "  p = ",p_value_outer))) +
  geom_text(aes(x = 2.5, y = 4.5, label = paste0("Slope = ",slope_full, "  R^2 = ",R2_full, "  p = ",p_value_full))) +
  geom_text(aes(x = 6.5, y = 7.4, label = paste0("Slope = ",slope_inner, "  R^2 = ",R2_inner, "  p = ",p_value_inner)))

ggsave(year_beta_bydomain_annotate, path = file.path("Figures"), filename = "year_beta_bydomain_annotate.jpg", height = 6, width = 8, unit = "in")

#####################
#######SUPPLEMENT
####################

bray_curtis_balanced_theilsen_full_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",])
bray_curtis_balanced_theilsen_inner_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",])
bray_curtis_balanced_theilsen_middle_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",])
bray_curtis_balanced_theilsen_outer_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",])

EBS.distances_dissimilarities_allyears[domain == "Full", b_divers_theilsen_pred := predict(bray_curtis_balanced_theilsen_full_mod)]
EBS.distances_dissimilarities_allyears[domain == "Inner", b_divers_theilsen_pred := predict(bray_curtis_balanced_theilsen_inner_mod)]
EBS.distances_dissimilarities_allyears[domain == "Outer", b_divers_theilsen_pred := predict(bray_curtis_balanced_theilsen_outer_mod)]
EBS.distances_dissimilarities_allyears[domain == "Middle", b_divers_theilsen_pred := predict(bray_curtis_balanced_theilsen_middle_mod)]


###Visualize change over time (THEILSEN)
ggplot(EBS.distances_dissimilarities_allyears) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean, color = domain)) +
  labs(color = "Domain", x = "Year", y = "β diversity") +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  theme_classic()

year_beta_bray_curtis_balanced_theilsenreg_bydomain <- ggplot(EBS.distances_dissimilarities_allyears) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean, color = domain)) +
  labs(color = "Domain", x = "Year", y = "β diversity (Bray Curtis balanced)") +
  geom_line(aes(x = year, y = b_divers_theilsen_pred), color = "red") +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  facet_wrap(~domain) +
  theme_classic() +
  theme(legend.position = "null")

saveRDS(year_beta_bray_curtis_balanced_theilsenreg_bydomain, file.path("Figures","Supplement","Theil_Sen","year_beta_bray_curtis_balanced_theilsenreg_bydomain.Rds"))

#slopes
slope_bray_curtis_balanced_theilsen_full <- round(coefficients(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_bray_curtis_balanced_theilsen_inner <- round(coefficients(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_bray_curtis_balanced_theilsen_middle <- round(coefficients(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_bray_curtis_balanced_theilsen_outer <- round(coefficients(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)

slope_bray_curtis_balanced_theilsen_full
slope_bray_curtis_balanced_theilsen_inner
slope_bray_curtis_balanced_theilsen_middle
slope_bray_curtis_balanced_theilsen_outer

#R^2 values summary(model)
R2_bray_curtis_balanced_theilsen_full <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_bray_curtis_balanced_theilsen_inner <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_bray_curtis_balanced_theilsen_middle <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_bray_curtis_balanced_theilsen_outer <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)

R2_bray_curtis_balanced_theilsen_full
R2_bray_curtis_balanced_theilsen_inner
R2_bray_curtis_balanced_theilsen_middle
R2_bray_curtis_balanced_theilsen_outer

#p_values summary(model)
p_value_bray_curtis_balanced_theilsen_full <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_bray_curtis_balanced_theilsen_inner <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_bray_curtis_balanced_theilsen_middle <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_bray_curtis_balanced_theilsen_outer <- round(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)

p_value_bray_curtis_balanced_theilsen_full
p_value_bray_curtis_balanced_theilsen_inner
p_value_bray_curtis_balanced_theilsen_middle
p_value_bray_curtis_balanced_theilsen_outer


#place model values on plot
year_beta_bray_curtis_balanced_theilsenreg_bydomain_annotate <- ggdraw(xlim = c(0,10), ylim = c(0,8)) +
  draw_plot(year_beta_bray_curtis_balanced_theilsenreg_bydomain, x = 0, y = 0, width = 10, height = 8) +
  geom_text(aes(x = 2.5, y = 0.8, label = paste0("Slope = ",slope_bray_curtis_balanced_theilsen_middle, "  R^2 = ",R2_bray_curtis_balanced_theilsen_middle, "  p = ",p_value_bray_curtis_balanced_theilsen_middle))) +
  geom_text(aes(x = 6.8, y = 0.8, label = paste0("Slope = ",slope_bray_curtis_balanced_theilsen_outer, "  R^2 = ",R2_bray_curtis_balanced_theilsen_outer, "  p = ",p_value_bray_curtis_balanced_theilsen_outer))) +
  geom_text(aes(x = 2.5, y = 4.5, label = paste0("Slope = ",slope_bray_curtis_balanced_theilsen_full, "  R^2 = ",R2_bray_curtis_balanced_theilsen_full, "  p = ",p_value_bray_curtis_balanced_theilsen_full))) +
  geom_text(aes(x = 6.9, y = 7.4, label = paste0("Slope = ",slope_bray_curtis_balanced_theilsen_inner, "  R^2 = ",R2_bray_curtis_balanced_theilsen_inner, "  p = ",p_value_bray_curtis_balanced_theilsen_inner)))

ggsave(year_beta_bray_curtis_balanced_theilsenreg_bydomain_annotate,
       path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "year_beta_bray_curtis_balanced_theilsenreg_bydomain_annotate.jpg", height = 6, width = 8, unit = "in")



###Visualize change over time (Biomass.gradient BC)
ggplot(EBS.distances_dissimilarities_allyears) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_gradient_mean, color = domain)) +
  labs(color = "Domain", x = "Year", y = "β diversity") +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  theme_classic()

year_beta_bray_curtis_biomassgradient_bydomain <- ggplot(EBS.distances_dissimilarities_allyears) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_gradient_mean, color = domain)) +
  labs(color = "Domain", x = "Year", y = "β diversity (Bray Curtis biomass gradient)") +
  geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), method = "lm", se = F, color = "red") +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  facet_wrap(~domain) +
  theme_classic() +
  theme(legend.position = "null")

saveRDS(year_beta_bray_curtis_biomassgradient_bydomain, file.path("Figures","Supplement","Biomass_gradient","year_beta_bray_curtis_biomassgradient_bydomain.Rds"))

#slopes
slope_bray_curtis_biomass_gradient_full <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_bray_curtis_biomass_gradient_inner <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_bray_curtis_biomass_gradient_middle <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_bray_curtis_biomass_gradient_outer <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)

slope_bray_curtis_biomass_gradient_full
slope_bray_curtis_biomass_gradient_inner
slope_bray_curtis_biomass_gradient_middle
slope_bray_curtis_biomass_gradient_outer

#R^2 values summary(model)
R2_bray_curtis_biomass_gradient_full <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_bray_curtis_biomass_gradient_inner <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_bray_curtis_biomass_gradient_middle <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_bray_curtis_biomass_gradient_outer <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)

R2_bray_curtis_biomass_gradient_full
R2_bray_curtis_biomass_gradient_inner
R2_bray_curtis_biomass_gradient_middle
R2_bray_curtis_biomass_gradient_outer

#p_values summary(model)
p_value_bray_curtis_biomass_gradient_full <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_bray_curtis_biomass_gradient_inner <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_bray_curtis_biomass_gradient_middle <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_bray_curtis_biomass_gradient_outer <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)

p_value_bray_curtis_biomass_gradient_full
p_value_bray_curtis_biomass_gradient_inner
p_value_bray_curtis_biomass_gradient_middle
p_value_bray_curtis_biomass_gradient_outer


#place model values on plot
year_beta_bray_curtis_biomassgradient_bydomain_annotate <- ggdraw(xlim = c(0,10), ylim = c(0,8)) +
  draw_plot(year_beta_bray_curtis_biomassgradient_bydomain, x = 0, y = 0, width = 10, height = 8) +
  geom_text(aes(x = 2.5, y = 0.8, label = paste0("Slope = ",slope_bray_curtis_biomass_gradient_middle, "  R^2 = ",R2_bray_curtis_biomass_gradient_middle, "  p = ",p_value_bray_curtis_biomass_gradient_middle))) +
  geom_text(aes(x = 7.1, y = 0.8, label = paste0("Slope = ",slope_bray_curtis_biomass_gradient_outer, "  R^2 = ",R2_bray_curtis_biomass_gradient_outer, "  p = ",p_value_bray_curtis_biomass_gradient_outer))) +
  geom_text(aes(x = 2.5, y = 4.5, label = paste0("Slope = ",slope_bray_curtis_biomass_gradient_full, "  R^2 = ",R2_bray_curtis_biomass_gradient_full, "  p = ",p_value_bray_curtis_biomass_gradient_full))) +
  geom_text(aes(x = 7.1, y = 7.4, label = paste0("Slope = ",slope_bray_curtis_biomass_gradient_inner, "  R^2 = ",R2_bray_curtis_biomass_gradient_inner, "  p = ",p_value_bray_curtis_biomass_gradient_inner)))

ggsave(year_beta_bray_curtis_biomassgradient_bydomain_annotate, path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "year_beta_bray_curtis_biomassgradient_bydomain_annotate.jpg", height = 6, width = 8, unit = "in")



###Visualize change over time (Jaccard)
ggplot(EBS.distances_dissimilarities_allyears) +
  geom_point(aes(x = year, y = jaccard_dissimilarity_turnover_mean, color = domain)) +
  labs(color = "Domain", x = "Year", y = "β diversity") +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  theme_classic()

year_beta_jaccard_bydomain <- ggplot(EBS.distances_dissimilarities_allyears) +
  geom_point(aes(x = year, y = jaccard_dissimilarity_turnover_mean, color = domain)) +
  labs(color = "Domain", x = "Year", y = "β diversity (Jaccard)") +
  geom_smooth(aes(x = year, y = jaccard_dissimilarity_turnover_mean), method = "lm", se = F, color = "red") +
  scale_color_manual(values = c("black", "#AA4499","#44AA99","#999933"), labels = c("Full EBS","Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  facet_wrap(~domain) +
  theme_classic() +
  theme(legend.position = "null")

saveRDS(year_beta_jaccard_bydomain, file.path("Figures","Supplement","Jaccard","year_beta_jaccard_bydomain.Rds"))

#slopes
slope_jaccard_full <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_jaccard_inner <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_jaccard_middle <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_jaccard_outer <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)

slope_jaccard_full
slope_jaccard_inner
slope_jaccard_middle
slope_jaccard_outer

#R^2 values summary(model)
R2_jaccard_full <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_jaccard_inner <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_jaccard_middle <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_jaccard_outer <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)

R2_jaccard_full
R2_jaccard_inner
R2_jaccard_middle
R2_jaccard_outer

#p_values summary(model)
p_value_jaccard_full <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_jaccard_inner <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_jaccard_middle <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_jaccard_outer <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)

p_value_jaccard_full
p_value_jaccard_inner
p_value_jaccard_middle
p_value_jaccard_outer


#place model values on plot
year_beta_jaccard_bydomain_annotate <- ggdraw(xlim = c(0,10), ylim = c(0,8)) +
  draw_plot(year_beta_jaccard_bydomain, x = 0, y = 0, width = 10, height = 8) +
  geom_text(aes(x = 2.5, y = 0.8, label = paste0("Slope = ",slope_jaccard_middle, "  R^2 = ",R2_jaccard_middle, "  p = ",p_value_jaccard_middle))) +
  geom_text(aes(x = 7.1, y = 0.8, label = paste0("Slope = ",slope_jaccard_outer, "  R^2 = ",R2_jaccard_outer, "  p = ",p_value_jaccard_outer))) +
  geom_text(aes(x = 2.5, y = 4.5, label = paste0("Slope = ",slope_jaccard_full, "  R^2 = ",R2_jaccard_full, "  p = ",p_value_jaccard_full))) +
  geom_text(aes(x = 7.1, y = 7.4, label = paste0("Slope = ",slope_jaccard_inner, "  R^2 = ",R2_jaccard_inner, "  p = ",p_value_jaccard_inner)))

ggsave(year_beta_jaccard_bydomain_annotate, path = file.path("Figures","Supplement","Jaccard"), filename = "year_beta_jaccard_bydomain_annotate.jpg", height = 6, width = 8, unit = "in")


#####################
##Save
#####################
fwrite(EBS.distances_dissimilarities_allyears, file.path("Output","EBS.distances_dissimilarities_allyears.csv"))
