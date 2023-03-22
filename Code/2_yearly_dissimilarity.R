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

#######################
##DATA##
#######################

#Pull in cleaned EBS data from 1_PullEBS_data.R

EBS.r <- fread(file.path("Data","EBS.r.csv"))

#######################
##ANALYSIS
#######################

#check that there is one species record for each haul id
stopifnot(
  nrow(unique(EBS.r[,.(accepted_name, haul_id)])) == nrow(EBS.r)
          )

##Calculate dissimilarity values


#list years
EBS.r[, year:= as.numeric(year)] #make numeric
setorder(EBS.r,  year)
years <- unique(EBS.r[, year])

#haul id keys
haul_ids <- unique(EBS.r[, haul_id])
haul_ids_key <- data.table(haul_id = haul_ids,  key_ID = seq(1, length(haul_ids),  by = 1))


#convert haul_ids to numeric key_ids
EBS.r <- EBS.r[haul_ids_key,  on = "haul_id"]

#set up data table
distances_dissimilarities_allyears <- data.table("haul_id1" = integer(), 
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
  EBS.singleyear <- EBS.r[year == years[j], ]
  
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
    distances_dissimilarities_allyears <- rbind(distances_dissimilarities_allyears,  dissimilarities_full)

    
    print(paste0(j, "/", length(years)))
    
}

##Take average for each year
distances_dissimilarities_allyears[, bray_curtis_dissimilarity_balanced_mean := 
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

EBS.distances_dissimilarities_allyears.r <- unique(distances_dissimilarities_allyears[,.(year,bray_curtis_dissimilarity_balanced_mean, bray_curtis_dissimilarity_gradient_mean, bray_curtis_dissimilarity_total_mean, jaccard_dissimilarity_turnover_mean, jaccard_dissimilarity_nestedness_mean, jaccard_dissimilarity_total_mean)])

#####################
##Save
#####################
fwrite(EBS.distances_dissimilarities_allyears.r, file.path("Output","EBS.distances_dissimilarities_allyears.r.csv"))


