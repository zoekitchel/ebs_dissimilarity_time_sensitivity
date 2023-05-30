###########################
#Subsample function adapted by Zoë Kitchel from Easton White for Eastern Bering Sea β diversity (dissimilarity) Through Time

#######################
##VERSIONS##
#R 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
#macOS Big Sur 11.7

#######################
##PACKAGES##
#######################
library(data.table)
library(ggplot2)
library(animation)
library(RobustLinearReg)
library(cowplot)
library(purrr)
library(dplyr)
library(ggridges)
#######################
##FUNCTION
#######################
source(file.path("Code","FUN_linefit.R"))

#######################
##DATA
#######################
EBS.distances_dissimilarities_allyears <- fread(file.path("Output","EBS.distances_dissimilarities_allyears.csv"))
EBS.distances_dissimilarities_allyears[,Domain := factor(domain, levels = c("Full","Outer","Middle","Inner"))]
#######################
##ANALYSIS
#######################


#num_years = 3 (for 1000 combos of 3 years)
#resample count = # of times to run sort(sample())

breakup_randomcombos <-function(data, num_years, sampling_period, resample_count = 1000, linear_model = "lm", beta_term = "bray_curtis_dissimilarity_total_mean"){ #window is the size of the window we want to use, linear_model could be however we want to regress year~dissimilarity
  remaining<-data #create dummy data set to operate on
 # output<-data.table(year=integer(0), #create empty data frame to put our output variables in
 #                    length=integer(0), 
 #                    years=integer(0),
 #                    slope=numeric(0), 
 #                    slope_se=numeric(0), 
 #                    p_value=numeric(0),
 #                    intercept=numeric(0), 
 #                    intercept_se=numeric(0), 
 #                    intercept_p_value=numeric(0),
 #                    r_square=numeric(0),
 #                    adj_r_square=numeric(0),
 #                    linear_model=character(0),
 #                    year_list=character(0))
  
  #vector of 1000 sets of 2 years to include

  subsample_years_model_output = function(data, beta_term = beta_term) {
      #randomly select X sequence of Y years
      max_start_year = max(data$year)-sampling_period+1 #pick max start year for given sampling period length and data set
      
      sampling_period_start = sample(seq(min(data$year),max_start_year,by = 1), size = 1) #randomly pick start year
      
      sampling_period_all_years = seq(sampling_period_start, sampling_period_start+sampling_period-1, by = 1)
      
      sampling_period_intermediate_years = sampling_period_all_years[2:(sampling_period-1)]
      
      years_subset_vector = sort(c(sampling_period_start,
                                   sample(x = sampling_period_intermediate_years, size = (num_years-2), replace = F),
                                   max(sampling_period_all_years))) #randomly select X years
      
      chunk=remaining[year %in% years_subset_vector,] #pull out a chunk with only *num_years* randomly sampled years
      
        out = eval(call("linefit",chunk, linear_model = linear_model,beta_term = beta_term))  #fit a linear model (either lm or theil) and get relevant statistics on chunk
                                 
        #add column with list of years included
       # out[,year_list := paste(years_subset_vector, collapse = ",")]
        
      
         # output<-rbind(output, out, use.names = F) #append the stats to the output data frame
          return(out)
  }
  
  reps <- replicate(resample_count, subsample_years_model_output(EBS.distances_dissimilarities_allyears[domain == "Full"], beta_term = beta_term), simplify = FALSE)

  mod_outputs_persubsample_count <- data.table(do.call(rbind, reps))
    
  names(mod_outputs_persubsample_count)<-c("start_year", "N_data", "study_duration", "slope", "slope_se","lower_CI","upper_CI", "p_value",
                   "intercept", "intercept_se", "intercept_p_value","r_square",
                   "adj_r_square", "linear_model")
  return(mod_outputs_persubsample_count)#output the data table
}

#empty data table to populate
subsampling_output_full <- data.table(start_year = as.numeric(), N_data = as.numeric(), study_duration = as.numeric(), slope = as.numeric(), slope_se = as.numeric(),
                                      p_value = as.numeric(),lower_CI = as.numeric(),upper_CI = as.numeric(),
                                         intercept = as.numeric(), intercept_se = as.numeric(), intercept_p_value = as.numeric(),r_square = as.numeric(),
                                         adj_r_square = as.numeric(), linear_model = as.numeric())

#2 years across 5, 10, 15, 20, 25, 30, 35 years
#3 years across 5, 10, 15, 20, 25, 30, 35 years
#5 years across 5, 10, 15, 20, 25, 30, 35 years
#10 years across 10, 15, 20, 25, 30, 35 years
#15 years across 15, 20, 25, 30, 35 years
#20 years across 20, 25, 30, 35 years
#25 years across 25, 30, 35 years
#30 years across 30, 35 years

for (i in c(2,3,5,10,15,20,25,30)){ #number of years of data
  for (j in c(5, 10, 15, 20, 25, 30, 35)) {#length of sampling period
    if(i > j) {
      next
      }else{
  final <- breakup_randomcombos(EBS.distances_dissimilarities_allyears[domain == "Full"],
                                num_years = i, sampling_period = j, resample_count = 5000, linear_model = "lm")
  
  subsampling_output_full <- rbind(subsampling_output_full, final)
  
      }
    print(paste0("Number of years of observations = ",i,", Sampling period = ",j))
  }
}

#save output as csv
fwrite(subsampling_output_full, file.path("Output","subsampling_output_full.csv"))

#subsampling_output_full<-fread(file.path("Output","subsampling_output_full.csv"))
#plot

#slope of full dataset
Full_slope <- coef(lm(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))[[2]]
Full_confint_lower <- confint(lm(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,1] #lower bound of CI for parameter
Full_confint_upper <- confint(lm(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,2] #lower bound of CI for parameter
Full_rsquared <- summary(lm(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$r.squared
Full_pvalue <- summary(lm(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$coefficients[2,4]

subsampling_output_full[,study_duration_label := factor(study_duration,
                  labels = c("Study duration: 5 years",  "Study duration: 10 years", "Study duration: 15 years", "Study duration: 20 years", "Study duration: 25 years", "Study duration: 30 years", "Study duration: 35 years"))]


slope_distribution_obs_duration_ggridge <- ggplot(subsampling_output_full, aes(x = slope, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "Slope") +
  geom_vline(xintercept = Full_slope, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(slope_distribution_obs_duration_ggridge, path = file.path("Figures"),
       filename = "slope_distribution_obs_duration_ggridge.jpg", height = 7, width = 10, unit = "in")

rsquared_distribution_obs_duration_ggridge <- ggplot(subsampling_output_full[N_data >2,], aes(x = r_square, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline",bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "R-squared") +
  geom_vline(xintercept = Full_rsquared, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(rsquared_distribution_obs_duration_ggridge, path = file.path("Figures"),
       filename = "rsquared_distribution_obs_duration_ggridge.jpg", height = 7, width = 10, unit = "in")

pvalue_distribution_obs_duration_ggridge <- ggplot(subsampling_output_full[N_data >2,], aes(x = p_value, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "P-value") +
  geom_vline(xintercept = Full_pvalue, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0.05, color = "red") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(pvalue_distribution_obs_duration_ggridge, path = file.path("Figures"),
       filename = "pvalue_distribution_obs_duration_ggridge.jpg", height = 7, width = 10, unit = "in")

#percent significant by N_year and study_duration
subsampling_output_full[,percent_match :=
                          sum(slope < Full_confint_upper & slope > Full_confint_lower, na.rm = T)/.N,
                        .(N_data,study_duration)][,mean_slope := mean(slope,na.rm. = T),.(N_data, study_duration)][,
                        slope_SD := sd(slope),.(N_data, study_duration)]

subsampling_output_full.u <- unique(subsampling_output_full[,.(N_data,study_duration, percent_match, mean_slope, slope_SD)])

###############
##heatmaps
###############

#heat map with # years,#obs, and %match value
heatmap_year_obs_pmatch <- ggplot(subsampling_output_full.u, aes(y = factor(N_data), x = factor(study_duration), fill= percent_match)) + 
  geom_tile() +
  #geom_text(aes(label= signif(percent_match,2)),  size =3) +
  scale_fill_gradientn(colors = c("#0d2bb5","#0B59B3","#1BA6D7", "#2EE8ED","#86FAF1"), guide = guide_colorbar(frame.colour = "black", ticks.colour = NA,
                     title.position = "top",
                     title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "P(match)") +
  theme_classic() +
  theme(legend.position = c(0.2,0.87),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_pmatch, path = file.path("Figures"),
       filename = "heatmap_year_obs_pmatch.jpg", height = 7, width = 9, unit = "in")




#heat map with # years,#obs, and slope value
heatmap_year_obs_slope <- ggplot(subsampling_output_full.u, aes(y = factor(N_data), x = factor(study_duration), fill= mean_slope)) + 
  geom_tile() +
  #geom_text(aes(label= signif(mean_slope,2)),  size =3) +
  scale_fill_gradientn(
    breaks = c(-0.0018,Full_slope,max(subsampling_output_full.u$mean_slope)),
                       colours = c("#65469c","#3478A2","#72D4AC","white","#F9DFCC","#F6B38F","#DE2D43"),
                       labels =c("-0.002",paste0("-0.00076\nlong\nterm\nslope"),"0.0008"),
                       guide = guide_colorbar(frame.colour = "black", 
                                              #ticks.colour = NA,
                                              title.position = "top",
                                              title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "Mean slope") +
  theme_classic() +
  theme(legend.position = c(0.2,0.85),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_slope, path = file.path("Figures"),
       filename = "heatmap_year_obs_slope.jpg", height = 7, width = 9, unit = "in")


#heat map with # years,#obs, and SD of slope value
heatmap_year_obs_slope_SD <- ggplot(subsampling_output_full.u, aes(y = factor(N_data), x = factor(study_duration), fill= slope_SD)) + 
  geom_tile() +
  #geom_text(aes(label= signif(slope_SD,2)),  size =3) +
  scale_fill_gradientn(colors = c("#FEFFD5","yellow","orange","darkred","#660000"),
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = NA,
                                            title.position = "top",
                                            title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "SD of slopes") +
  theme_classic() +
  theme(legend.position = c(0.2,0.88),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"),
        legend.text = element_text(angle = 90, hjust = 1))

ggsave(heatmap_year_obs_slope_SD, path = file.path("Figures"),
       filename = "heatmap_year_obs_slope_SD.jpg", height = 7, width = 9, unit = "in")

###############
#subsample visuals
###############

f3_d5_slope <- signif(coefficients(lm(bray_curtis_dissimilarity_total_mean~year,
           EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))[2],2)

f3_d5_p_value <- signif(summary(lm(bray_curtis_dissimilarity_total_mean~year,
                               EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))$coefficients[2,4],2)

subsample_viz_f3_d5 <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),], 
              aes(x = year, y = bray_curtis_dissimilarity_total_mean), color = "black",
              method = "lm", se = F, linetype = "longdash") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
 # geom_text(aes(x = 1989, y = 0.43, label = paste0("3 obs over 5 years    ","Slope = ",f3_d5_slope,"\np-value = ",f3_d5_p_value)), size= 3) +
  ggtitle(paste0("3 obs over 5 years      ","slope = ",f3_d5_slope,"      p-value = ",f3_d5_p_value)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f10_d20_slope <- signif(coefficients(lm(bray_curtis_dissimilarity_total_mean~year,
                                      EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))[2],2)

f10_d20_p_value <- signif(summary(lm(bray_curtis_dissimilarity_total_mean~year,
                                   EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))$coefficients[2,4],2)

subsample_viz_f10_d20 <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),], 
              aes(x = year, y = bray_curtis_dissimilarity_total_mean), color = "black",
              method = "lm", se = F, linetype ="dashed") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
#  geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f10_d20_slope,"\np-value = ",f10_d20_p_value)), size= 3) +
  ggtitle(paste0("10 obs over 20 years      ","slope = ",f10_d20_slope,"      p-value = ",f10_d20_p_value)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))



f30_d30_slope <- signif(coefficients(lm(bray_curtis_dissimilarity_total_mean~year,
                                        EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))[2],2)

f30_d30_p_value <- signif(summary(lm(bray_curtis_dissimilarity_total_mean~year,
                                     EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))$coefficients[2,4],2)


subsample_viz_f30_d30 <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),], 
              aes(x = year, y = bray_curtis_dissimilarity_total_mean), color = "black",
              method = "lm", se = F) +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("30 obs over 30 years      ","slope = ",f30_d30_slope,"      p-value = ",f30_d30_p_value)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f2_d30_slope <- signif(coefficients(lm(bray_curtis_dissimilarity_total_mean~year,
                                        EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))[2],2)

f2_d30_p_value <- signif(summary(lm(bray_curtis_dissimilarity_total_mean~year,
                                     EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))$coefficients[2,4],2)

subsample_viz_f2_d30 <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),], 
              aes(x = year, y = bray_curtis_dissimilarity_total_mean), color = "black",
              method = "lm", se = F, linetype = "dashed") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  #geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f2_d30_slope,"\np-value = NA")), size= 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("2 obs over 30 years      ","slope = ",f2_d30_slope,"      p-value = NA")) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

#merge examples
subsample_viz_merge <- plot_grid(subsample_viz_f3_d5 + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                 subsample_viz_f10_d20 + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                 subsample_viz_f30_d30 + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                 subsample_viz_f2_d30 + theme(axis.title.x = element_blank(), text = element_text(size = 20)),
                                 nrow = 4, ncol = 1, labels = c("a.","b.","c.","d."), hjust = 0.5, label_size = 18, label_fontface = "bold", vjust = 1.8)

ggsave(subsample_viz_merge,  path = file.path("Figures"),
       filename = "subsample_viz_merge.jpg", height = 11, width = 4, units = "in")

############################
#examples and heatmap merge
############################

#heatmaps with annotations
heatmap_year_obs_slope_annotate <- heatmap_year_obs_slope + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_pmatch_annotate <- heatmap_year_obs_pmatch + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_slope_SD_annotate <- heatmap_year_obs_slope_SD + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "gray70", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "gray41", size = 8)

#examples and heat map
subsample_heatmap_viz_merge <- 
  plot_grid(subsample_viz_merge, 
            heatmap_year_obs_pmatch_annotate + 
              ggtitle("e.") +
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20),
                    axis.title.x = element_blank()),
            heatmap_year_obs_slope_annotate + ggtitle("f.") + 
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20)), 
            heatmap_year_obs_slope_SD_annotate  + 
              ggtitle("g.") +
              theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                    plot.title = element_text(size = 18, face = "bold"), text = element_text(size = 20)),
            nrow = 2, ncol = 2, align = "hv")

ggsave(subsample_heatmap_viz_merge,  path = file.path("Figures"),
       filename = "subsample_heatmap_viz_merge.jpg", height = 17, width = 16, units = "in")

#######################
##SUPPLEMENT
#######################


#############
##THEIL SEN
############
#empty data table to populate
subsampling_output_full_theilsen <- data.table(start_year = as.numeric(), N_data = as.numeric(), study_duration = as.numeric(), slope = as.numeric(), slope_se = as.numeric(),
                                      p_value = as.numeric(),lower_CI = as.numeric(),upper_CI = as.numeric(),
                                      intercept = as.numeric(), intercept_se = as.numeric(), intercept_p_value = as.numeric(),r_square = as.numeric(),
                                      adj_r_square = as.numeric(), linear_model = as.numeric())

#2 years across 5, 10, 15, 20, 25, 30, 35 years
#3 years across 5, 10, 15, 20, 25, 30, 35 years
#5 years across 5, 10, 15, 20, 25, 30, 35 years
#10 years across 10, 15, 20, 25, 30, 35 years
#15 years across 15, 20, 25, 30, 35 years
#20 years across 20, 25, 30, 35 years
#25 years across 25, 30, 35 years
#30 years across 30, 35 years

for (i in c(2,3,5,10,15,20,25,30)){ #number of years of data
  for (j in c(5, 10, 15, 20, 25, 30, 35)) {#length of sampling period
    if(i > j) {
      next
    }else{
      final <- breakup_randomcombos(EBS.distances_dissimilarities_allyears[domain == "Full"],
                                    num_years = i, sampling_period = j, resample_count = 5000,
                                    linear_model = "theil_sen_regression")
      
      subsampling_output_full_theilsen <- rbind(subsampling_output_full_theilsen, final)
      
    }
    print(paste0("Number of years of observations = ",i,", Sampling period = ",j))
  }
}



#save output as csv
fwrite(subsampling_output_full_theilsen, file.path("Output","Supplement","Theil_Sen","subsampling_output_full_theilsen.csv"))
subsampling_output_full_theilsen <- fread(file.path("Output","Supplement","Theil_Sen","subsampling_output_full_theilsen.csv"))

#plot

#slope of full dataset
Full_slope_theilsen <- coef(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))[[2]]
Full_confint_lower_theilsen <- confint(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,1] #lower bound of CI for parameter
Full_confint_upper_theilsen <- confint(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,2] #lower bound of CI for parameter
Full_rsquared_theilsen <- summary(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$r.squared
Full_pvalue_theilsen <- summary(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$coefficients[2,4]

subsampling_output_full_theilsen[,study_duration_label := factor(study_duration,
                                                        labels = c("Study duration: 5 years",  "Study duration: 10 years", "Study duration: 15 years", "Study duration: 20 years", "Study duration: 25 years", "Study duration: 30 years", "Study duration: 35 years"))]


slope_distributisubsampling_output_full_theilsen <- ggplot(subsampling_output_full_theilsen, aes(x = slope, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "Slope") +
  geom_vline(xintercept = Full_slope_theilsen, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(slope_distributisubsampling_output_full_theilsen, path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "slope_distributisubsampling_output_full_theilsen.jpg", height = 7, width = 10, unit = "in")

rsquared_distribution_obs_duration_ggridge_theilsen <- ggplot(subsampling_output_full_theilsen[N_data >2,], aes(x = r_square, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline",bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "R-squared") +
  geom_vline(xintercept = Full_rsquared_theilsen, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(rsquared_distribution_obs_duration_ggridge_theilsen, path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "rsquared_distribution_obs_duration_ggridge_theilsen.jpg", height = 7, width = 10, unit = "in")

pvalue_distribution_obs_duration_ggridge_theilsen <- ggplot(subsampling_output_full_theilsen[N_data >2,], aes(x = p_value, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "P-value") +
  geom_vline(xintercept = Full_pvalue_theilsen, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0.05, color = "red") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(pvalue_distribution_obs_duration_ggridge_theilsen, path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "pvalue_distribution_obs_duration_ggridge_theilsen.jpg", height = 7, width = 10, unit = "in")

#percent significant by N_year and study_duration
subsampling_output_full_theilsen[,percent_match :=
                          sum(slope < Full_confint_upper_theilsen & slope > Full_confint_lower_theilsen, na.rm = T)/.N,
                        .(N_data,study_duration)][,mean_slope := mean(slope,na.rm. = T),.(N_data, study_duration)][,
                                                                                                                   slope_SD := sd(slope),.(N_data, study_duration)]

subsampling_output_full_theilsen.u <- unique(subsampling_output_full_theilsen[,.(N_data,study_duration, percent_match, mean_slope, slope_SD)])

###############
##heatmaps
###############

#heat map with # years,#obs, and %match value
heatmap_year_obs_pmatch_theilsen <- ggplot(subsampling_output_full_theilsen.u, aes(y = factor(N_data),
                                                                                   x = factor(study_duration), fill= percent_match)) + 
  geom_tile() +
  #geom_text(aes(label= signif(percent_match,2)),  size =3) +
  scale_fill_gradientn(colors = c("#0d2bb5","#0B59B3","#1BA6D7", "#2EE8ED","#86FAF1"), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = NA,
                                                                                                              title.position = "top",
                                                                                                              title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "P(match)") +
  theme_classic() +
  theme(legend.position = c(0.2,0.87),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_pmatch_theilsen, path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "heatmap_year_obs_pmatch_theilsen.jpg", height = 7, width = 9, unit = "in")



#EDIT BREAKS AND LABELS
#heat map with # years,#obs, and slope value
heatmap_year_obs_slope_theilsen <- ggplot(subsampling_output_full_theilsen.u, aes(y = factor(N_data), x = factor(study_duration), fill= mean_slope)) + 
  geom_tile() +
  #geom_text(aes(label= signif(mean_slope,2)),  size =3) +
  scale_fill_gradientn(
    breaks = c(-0.0017,Full_slope_theilsen,max(subsampling_output_full_theilsen.u$mean_slope)),
    colours = c("#65469c","#3478A2","#72D4AC","white","#F6B38F","#DE2D43"),
    labels =c("-0.002",paste0("-0.00034\nlong\nterm\nslope"),".0011"),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",
                           title.position = "top",
                           title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "Mean slope") +
  theme_classic() +
  theme(legend.position = c(0.2,0.85),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_slope_theilsen, path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "heatmap_year_obs_slope_theilsen.jpg", height = 7, width = 9, unit = "in")


#heat map with # years,#obs, and SD of slope value
heatmap_year_obs_slope_SD_theilsen <- ggplot(subsampling_output_full_theilsen.u, aes(y = factor(N_data), x = factor(study_duration), fill= slope_SD)) + 
  geom_tile() +
  #geom_text(aes(label= signif(slope_SD,2)),  size =3) +
  scale_fill_gradientn(colors = c("#FEFFD5","yellow","orange","darkred","#660000"),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = NA,
                                              title.position = "top",
                                              title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "SD of slopes") +
  theme_classic() +
  theme(legend.position = c(0.2,0.88),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"),
        legend.text=element_text(angle = 90, hjust = 1.2))

ggsave(heatmap_year_obs_slope_SD_theilsen, path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "heatmap_year_obs_slope_SD_theilsen.jpg", height = 7, width = 9, unit = "in")

###############
#subsample visuals
###############

######NEXT: fix these to show Theil Sen trends instead of basic lm

f3_d5_slope_theilsen <- signif(coefficients(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                      EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))[2],2)

f3_d5_p_value_theilsen <- signif(summary(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                   EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))$coefficients[2,4],2)

subset_f3_d5 <- EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]
subset_f3_d5[,theilsen_estimate := predict(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))]

subsample_viz_f3_d5_theilsen <-  ggplot() +
  geom_line(data = subset_f3_d5, 
              aes(x = year, y = theilsen_estimate), color = "black",linetype = "longdash") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  # geom_text(aes(x = 1989, y = 0.43, label = paste0("3 obs over 5 years    ","Slope = ",f3_d5_slope,"\np-value = ",f3_d5_p_value)), size= 3) +
  ggtitle(paste0("3 obs over 5 years      ","slope = ",f3_d5_slope_theilsen,"      p-value = ",f3_d5_p_value_theilsen)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f10_d20_slope_theilsen <- signif(coefficients(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                        EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))[2],2)

f10_d20_p_value_theilsen <- signif(summary(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                     EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))$coefficients[2,4],2)

subset_f10_d20 <- EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]
subset_f10_d20[,theilsen_estimate := predict(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))]

subsample_viz_f10_d20_theilsen <-  ggplot() +
  geom_line(data = subset_f10_d20, 
              aes(x = year, y = theilsen_estimate), color = "black",linetype ="dashed") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  #  geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f10_d20_slope,"\np-value = ",f10_d20_p_value)), size= 3) +
  ggtitle(paste0("10 obs over 20 years      ","slope = ",f10_d20_slope_theilsen,"      p-value = ",f10_d20_p_value_theilsen)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))



f30_d30_slope_theilsen <- signif(coefficients(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                        EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))[2],2)

f30_d30_p_value_theilsen <- signif(summary(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                     EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))$coefficients[2,4],2)

subset_f30_d30 <- EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]
subset_f30_d30[,theilsen_estimate := predict(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))]


subsample_viz_f30_d30_theilsen <-  ggplot() +
  geom_line(data = subset_f30_d30, 
              aes(x = year, y = theilsen_estimate), color = "black") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("30 obs over 30 years      ","slope = ",f30_d30_slope_theilsen,"      p-value = ",f30_d30_p_value_theilsen)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f2_d30_slope_theilsen <- signif(coefficients(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                       EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))[2],2)

f2_d30_p_value_theilsen <- signif(summary(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year,
                                    EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))$coefficients[2,4],2)


subset_f2_d30 <- EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]
subset_f2_d30[,theilsen_estimate := predict(theil_sen_regression(bray_curtis_dissimilarity_total_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))]


subsample_viz_f2_d30_theilsen <-  ggplot() +
  geom_line(data = subset_f2_d30, 
              aes(x = year, y = theilsen_estimate), color = "black",linetype = "dashed") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),],
             aes(x = year, y = bray_curtis_dissimilarity_total_mean), size = 3) +
  #geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f2_d30_slope,"\np-value = NA")), size= 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("2 obs over 30 years      ","slope = ",f2_d30_slope_theilsen,"      p-value = NA")) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

#merge examples
subsample_viz_merge_theilsen <- plot_grid(subsample_viz_f3_d5_theilsen + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                 subsample_viz_f10_d20_theilsen + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                 subsample_viz_f30_d30_theilsen + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                 subsample_viz_f2_d30_theilsen + theme(axis.title.x = element_blank(), text = element_text(size = 20)),
                                 nrow = 4, ncol = 1, labels = c("a.","b.","c.","d."), hjust = 0.5, label_size = 18, label_fontface = "bold", vjust = 1.8)

ggsave(subsample_viz_merge_theilsen,  path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "subsample_viz_merge_theilsen.jpg", height = 11, width = 4, units = "in")

############################
#examples and heatmap merge
############################

#heatmaps with annotations
heatmap_year_obs_slope_annotate_theilsen <- heatmap_year_obs_slope_theilsen + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_pmatch_annotate_theilsen <- heatmap_year_obs_pmatch_theilsen + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_slope_SD_annotate_theilsen <- heatmap_year_obs_slope_SD_theilsen + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "gray70", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "gray41", size = 8)

#examples and heat map
subsample_heatmap_viz_merge_theilsen <- 
  plot_grid(subsample_viz_merge_theilsen, 
            heatmap_year_obs_pmatch_annotate_theilsen + 
              ggtitle("e.") +
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20),
                    axis.title.x = element_blank()),
            heatmap_year_obs_slope_annotate_theilsen + ggtitle("f.") + 
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20)), 
            heatmap_year_obs_slope_SD_annotate_theilsen  + 
              ggtitle("g.") +
              theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                    plot.title = element_text(size = 18, face = "bold"), text = element_text(size = 20)),
            nrow = 2, ncol = 2, align = "hv")

ggsave(subsample_heatmap_viz_merge_theilsen,  path = file.path("Figures","Supplement","Theil_Sen"),
       filename = "subsample_heatmap_viz_merge_theilsen.jpg", height = 17, width = 16, units = "in")


############
#Biomass Gradient
############

#empty data table to populate
subsampling_output_full_biomass_gradient <- data.table(start_year = as.numeric(), N_data = as.numeric(), study_duration = as.numeric(), slope = as.numeric(), slope_se = as.numeric(),
                                               p_value = as.numeric(),lower_CI = as.numeric(),upper_CI = as.numeric(),
                                               intercept = as.numeric(), intercept_se = as.numeric(), intercept_p_value = as.numeric(),r_square = as.numeric(),
                                               adj_r_square = as.numeric(), linear_model = as.numeric())

#2 years across 5, 10, 15, 20, 25, 30, 35 years
#3 years across 5, 10, 15, 20, 25, 30, 35 years
#5 years across 5, 10, 15, 20, 25, 30, 35 years
#10 years across 10, 15, 20, 25, 30, 35 years
#15 years across 15, 20, 25, 30, 35 years
#20 years across 20, 25, 30, 35 years
#25 years across 25, 30, 35 years
#30 years across 30, 35 years

for (i in c(2,3,5,10,15,20,25,30)){ #number of years of data
  for (j in c(5, 10, 15, 20, 25, 30, 35)) {#length of sampling period
    if(i > j) {
      next
    }else{
      final <- breakup_randomcombos(EBS.distances_dissimilarities_allyears[domain == "Full"],
                                    num_years = i, sampling_period = j, resample_count = 5000,
                                    beta_term =  "bray_curtis_dissimilarity_gradient_mean")
      
      subsampling_output_full_biomass_gradient <- rbind(subsampling_output_full_biomass_gradient, final)
      
    }
    print(paste0("Number of years of observations = ",i,", Sampling period = ",j))
  }
}



#save output as csv
fwrite(subsampling_output_full_biomass_gradient, file.path("Output","Supplement","Biomass_gradient","subsampling_output_full_biomass_gradient.csv"))

#plot

#slope of full dataset
Full_slope_biomass_gradient <- coef(lm(bray_curtis_dissimilarity_gradient_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))[[2]]
Full_confint_lower_biomass_gradient <- confint(lm(bray_curtis_dissimilarity_gradient_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,1] #lower bound of CI for parameter
Full_confint_upper_biomass_gradient <- confint(lm(bray_curtis_dissimilarity_gradient_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,2] #lower bound of CI for parameter
Full_rsquared_biomass_gradient <- summary(lm(bray_curtis_dissimilarity_gradient_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$r.squared
Full_pvalue_biomass_gradient <- summary(lm(bray_curtis_dissimilarity_gradient_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$coefficients[2,4]

subsampling_output_full_biomass_gradient[,study_duration_label := factor(study_duration,
                                                                 labels = c("Study duration: 5 years",  "Study duration: 10 years", "Study duration: 15 years", "Study duration: 20 years", "Study duration: 25 years", "Study duration: 30 years", "Study duration: 35 years"))]


slope_distributisubsampling_output_full_biomass_gradient <- ggplot(subsampling_output_full_biomass_gradient, aes(x = slope, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "Slope") +
  geom_vline(xintercept = Full_slope_biomass_gradient, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(slope_distributisubsampling_output_full_biomass_gradient, path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "slope_distributisubsampling_output_full_biomass_gradient.jpg", height = 7, width = 10, unit = "in")

rsquared_distribution_obs_duration_ggridge_biomass_gradient <- ggplot(subsampling_output_full_biomass_gradient[N_data >2,], aes(x = r_square, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline",bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "R-squared") +
  geom_vline(xintercept = Full_rsquared_biomass_gradient, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(rsquared_distribution_obs_duration_ggridge_biomass_gradient, path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "rsquared_distribution_obs_duration_ggridge_biomass_gradient.jpg", height = 7, width = 10, unit = "in")

pvalue_distribution_obs_duration_ggridge_biomass_gradient <- ggplot(subsampling_output_full_biomass_gradient[N_data >2,], aes(x = p_value, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "P-value") +
  geom_vline(xintercept = Full_pvalue_biomass_gradient, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0.05, color = "red") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(pvalue_distribution_obs_duration_ggridge_biomass_gradient, path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "pvalue_distribution_obs_duration_ggridge_biomass_gradient.jpg", height = 7, width = 10, unit = "in")

#percent significant by N_year and study_duration
subsampling_output_full_biomass_gradient[,percent_match :=
                                   sum(slope < Full_confint_upper_biomass_gradient & slope > Full_confint_lower_biomass_gradient, na.rm = T)/.N,
                                 .(N_data,study_duration)][,mean_slope := mean(slope,na.rm. = T),.(N_data, study_duration)][,
                                                                                                                            slope_SD := sd(slope),.(N_data, study_duration)]

subsampling_output_full_biomass_gradient.u <- unique(subsampling_output_full_biomass_gradient[,.(N_data,study_duration, percent_match, mean_slope, slope_SD)])

###############
##heatmaps
###############

#heat map with # years,#obs, and %match value
heatmap_year_obs_pmatch_biomass_gradient <- ggplot(subsampling_output_full_biomass_gradient.u, aes(y = factor(N_data),
                                                                                   x = factor(study_duration), fill= percent_match)) + 
  geom_tile() +
  #geom_text(aes(label= signif(percent_match,2)),  size =3) +
  scale_fill_gradientn(colors = c("#0d2bb5","#0B59B3","#1BA6D7", "#2EE8ED","#86FAF1"), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = NA,
                                              title.position = "top",
                                              title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "P(match)") +
  theme_classic() +
  theme(legend.position = c(0.2,0.87),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_pmatch_biomass_gradient, path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "heatmap_year_obs_pmatch_biomass_gradient.jpg", height = 7, width = 9, unit = "in")



#heat map with # years,#obs, and slope value
heatmap_year_obs_slope_biomass_gradient <- ggplot(subsampling_output_full_biomass_gradient.u, aes(y = factor(N_data), x = factor(study_duration), fill= mean_slope)) + 
  geom_tile() +
  #geom_text(aes(label= signif(mean_slope,2)),  size =3) +
  scale_fill_gradientn(
    breaks = c(-0.0003,Full_slope_biomass_gradient,max(subsampling_output_full_biomass_gradient.u$mean_slope)),
    colours = c("#65469c","#3478A2","#72D4AC","white","#F6B38F","#DE2D43","#7F1F5A"),
    labels =c("-0.0003",paste0("-8.6e5\nlong\nterm\nslope"),"0.00015"),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",
                           title.position = "top",
                           title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "Mean slope") +
  theme_classic() +
  theme(legend.position = c(0.2,0.85),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_slope_biomass_gradient, path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "heatmap_year_obs_slope_biomass_gradient.jpg", height = 7, width = 9, unit = "in")


#heat map with # years,#obs, and SD of slope value
heatmap_year_obs_slope_SD_biomass_gradient <- ggplot(subsampling_output_full_biomass_gradient.u, aes(y = factor(N_data), x = factor(study_duration), fill= slope_SD)) + 
  geom_tile() +
  #geom_text(aes(label= signif(slope_SD,2)),  size =3) +
  scale_fill_gradientn(colors = c("#FEFFD5","yellow","orange","darkred","#660000"),
                       breaks = seq(0,0.006,0.001),
                       labels = c(0,"",0.002,"",0.004,"",0.006),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = NA,
                                              title.position = "top",
                                              title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "SD of slopes") +
  theme_classic() +
  theme(legend.position = c(0.2,0.88),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_slope_SD_biomass_gradient, path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "heatmap_year_obs_slope_SD_biomass_gradient.jpg", height = 7, width = 9, unit = "in")

###############
#subsample visuals
###############

#CHECK SIGNIFICANCE

f3_d5_slope_biomass_gradient <- signif(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                                 EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))[2],2)

f3_d5_p_value_biomass_gradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                              EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))$coefficients[2,4],2)

subsample_viz_f3_d5_biomass_gradient <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),], 
            aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), color = "black",linetype = "longdash", method = "lm", se = F) +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),],
             aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  # geom_text(aes(x = 1989, y = 0.43, label = paste0("3 obs over 5 years    ","Slope = ",f3_d5_slope,"\np-value = ",f3_d5_p_value)), size= 3) +
  ggtitle(paste0("3 obs over 5 years      ","slope = ",f3_d5_slope_biomass_gradient,"      p-value = ",f3_d5_p_value_biomass_gradient)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f10_d20_slope_biomass_gradient <- signif(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                                   EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))[2],2)

f10_d20_p_value_biomass_gradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                                EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))$coefficients[2,4],2)


subsample_viz_f10_d20_biomass_gradient <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),], 
            aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), color = "black",linetype ="dashed", se = F, method = "lm") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),],
             aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  #  geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f10_d20_slope,"\np-value = ",f10_d20_p_value)), size= 3) +
  ggtitle(paste0("10 obs over 20 years      ","slope = ",f10_d20_slope_biomass_gradient,"      p-value = ",f10_d20_p_value_biomass_gradient)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))



f30_d30_slope_biomass_gradient <- signif(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                                   EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))[2],2)

f30_d30_p_value_biomass_gradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                                EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))$coefficients[2,4],2)


subsample_viz_f30_d30_biomass_gradient <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),], 
            aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), color = "black", method = "lm", se = F, linetype = "longdash") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),],
             aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("30 obs over 30 years      ","slope = ",f30_d30_slope_biomass_gradient,"      p-value = ",f30_d30_p_value_biomass_gradient)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f2_d30_slope_biomass_gradient <- signif(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                                  EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))[2],2)

f2_d30_p_value_biomass_gradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                                               EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))$coefficients[2,4],2)


subsample_viz_f2_d30_biomass_gradient <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),], 
            aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), color = "black",linetype = "dashed", method = "lm", se = F) +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),],
             aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), size = 3) +
  #geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f2_d30_slope,"\np-value = NA")), size= 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("2 obs over 30 years      ","slope = ",f2_d30_slope_biomass_gradient,"      p-value = NA")) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

#merge examples
subsample_viz_merge_biomass_gradient <- plot_grid(subsample_viz_f3_d5_biomass_gradient + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                          subsample_viz_f10_d20_biomass_gradient + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                          subsample_viz_f30_d30_biomass_gradient + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                          subsample_viz_f2_d30_biomass_gradient + theme(axis.title.x = element_blank(), text = element_text(size = 20)),
                                          nrow = 4, ncol = 1, labels = c("a.","b.","c.","d."), hjust = 0.5, label_size = 18, label_fontface = "bold", vjust = 1.8)

ggsave(subsample_viz_merge_biomass_gradient,  path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "subsample_viz_merge_biomass_gradient.jpg", height = 11, width = 4, units = "in")

############################
#examples and heatmap merge
############################

#heatmaps with annotations
heatmap_year_obs_slope_annotate_biomass_gradient <- heatmap_year_obs_slope_biomass_gradient + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_pmatch_annotate_biomass_gradient <- heatmap_year_obs_pmatch_biomass_gradient + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_slope_SD_annotate_biomass_gradient <- heatmap_year_obs_slope_SD_biomass_gradient + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "gray70", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "gray41", size = 8)

#examples and heat map
subsample_heatmap_viz_merge_biomass_gradient <- 
  plot_grid(subsample_viz_merge_biomass_gradient, 
            heatmap_year_obs_pmatch_annotate_biomass_gradient + 
              ggtitle("e.") +
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20),
                    axis.title.x = element_blank()),
            heatmap_year_obs_slope_annotate_biomass_gradient + ggtitle("f.") + 
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20)), 
            heatmap_year_obs_slope_SD_annotate_biomass_gradient  + 
              ggtitle("g.") +
              theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                    plot.title = element_text(size = 18, face = "bold"), text = element_text(size = 20)),
            nrow = 2, ncol = 2, align = "hv")

ggsave(subsample_heatmap_viz_merge_biomass_gradient,  path = file.path("Figures","Supplement","Biomass_gradient"),
       filename = "subsample_heatmap_viz_merge_biomass_gradient.jpg", height = 17, width = 16, units = "in")



############
#Jaccard
############

#empty data table to populate
subsampling_output_full_jaccard <- data.table(start_year = as.numeric(), N_data = as.numeric(), study_duration = as.numeric(), slope = as.numeric(), slope_se = as.numeric(),
                                                       p_value = as.numeric(),lower_CI = as.numeric(),upper_CI = as.numeric(),
                                                       intercept = as.numeric(), intercept_se = as.numeric(), intercept_p_value = as.numeric(),r_square = as.numeric(),
                                                       adj_r_square = as.numeric(), linear_model = as.numeric())

#2 years across 5, 10, 15, 20, 25, 30, 35 years
#3 years across 5, 10, 15, 20, 25, 30, 35 years
#5 years across 5, 10, 15, 20, 25, 30, 35 years
#10 years across 10, 15, 20, 25, 30, 35 years
#15 years across 15, 20, 25, 30, 35 years
#20 years across 20, 25, 30, 35 years
#25 years across 25, 30, 35 years
#30 years across 30, 35 years

for (i in c(2,3,5,10,15,20,25,30)){ #number of years of data
  for (j in c(5, 10, 15, 20, 25, 30, 35)) {#length of sampling period
    if(i > j) {
      next
    }else{
      final <- breakup_randomcombos(EBS.distances_dissimilarities_allyears[domain == "Full"],
                                    num_years = i, sampling_period = j, resample_count = 5000,
                                    beta_term =  "jaccard_dissimilarity_total_mean")
      
      subsampling_output_full_jaccard <- rbind(subsampling_output_full_jaccard, final)
      
    }
    print(paste0("Number of years of observations = ",i,", Sampling period = ",j))
  }
}



#save output as csv
fwrite(subsampling_output_full_jaccard, file.path("Output","Supplement","Jaccard","subsampling_output_full_jaccard.csv"))
#subsampling_output_full_jaccard <- fread(file.path("Output","Supplement","Jaccard","subsampling_output_full_jaccard.csv"))


#plot

#slope of full dataset
Full_slope_jaccard <- coef(lm(jaccard_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))[[2]]
Full_confint_lower_jaccard <- confint(lm(jaccard_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,1] #lower bound of CI for parameter
Full_confint_upper_jaccard <- confint(lm(jaccard_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]), level = 0.95)[2,2] #lower bound of CI for parameter
Full_rsquared_jaccard <- summary(lm(jaccard_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$r.squared
Full_pvalue_jaccard <- summary(lm(jaccard_dissimilarity_total_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$coefficients[2,4]

subsampling_output_full_jaccard[,study_duration_label := factor(study_duration,
                                                                         labels = c("Study duration: 5 years",  "Study duration: 10 years", "Study duration: 15 years", "Study duration: 20 years", "Study duration: 25 years", "Study duration: 30 years", "Study duration: 35 years"))]


slope_distributisubsampling_output_full_jaccard <- ggplot(subsampling_output_full_jaccard, aes(x = slope, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "Slope") +
  geom_vline(xintercept = Full_slope_jaccard, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(slope_distributisubsampling_output_full_jaccard, path = file.path("Figures","Supplement","Jaccard"),
       filename = "slope_distributisubsampling_output_full_jaccard.jpg", height = 7, width = 10, unit = "in")

rsquared_distribution_obs_duration_ggridge_jaccard <- ggplot(subsampling_output_full_jaccard[N_data >2,], aes(x = r_square, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline",bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "R-squared") +
  geom_vline(xintercept = Full_rsquared_jaccard, color = "grey", linetype = "dashed") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(rsquared_distribution_obs_duration_ggridge_jaccard, path = file.path("Figures","Supplement","Jaccard"),
       filename = "rsquared_distribution_obs_duration_ggridge_jaccard.jpg", height = 7, width = 10, unit = "in")

pvalue_distribution_obs_duration_ggridge_jaccard <- ggplot(subsampling_output_full_jaccard[N_data >2,], aes(x = p_value, y = as.factor(N_data), group = as.factor(N_data), fill = as.factor(N_data))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "P-value") +
  geom_vline(xintercept = Full_pvalue_jaccard, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0.05, color = "red") +
  facet_wrap(~study_duration_label) +
  theme_classic() +
  theme(legend.position = "null")

ggsave(pvalue_distribution_obs_duration_ggridge_jaccard, path = file.path("Figures","Supplement","Jaccard"),
       filename = "pvalue_distribution_obs_duration_ggridge_jaccard.jpg", height = 7, width = 10, unit = "in")

#percent significant by N_year and study_duration
subsampling_output_full_jaccard[,percent_match :=
                                           sum(slope < Full_confint_upper_jaccard & slope > Full_confint_lower_jaccard, na.rm = T)/.N,
                                         .(N_data,study_duration)][,mean_slope := mean(slope,na.rm. = T),.(N_data, study_duration)][,
                                                                                                                                    slope_SD := sd(slope),.(N_data, study_duration)]

subsampling_output_full_jaccard.u <- unique(subsampling_output_full_jaccard[,.(N_data,study_duration, percent_match, mean_slope, slope_SD)])

###############
##heatmaps
###############

#heat map with # years,#obs, and %match value
heatmap_year_obs_pmatch_jaccard <- ggplot(subsampling_output_full_jaccard.u, aes(y = factor(N_data),
                                                                                                   x = factor(study_duration), fill= percent_match)) + 
  geom_tile() +
  #geom_text(aes(label= signif(percent_match,2)),  size =3) +
  scale_fill_gradientn(colors = c("#0d2bb5","#0B59B3","#1BA6D7", "#2EE8ED","#86FAF1"), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = NA,
                                              title.position = "top",
                                              title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "P(match)") +
  theme_classic() +
  theme(legend.position = c(0.2,0.87),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_pmatch_jaccard, path = file.path("Figures","Supplement","Jaccard"),
       filename = "heatmap_year_obs_pmatch_jaccard.jpg", height = 7, width = 9, unit = "in")



#heat map with # years,#obs, and slope value
heatmap_year_obs_slope_jaccard <- ggplot(subsampling_output_full_jaccard.u, aes(y = factor(N_data), x = factor(study_duration), fill= mean_slope)) + 
  geom_tile() +
  #geom_text(aes(label= signif(mean_slope,2)),  size =3) +
  scale_fill_gradientn(
    breaks = c(-0.001,Full_slope_jaccard,max(subsampling_output_full_jaccard.u$mean_slope)),
      colours = c("#65469c","#3478A2","#72D4AC","#c3fae3","#ebfcf5","white","#F9DFCC","#F6B38F","#DE2D43","#7F1F5A"),
       labels =c("-0.001",paste0("-1.5e-4\nlong\nterm\nslope"),"0.0008"),
     guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",
    title.position = "top",
    title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "Mean slope") +
  theme_classic() +
  theme(legend.position = c(0.2,0.85),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_slope_jaccard, path = file.path("Figures","Supplement","Jaccard"),
       filename = "heatmap_year_obs_slope_jaccard.jpg", height = 7, width = 9, unit = "in")


#heat map with # years,#obs, and SD of slope value
heatmap_year_obs_slope_SD_jaccard <- ggplot(subsampling_output_full_jaccard.u, aes(y = factor(N_data), x = factor(study_duration), fill= slope_SD)) + 
  geom_tile() +
  #geom_text(aes(label= signif(slope_SD,2)),  size =3) +
  scale_fill_gradientn(colors = c("#FEFFD5","yellow","orange","darkred","#660000"),
                       labels = c("0","0.002","0.004","0.006","0.008"),
                       breaks = seq(0.000,0.008,0.002),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = NA,
                                              title.position = "top",
                                              title.hjust = 0.5)) +
  labs(x = "Study duration", y = "Number of years of data", fill = "SD of slopes") +
  theme_classic() +
  theme(legend.position = c(0.2,0.88),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"))

ggsave(heatmap_year_obs_slope_SD_jaccard, path = file.path("Figures","Supplement","Jaccard"),
       filename = "heatmap_year_obs_slope_SD_jaccard.jpg", height = 7, width = 9, unit = "in")

###############
#subsample visuals
###############

#CHECK SIGNIFICANCE

f3_d5_slope_jaccard <- signif(coefficients(lm(jaccard_dissimilarity_total_mean~year,
                                                       EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))[2],2)

f3_d5_p_value_jaccard <- signif(summary(lm(jaccard_dissimilarity_total_mean~year,
                                                    EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),]))$coefficients[2,4],2)


subsample_viz_f3_d5_jaccard <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),], 
            aes(x = year, y = jaccard_dissimilarity_total_mean), color = "black",linetype = "longdash", method = "lm", se = F) +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = jaccard_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1991,1994,1996),],
             aes(x = year, y = jaccard_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  # geom_text(aes(x = 1989, y = 0.43, label = paste0("3 obs over 5 years    ","Slope = ",f3_d5_slope,"\np-value = ",f3_d5_p_value)), size= 3) +
  ggtitle(paste0("3 obs over 5 years      ","slope = ",f3_d5_slope_jaccard,"      p-value = ",f3_d5_p_value_jaccard)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f10_d20_slope_jaccard <- signif(coefficients(lm(jaccard_dissimilarity_total_mean~year,
                                                         EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))[2],2)

f10_d20_p_value_jaccard <- signif(summary(lm(jaccard_dissimilarity_total_mean~year,
                                                      EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),]))$coefficients[2,4],2)

subsample_viz_f10_d20_jaccard <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),], 
            aes(x = year, y = jaccard_dissimilarity_total_mean), color = "black",linetype ="solid", se = F, method = "lm") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = jaccard_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1992, 1993, 1995, 2000, 2001, 2003, 2008, 2009, 2011, 2012),],
             aes(x = year, y = jaccard_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  #  geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f10_d20_slope,"\np-value = ",f10_d20_p_value)), size= 3) +
  ggtitle(paste0("10 obs over 20 years      ","slope = ",f10_d20_slope_jaccard,"      p-value = ",f10_d20_p_value_jaccard)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))



f30_d30_slope_jaccard <- signif(coefficients(lm(jaccard_dissimilarity_total_mean~year,
                                                         EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))[2],2)

f30_d30_p_value_jaccard <- signif(summary(lm(jaccard_dissimilarity_total_mean~year,
                                                      EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),]))$coefficients[2,4],2)

subsample_viz_f30_d30_jaccard <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),], 
            aes(x = year, y = jaccard_dissimilarity_total_mean), color = "black", linetype = "solid", method = "lm", se = F) +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = jaccard_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% seq(1984,2014,1),],
             aes(x = year, y = jaccard_dissimilarity_total_mean), size = 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("30 obs over 30 years      ","slope = ",f30_d30_slope_jaccard,"      p-value = ",f30_d30_p_value_jaccard)) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

f2_d30_slope_jaccard <- signif(coefficients(lm(jaccard_dissimilarity_total_mean~year,
                                                        EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))[2],2)

f2_d30_p_value_jaccard <- signif(summary(lm(jaccard_dissimilarity_total_mean~year,
                                                     EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),]))$coefficients[2,4],2)



subsample_viz_f2_d30_jaccard <-  ggplot() +
  geom_smooth(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),], 
            aes(x = year, y = jaccard_dissimilarity_total_mean), color = "black",linetype = "dashed", se = F, method = "lm") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full",], aes(x = year, y = jaccard_dissimilarity_total_mean), shape = 21, fill = "white", color = "black",  size = 3) +
  labs(x = "Year", y = "β diversity") +
  geom_point(data = EBS.distances_dissimilarities_allyears[Domain == "Full" & year %in% c(1985,2015),],
             aes(x = year, y = jaccard_dissimilarity_total_mean), size = 3) +
  #geom_text(aes(x = 1989, y = 0.43, label = paste0("Slope = ",f2_d30_slope,"\np-value = NA")), size= 3) +
  labs(x = "Year", y = "β diversity") +
  ggtitle(paste0("2 obs over 30 years      ","slope = ",f2_d30_slope_jaccard,"      p-value = NA")) +
  theme_classic() +
  theme(legend.position = "null", plot.title = element_text(size = 14))

#merge examples
subsample_viz_merge_jaccard <- plot_grid(subsample_viz_f3_d5_jaccard + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                                  subsample_viz_f10_d20_jaccard + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                                  subsample_viz_f30_d30_jaccard + theme(axis.title.x = element_blank(),axis.text.x = element_blank(), text = element_text(size = 20)),
                                                  subsample_viz_f2_d30_jaccard + theme(axis.title.x = element_blank(), text = element_text(size = 20)),
                                                  nrow = 4, ncol = 1, labels = c("a.","b.","c.","d."), hjust = 0.5, label_size = 18, label_fontface = "bold", vjust = 1.8)

ggsave(subsample_viz_merge_jaccard,  path = file.path("Figures","Supplement","Jaccard"),
       filename = "subsample_viz_merge_jaccard.jpg", height = 11, width = 4, units = "in")

############################
#examples and heatmap merge
############################

#heatmaps with annotations
heatmap_year_obs_slope_annotate_jaccard <- heatmap_year_obs_slope_jaccard + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_pmatch_annotate_jaccard <- heatmap_year_obs_pmatch_jaccard + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "white", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "white", size = 8)

heatmap_year_obs_slope_SD_annotate_jaccard <- heatmap_year_obs_slope_SD_jaccard + 
  annotate(geom = "text", x = 1,y = 2, label = "a", fontface = "bold", color = "gray70", size = 8) +
  annotate(geom = "text", x = 4,y = 4, label = "b", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 8, label = "c", fontface = "bold", color = "gray41", size = 8) +
  annotate(geom = "text", x = 6,y = 1, label = "d", fontface = "bold", color = "gray41", size = 8)

#examples and heat map
subsample_heatmap_viz_merge_jaccard <- 
  plot_grid(subsample_viz_merge_jaccard, 
            heatmap_year_obs_pmatch_annotate_jaccard + 
              ggtitle("e.") +
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20),
                    axis.title.x = element_blank()),
            heatmap_year_obs_slope_annotate_jaccard + ggtitle("f.") + 
              theme(plot.title = element_text(size = 18, face = "bold"),text = element_text(size = 20)), 
            heatmap_year_obs_slope_SD_annotate_jaccard  + 
              ggtitle("g.") +
              theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                    plot.title = element_text(size = 18, face = "bold"), text = element_text(size = 20)),
            nrow = 2, ncol = 2, align = "hv")

ggsave(subsample_heatmap_viz_merge_jaccard,  path = file.path("Figures","Supplement","Jaccard"),
       filename = "subsample_heatmap_viz_merge_jaccard.jpg", height = 17, width = 16, units = "in")

