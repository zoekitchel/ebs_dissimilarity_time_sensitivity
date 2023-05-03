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

#######################
##ANALYSIS
#######################


#num_years = 3 (for 1000 combos of 3 years)
#resample count = # of times to run sort(sample())

breakup_randomcombos <-function(data, num_years, sampling_period, resample_count = 1000, linear_model = "lm"){ #window is the size of the window we want to use, linear_model could be however we want to regress year~dissimilarity
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

  subsample_years_model_output = function(data) {
      #randomly select X sequence of Y years
      max_start_year = max(data$year)-sampling_period+1 #pick max start year for given sampling period length and data set
      
      sampling_period_start = sample(seq(min(data$year),max_start_year,by = 1), size = 1) #randomly pick start year
      
      sampling_period_all_years = seq(sampling_period_start, sampling_period_start+sampling_period-1, by = 1)
      
      sampling_period_intermediate_years = sampling_period_all_years[2:(sampling_period-1)]
      
      years_subset_vector = sort(c(sampling_period_start,
                                   sample(x = sampling_period_intermediate_years, size = (num_years-2), replace = F),
                                   max(sampling_period_all_years))) #randomly select X years
      
      chunk=remaining[year %in% years_subset_vector,] #pull out a chunk with only *num_years* randomly sampled years
      
        out = eval(call("linefit",chunk, linear_model = "lm"))  #fit a linear model (either lm or theil) and get relevant statistics on chunk
                                 
        #add column with list of years included
       # out[,year_list := paste(years_subset_vector, collapse = ",")]
        
      
         # output<-rbind(output, out, use.names = F) #append the stats to the output data frame
          return(out)
  }
  
  reps <- replicate(1000, subsample_years_model_output(EBS.distances_dissimilarities_allyears[domain == "Full"]), simplify = FALSE)

  mod_outputs_persubsample_count <- data.table(do.call(rbind, reps))
    
  names(mod_outputs_persubsample_count)<-c("start_year", "N_data", "study_duration", "slope", "slope_se", "p_value",
                   "intercept", "intercept_se", "intercept_p_value","r_square",
                   "adj_r_square", "linear_model")
  return(mod_outputs_persubsample_count)#output the data table
}

#empty data table to populate
subsampling_output_full <- data.table(start_year = as.numeric(), N_data = as.numeric(), study_duration = as.numeric(), slope = as.numeric(), slope_se = as.numeric(), p_value = as.numeric(),
                                         intercept = as.numeric(), intercept_se = as.numeric(), intercept_p_value = as.numeric(),r_square = as.numeric(),
                                         adj_r_square = as.numeric(), linear_model = as.numeric())

#2 years across 5, 10, 15, 20, 25, 30, 35 years
#3 years across 5, 10, 15, 20, 25, 30, 35 years
#5 years across 5, 10, 15, 20, 25, 30, 35 years
#10 years across 10, 15, 20, 25, 30, 35 years
#20 years across 20, 25, 30, 35 years

for (i in c(2,3,5,10,20)){
  for (j in c(5, 10, 15, 20, 25, 30, 35)) {
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

#plot

#slope of full dataset
Full_slope <- coef(lm(bray_curtis_dissimilarity_balanced_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))[[2]]
Full_rsquared <- summary(lm(bray_curtis_dissimilarity_balanced_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$r.squared
Full_pvalue <- summary(lm(bray_curtis_dissimilarity_balanced_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))$coefficients[2,4]

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
subsampling_output_full[,percent_significant := sum(p_value <= 0.05)/.N,.(N_data,study_duration)]

subsampling_output_full.u <- unique(subsampling_output_full[,.(N_data,study_duration, percent_significant)])

ggplot(subsampling_output_full[N_data >2,]) +
  geom_point(aes(x = study_duration, y = percent_significant)) +
  facet_wrap(~N_data) +
  labs(x = "Study Duration", y = "Percent of models p < 0.05") +
  theme_classic()
