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
#######################
##FUNCTION
#######################

#1000 random combinations of 2 years
sort(sample(years, 2, replace = F))
#1000 random combinations of 3 years
sort(sample(years, 3, replace = F))
#1000 random combinations of 5 years
sort(sample(years, 5, replace = F))
#1000 random combinations of 10 years
sort(sample(years, 10, replace = F))
#1000 random combinations of 15 years
sort(sample(years, 15, replace = F))
#1000 random combinations of 20 years
sort(sample(years, 20, replace = F))
#1000 random combinations of 25 years
sort(sample(years, 25, replace = F))
#1000 random combinations of 30 years
sort(sample(years, 30, replace = F))

#num_years = 3 (for 1000 combos of 3 years)
#resample count = # of times to run sort(sample())

breakup_randomcombos <-function(data, num_years, resample_count = 1000, linear_model = "lm"){ #window is the size of the window we want to use, linear_model could be however we want to regress year~dissimilarity
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
  numyears<-length(unique(data$year))

  #vector of 1000 sets of 2 years to include

  subsample_years_model_output = function(data) {
      years_subset_vector = sort(sample(data$year, num_years, replace = F))
      
      chunk=remaining[year %in% years_subset_vector,] #pull out a chunk with only *num_years* randomly sampled years
      
        out = eval(call("linefit",chunk, linear_model = linear_model))  #fit a linear model (either lm or theil) and get relevant statistics on chunk
                                 
        #add column with list of years included
       # out[,year_list := paste(years_subset_vector, collapse = ",")]
        
      
         # output<-rbind(output, out, use.names = F) #append the stats to the output data frame
          return(out)
  }
  
  reps <- replicate(1000, subsample_years_model_output(EBS.distances_dissimilarities_allyears[domain == "Full"]), simplify = FALSE)

  mod_outputs_persubsample_count <- data.table(do.call(rbind, reps))
    
  numyears<-length(unique(remaining$year))

  names(mod_outputs_persubsample_count)<-c("start_year", "N_data", "N_years", "slope", "slope_se", "p_value",
                   "intercept", "intercept_se", "intercept_p_value","r_square",
                   "adj_r_square", "linear_model")
  return(mod_outputs_persubsample_count)#output the data table
}

for (i in c(2, 3, 5, 10, 15, 20, 25, 30)) {
final <- breakup_randomcombos(EBS.distances_dissimilarities_allyears[domain == "Full"],
                              num_years = i, resample_count = 1000, linear_model = "lm")

final.final <- rbind(final.final, final)
}

#plot
library(ggridges)

#slope of full dataset
Full_slope <- coef(lm(bray_curtis_dissimilarity_balanced_mean~year,data = EBS.distances_dissimilarities_allyears[domain == "Full"]))[[2]]

ggplot(final.final, aes(x = slope, y = as.factor(N_years), group = as.factor(N_years), fill = as.factor(N_years))) +
  geom_density_ridges(stat = "binline", bins = 50, scale = 1.3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Number of years of data", x = "Slope") +
  geom_vline(xintercept = Full_slope, color = "grey", linetype = "dashed") +
  theme_classic() +
  theme(legend.position = "null")

ggplot(aes(y = bray_curtis_dissimilarity_balanced_mean, x=year),
       data = EBS.distances_dissimilarities_allyears[domain == "Full"]) +
  geom_point() +
  labs(y = "BC dissimilarity", x = "Year") +
  theme_classic()
    
