###########################
#Written by Zoë Kitchel for Eastern Bering Sea β diversity (dissimilarity) Through Time

# How does when you look, and how long, affect the conclusions you reach about your data?
# Are short term studies more likely to yield significant results?
# Are short term studies more likely to find *erroneous* significant trends?
# This script will perform a simple moving window analysis to answer these questions
# for long term data across a variety of domains- essentially, data gets subsetted, 
# we run a simple linear regression on each subset and record summary stats, trends

#assume data is coming in in the form year, response variable
#use this test data set to build stuff

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
#######################
##FUNCTIONS
#######################
#pull in neccessary functions
source('code/FUN_subsamples.R')
source('code/FUN_linefit.R')
#######################
##DATA
#######################
EBS.distances_dissimilarities_allyears <- fread(file.path("Output","EBS.distances_dissimilarities_allyears.csv"))


#simplified for BC balanced variation only

EBS.dissim.simp <- EBS.distances_dissimilarities_allyears[,.(year, bray_curtis_dissimilarity_balanced_mean, domain)]

# next we need a function that runs a simple linear model of x=year, y=response variable

linefit<-function (data, linear_model = "lm"){
  #fit the model
  model<-eval(call(linear_model, bray_curtis_dissimilarity_balanced_mean~year, data=data))
  #create a vector of relevant outputs. We want slope, error, P value
  output<-c(min(data$year), #year the analysis started on
            nrow(data), #number of data points the analysis includes
            length(unique(data$year)), #number of years the analysis includes
            summary(model)$coefficients[2,1], # slope
            summary(model)$coefficients[2,2], # se for slope
            summary(model)$coefficients[2,4], #p value slope
            summary(model)$coefficients[1,1], # intercept
            summary(model)$coefficients[1,2], # se for intercept
            summary(model)$coefficients[1,4], # p value for intercept
            summary(model)$r.squared, #r-squared
            summary(model)$adj.r.squared, #adjusted r-squared
            ifelse(linear_model == "lm",0,1)) 
  return(output)
}

#and try this on test data
#linear model
linefit(EBS.dissim.simp[domain == "Full"])
linefit(EBS.dissim.simp[domain == "Outer"])
linefit(EBS.dissim.simp[domain == "Middle"])
linefit(EBS.dissim.simp[domain == "Inner"])

# functional!

#alternative line fit using Theil–Sen (what Olaf recommended) https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator

#Plot a line between all the points in your data
#Calculate the slope for each line
#The median slope is your regression slope
#Calculating the slope this way happens to be quite robust.
#And when the errors are normally distributed and you have no outliers, the slope is very similar to OLS

#theil-sen
linefit(EBS.dissim.simp[domain == "Full"], linear_model = "theil_sen_regression")
linefit(EBS.dissim.simp[domain == "Outer"], linear_model = "theil_sen_regression")
linefit(EBS.dissim.simp[domain == "Middle"], linear_model = "theil_sen_regression")
linefit(EBS.dissim.simp[domain == "Inner"], linear_model = "theil_sen_regression")


#now we need to think about how to iterate through the dataset. We want a
#function that starts at the first year, counts the number of rows specified
#and then feeds that resultant data frame to the fitting function. 
#then we want to discard the first row of the data set, and repeat until fewer than
#the number of rows specified remains

breakup<-function(data, window, linear_model = "lm"){ #window is the size of the window we want to use, linear_model could be however we want to regress year~dissimilarity
  remaining<-data #create dummy data set to operate on
  output<-data.frame(year=integer(0), #create empty data frame to put our output variables in
                     length=integer(0), 
                     years=integer(0),
                     slope=numeric(0), 
                     slope_se=numeric(0), 
                     p_value=numeric(0),
                     intercept=numeric(0), 
                     intercept_se=numeric(0), 
                     intercept_p_value=numeric(0),
                     r_square=numeric(0),
                     adj_r_square=numeric(0),
                     linear_model=character(0))
  numyears<-length(unique(data$year))
  while (numyears>(window-1)){ #while there's still more years of data than in the window
    chunk<-subset(remaining, year<(min(year)+window)) #pull out a chunk as big as the window from the top of the data
    
    out <- eval(call("linefit",chunk, linear_model = linear_model)) #fit a linear model (either lm or theil) and get relevant statistics on chunk
    #add a conditional so that if there's missing data, it's not included in output
    if (window==length(unique(chunk$year))){
      output<-rbind(output, out) #append the stats to the output data frame
    }else{
      output<-output #leave it out if it has missing data
    }
    
    remaining<-subset(remaining, year>min(year)) #cut out the first year of the remaining data + repeat
    numyears<-length(unique(remaining$year))
  }
  names(output)<-c("start_year", "N_data", "N_years", "slope", "slope_se", "p_value",
                   "intercept", "intercept_se", "intercept_p_value","r_square", "adj_r_square", "linear_model")
  return(output)#output the data frame
}

#and now try this on the test data
#linear model
breakup(EBS.dissim.simp[domain == "Full"], 3)
breakup(EBS.dissim.simp[domain == "Inner"], 3)
breakup(EBS.dissim.simp[domain == "Outer"], 3)
breakup(EBS.dissim.simp[domain == "Middle"], 3)

#theil
breakup(EBS.dissim.simp[domain == "Full"], 3, linear_model = "theil_sen_regression")
breakup(EBS.dissim.simp[domain == "Inner"], 3, linear_model = "theil_sen_regression")
breakup(EBS.dissim.simp[domain == "Outer"], 3, linear_model = "theil_sen_regression")
breakup(EBS.dissim.simp[domain == "Middle"], 3, linear_model = "theil_sen_regression")

# now time to write the function that will iterate through our targetted windows

multiple_breakups<-function(data, linear_model = "lm"){
  count<-length(data$year)
  output<-data.frame(year=integer(0), #create empty data frame to put our output variables in
                     length=integer(0), 
                     years=integer(0),
                     slope=numeric(0), 
                     slope_se=numeric(0), 
                     p_value=numeric(0),
                     intercept=numeric(0), 
                     intercept_se=numeric(0), 
                     intercept_p_value=numeric(0),
                     r_square=numeric(0),
                     adj_r_square=numeric(0),
                     linear_model=character(0))
  for(i in 3:(count)){
    outeach<-breakup(data, i, linear_model = linear_model) #fit at each window length
    output<-rbind(output, outeach)#bind it to the frame
  }
  out<-output
  return(out)
}

#try on data
#linear model
multiple_breakups(EBS.dissim.simp[domain == "Full"])
multiple_breakups(EBS.dissim.simp[domain == "Inner"])
multiple_breakups(EBS.dissim.simp[domain == "Outer"])
multiple_breakups(EBS.dissim.simp[domain == "Middle"])

#theil
multiple_breakups(EBS.dissim.simp[domain == "Full"], linear_model = "theil_sen_regression")
multiple_breakups(EBS.dissim.simp[domain == "Inner"], linear_model = "theil_sen_regression")
multiple_breakups(EBS.dissim.simp[domain == "Outer"], linear_model = "theil_sen_regression")
multiple_breakups(EBS.dissim.simp[domain == "Middle"], linear_model = "theil_sen_regression")
#fan-flipping-tastic! it looks like that works


#let's create a plotting function

pyramid_plot<- function(data, title="", significance=0.05, plot_insig=TRUE, rsq_points=FALSE, linear_model = "lm"){
  out<-multiple_breakups(data, linear_model = linear_model)
  years<-length(unique(out$start_year))
  maxyears<-max(out$N_years)
  count<-nrow(out)
  #compute mean and sd of longest series for vertical lines
  true_slope<-out[count,4] #find the slope of the longest series
  #remember to convert standard error to standard deviation
  true_error<-(out[count,5])*(sqrt(out[count, 2]))#find the error of the longest series
  max_true<-true_slope+true_error #compute max and min values for slopes we are calling true
  min_true<-true_slope-true_error
  out$significance<-ifelse(out$p_value<significance, "Yes", "No")
  if(rsq_points==TRUE){
    point_scale<-10*out$r_square
    yespt<-1
  }else{
    point_scale<-2
    yespt<-16
  }
  if(plot_insig==FALSE){
    out<-out[which(out$p_value<significance),]
  }
  plot<- ggplot(out) +
    theme_classic() +
    geom_hline(yintercept = true_slope, linetype = 2) +
    geom_hline(yintercept = max_true, linetype = 3, color="grey38") + #confidence interval (true slope + true_error for slopes we're calling true)
    geom_hline(yintercept = min_true, linetype = 3, color="grey38") + #confidence interval (true slope + true_error for slopes we're calling true)
    geom_vline(xintercept = 10, linetype = "dashed", color = "grey") + #helpful gridlines
    geom_vline(xintercept = 20, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 30, linetype = "dashed", color = "grey") +
    aes(y = slope, x = N_years,  ymin = (slope-slope_se), 
        ymax = (slope+slope_se), shape=significance, color=significance) +
    geom_linerange(show.legend = F)+ 
    geom_point(size=point_scale)+ ggtitle(title)+
    scale_shape_manual(values=c("No"=4,"Yes"=yespt))+
    scale_color_manual(values=c("No"="red","Yes"="black"))+
    labs(x = "Number of years in window", y = "Slope", shape = "p-value < 0.05", color = "p-value < 0.05")+
    scale_x_continuous(lim=c(3, maxyears))+
    coord_flip()
  return(plot)
}
#linear model
pyramid_plot_full_lm <- pyramid_plot(EBS.dissim.simp[domain == "Full"], title="Full EBS, lm()", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_full_lm, path = file.path("Figures"), filename = "pyramid_plot_full_lm.jpg")
pyramid_plot_inner_lm <- pyramid_plot(EBS.dissim.simp[domain == "Inner"], title="Inner EBS, lm()", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_inner_lm, path = file.path("Figures"), filename = "pyramid_plot_inner_lm.jpg")
pyramid_plot_outer_lm <- pyramid_plot(EBS.dissim.simp[domain == "Outer"], title="Outer EBS, lm()", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_outer_lm, path = file.path("Figures"), filename = "pyramid_plot_outer_lm.jpg")
pyramid_plot_middle_lm <- pyramid_plot(EBS.dissim.simp[domain == "Middle"], title="Middle EBS, lm()", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_middle_lm, path = file.path("Figures"), filename = "pyramid_plot_middle_lm.jpg")

#theil sein
pyramid_plot_full_theil_sen <- pyramid_plot(EBS.dissim.simp[domain == "Full"], linear_model = "theil_sen_regression", title="Full EBS, theil sen", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_full_theil_sen, path = file.path("Figures"), filename = "pyramid_plot_full_theil_sen.jpg")
pyramid_plot_inner_theil_sen <- pyramid_plot(EBS.dissim.simp[domain == "Inner"], linear_model = "theil_sen_regression", title="Inner EBS, theil sen", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_inner_theil_sen, path = file.path("Figures"), filename = "pyramid_plot_inner_theil_sen.jpg")
pyramid_plot_outer_theil_sen <- pyramid_plot(EBS.dissim.simp[domain == "Outer"], linear_model = "theil_sen_regression", title="Outer EBS, theil sen", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_outer_theil_sen, path = file.path("Figures"), filename = "pyramid_plot_outer_theil_sen.jpg")
pyramid_plot_middle_theil_sen <- pyramid_plot(EBS.dissim.simp[domain == "Middle"], linear_model = "theil_sen_regression", title="Middle EBS, theil sen", plot_insig = TRUE, significance=0.05, rsq_points =TRUE)
ggsave(pyramid_plot_middle_theil_sen, path = file.path("Figures"), filename = "pyramid_plot_middle_theil_sen.jpg")

#merge into one figure
pyramid_plot_merge <- plot_grid(
  #lm (top left)
  pyramid_plot_full_lm + lims(y = c(-0.12, 0.06)) +  
    theme(plot.title = element_text(face = "bold", color = "black"), 
          legend.position = c(0.3, 0.65), legend.key.size = unit(0.5, unit = "cm"), legend.title = element_text(size = 16), legend.text = element_text(size = 16),
          axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16)), 
  
  
  pyramid_plot_inner_lm + lims(y = c(-0.12, 0.06)) + 
    theme(plot.title = element_text(face = "bold", color = "#AA4499"), 
          legend.position = "null", 
          axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 16)), 
  
  
  pyramid_plot_middle_lm + lims(y = c(-0.12, 0.06)) + 
    theme(plot.title = element_text(face = "bold", color = "#44AA99"), 
          legend.position = "null", 
          axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 16)), 
  
  
  pyramid_plot_outer_lm + lims(y = c(-0.12, 0.06)) + 
    theme(plot.title = element_text(face = "bold", color = "#999933"), 
          legend.position = "null", 
          axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 16)), 
  
  #theil sen (bottom left)
  pyramid_plot_full_theil_sen + lims(y = c(-0.12, 0.06)) + 
    theme(plot.title = element_text(face = "bold", color = "black"), 
          legend.position = "null", 
          axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)), 
  
  
  pyramid_plot_inner_theil_sen + lims(y = c(-0.12, 0.06)) + 
    theme(plot.title = element_text(face = "bold", color = "#AA4499"), 
          legend.position = "null", 
          axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 16)), 
  
  
  pyramid_plot_middle_theil_sen + lims(y = c(-0.12, 0.06)) + 
    theme(plot.title = element_text(face = "bold", color = "#44AA99"), 
          legend.position = "null", 
          axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 16)),
  
  #bottom right
  pyramid_plot_outer_theil_sen + lims(y = c(-0.12, 0.06)) + 
    theme(plot.title = element_text(face = "bold", color = "#999933"), 
          legend.position = "null", 
          axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 16)), 
  
  
  ncol = 4, nrow = 2)

ggsave(pyramid_plot_merge, path = file.path("Figures"),
       filename = "pyramid_plot_merge.jpg", height = 8, width = 15, unit = "in")


#now that we have visualization, we need a way to pull relevant metrics out of the computation
#so let's say our longest series is our 'truth', and we want to know how many years it takes 
#to reach 'stability'-so let's define stability as >(some percentage of slopes) occuring within 
#the standard deviation of the slope of the longest series, for a given window length, allow user to change # of SEs

#For INSIGNIFICANT, we are just interested in time it takes for an insignificant yet negative slope to be most likely

stability_time <-function(data, min_percent=95, error_multiplyer=1, linear_model = "lm"){#returns a number 
  test<-multiple_breakups(data, linear_model = linear_model)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  #remember to convert standard error to standard deviation
  true_error<-(test[count,5])*(sqrt(test[count, 2]))*error_multiplyer #find the error of the longest series
  max_true<-true_slope+true_error #compute max and min values for slopes we are calling true
  min_true<-true_slope-true_error
  windows<-unique(test$N_years)#get a list of unique window lengths
  stability<-max(windows) #start with the assumption that the longest window is the only stable one
  for(i in 1:length(windows)){#for each window length, compute proportion 'correct'
    window_length<-windows[i]
    test_subset<-test[which(test$N_years==window_length),]
    number_of_windows<-nrow(test_subset)#how many windows
    correct_subset<-test_subset[which((test_subset$slope<max_true) & (test_subset$slope>min_true)),]
    number_of_correct<-nrow(correct_subset)#how many windows give the right answer
    percentage_correct<-100*number_of_correct/number_of_windows
    if(percentage_correct > min_percent){
      if(window_length < stability){
        stability<-window_length
      }
    }
  }
  return(stability)
}


#linear model
stability_time(EBS.dissim.simp[domain == "Full"], error_multiplyer = 1) #21
stability_time(EBS.dissim.simp[domain == "Inner"], error_multiplyer = 1) #22
stability_time(EBS.dissim.simp[domain == "Outer"], error_multiplyer = 1) #26
stability_time(EBS.dissim.simp[domain == "Middle"], error_multiplyer = 1) #22

#theil
stability_time(EBS.dissim.simp[domain == "Full"], error_multiplyer = 1, linear_model = "theil_sen_regression") #21
stability_time(EBS.dissim.simp[domain == "Inner"], error_multiplyer = 1, linear_model = "theil_sen_regression") #21
stability_time(EBS.dissim.simp[domain == "Outer"], error_multiplyer = 1, linear_model = "theil_sen_regression") #26
stability_time(EBS.dissim.simp[domain == "Middle"], error_multiplyer = 1, linear_model = "theil_sen_regression") #22


#21-26 years of observation in EBS to get reliable dissimilarity trend

#now a function that finds the absolute range of findings, and the absolute 
#range of significant findings

abs_range<- function(data, only_significant=FALSE, significance=0.05, linear_model = "lm"){#returns a two unit vector with the max and min slopes
  test<-multiple_breakups(data, linear_model = linear_model)
  if(only_significant== TRUE){ #if user specifies only significant values wanted, pull those
    test1<-test[which(test$p_value<significance),]
  }else{
    test1<-test
  }
  max_slope<-max(test1$slope)
  min_slope<-min(test1$slope)
  sloperange<-c(min_slope, max_slope)
  return(sloperange)
  
}

#and try it out

#linear model
abs_range(EBS.dissim.simp[domain == "Full"], only_significant = F, significance = 0.05)
abs_range(EBS.dissim.simp[domain == "Inner"], only_significant = F, significance = 0.05)
abs_range(EBS.dissim.simp[domain == "Outer"], only_significant = F, significance = 0.05)
abs_range(EBS.dissim.simp[domain == "Middle"], only_significant = F, significance = 0.05)

#theil
abs_range(EBS.dissim.simp[domain == "Full"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
abs_range(EBS.dissim.simp[domain == "Inner"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
abs_range(EBS.dissim.simp[domain == "Outer"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
abs_range(EBS.dissim.simp[domain == "Middle"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")





#now we want to find the absolute over and under estimate compared to the slope of the 
#longest series

relative_range<- function(data, only_significant=FALSE, significance=0.05, linear_model = "lm"){#returns a two unit vector with the max and min slopes
  test<-multiple_breakups(data, linear_model = linear_model)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  if(only_significant== TRUE){ #if user specifies only significant values wanted, pull those
    test1<-test[which(test$p_value<significance),]
  }else{
    test1<-test
  }
  max_slope<-max(test1$slope)-true_slope
  min_slope<-min(test1$slope)-true_slope
  sloperange<-c(min_slope, max_slope)
  return(sloperange)
  
}


#linear model
relative_range(EBS.dissim.simp[domain == "Full"], only_significant = F, significance = 0.05)
relative_range(EBS.dissim.simp[domain == "Inner"], only_significant = F, significance = 0.05)
relative_range(EBS.dissim.simp[domain == "Outer"], only_significant = F, significance = 0.05)
relative_range(EBS.dissim.simp[domain == "Middle"], only_significant = F, significance = 0.05)

#theil
relative_range(EBS.dissim.simp[domain == "Full"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
relative_range(EBS.dissim.simp[domain == "Inner"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
relative_range(EBS.dissim.simp[domain == "Outer"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
relative_range(EBS.dissim.simp[domain == "Middle"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")



relative_range_after_stability<- function(data, only_significant=FALSE, significance=0.05, linear_model = "lm"){#returns a two unit vector with the max and min slopes
  test<-multiple_breakups(data, linear_model = linear_model)
  stime<-stability_time(data)
  stest<-test[which(test$N_years>=stime),]
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  if(only_significant== TRUE){ #if user specifies only significant values wanted, pull those
    test1<-stest[which(test$p_value<significance),]
  }else{
    test1<-stest
  }
  max_slope<-max(test1$slope)-true_slope
  min_slope<-min(test1$slope)-true_slope
  sloperange<-c(min_slope, max_slope)
  return(sloperange)
  
}

##returns a two unit vector with the max and min slopes AFTER stability is reached
#linear model
relative_range_after_stability(EBS.dissim.simp[domain == "Full"], only_significant = F, significance = 0.05)
relative_range_after_stability(EBS.dissim.simp[domain == "Inner"], only_significant = F, significance = 0.05)
relative_range_after_stability(EBS.dissim.simp[domain == "Outer"], only_significant = F, significance = 0.05)
relative_range_after_stability(EBS.dissim.simp[domain == "Middle"], only_significant = F, significance = 0.05)

#theil
relative_range_after_stability(EBS.dissim.simp[domain == "Full"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
relative_range_after_stability(EBS.dissim.simp[domain == "Inner"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
relative_range_after_stability(EBS.dissim.simp[domain == "Outer"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")
relative_range_after_stability(EBS.dissim.simp[domain == "Middle"], only_significant = F, significance = 0.05, linear_model = "theil_sen_regression")



#proportion significant- finds the proportion of total windows with statistically significant values

proportion_significant<- function(data, significance=0.05, linear_model = "lm"){#returns a single value between 0 and 1
  test<-multiple_breakups(data, linear_model = linear_model)
  count<-nrow(test)
  significant_regressions<-test[which(test$p_value<significance),]
  count_sig<-nrow(significant_regressions)
  proportion<-count_sig/count
  return(proportion)
  
}

#linear model
proportion_significant(EBS.dissim.simp[domain == "Full"], significance = 0.05)
proportion_significant(EBS.dissim.simp[domain == "Inner"], significance = 0.05)
proportion_significant(EBS.dissim.simp[domain == "Outer"], significance = 0.05)
proportion_significant(EBS.dissim.simp[domain == "Middle"], significance = 0.05)

#theil
proportion_significant(EBS.dissim.simp[domain == "Full"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_significant(EBS.dissim.simp[domain == "Inner"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_significant(EBS.dissim.simp[domain == "Outer"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_significant(EBS.dissim.simp[domain == "Middle"], significance = 0.05, linear_model = "theil_sen_regression")

#proportion significantly wrong- we're going to define this as 'directionally wrong'
#where there is a significant relationship that does not match the direction of the true slope


proportion_wrong <- function(data, significance=0.05, linear_model = "lm"){#returns a single value between 0 and 1
  test<-multiple_breakups(data, linear_model = linear_model)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  true_p<-test[count,6]
  #case 1: true slope is not significant
  if (true_p>significance){
    wrong_windows<-test[which(test$p_value<significance),]
  }else{ #true slope is significant
    if(true_slope>0){#true slope is positive
      wrong_windows<-test[which(test$slope<0|test$p_value>significance),]#wrong means the slope is the wrong sign or 0
    }else{#true slope is negative
      wrong_windows<-test[which(test$slope>0|test$p_value>significance),]#wrong means the slope is the wrong sign or 0
    }
  }
  count_wrong<-nrow(wrong_windows)
  proportion<-count_wrong/count
  return(proportion)
  
}

#linear model
proportion_wrong(EBS.dissim.simp[domain == "Full"], significance = 0.05)
proportion_wrong(EBS.dissim.simp[domain == "Inner"], significance = 0.05)
proportion_wrong(EBS.dissim.simp[domain == "Outer"], significance = 0.05)
proportion_wrong(EBS.dissim.simp[domain == "Middle"], significance = 0.05)

#theil
proportion_wrong(EBS.dissim.simp[domain == "Full"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_wrong(EBS.dissim.simp[domain == "Inner"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_wrong(EBS.dissim.simp[domain == "Outer"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_wrong(EBS.dissim.simp[domain == "Middle"], significance = 0.05, linear_model = "theil_sen_regression")


#14% - 38% of the time, the slope is "significant" but is not the same as the "true" slope for 36 years

#proportion wrong by series length- basically the same thing as proportion wrong but looped 
#over all the unique window lengths. Will output a data frame with a window length and proportion
#of outputs are significantly misleading, plus average r square for that window length

proportion_wrong_series<- function(data, significance=0.05, linear_model = "lm"){#returns a single value between 0 and 1
  test<-multiple_breakups(data, linear_model = linear_model)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  true_p<-test[count,6]
  windows<-unique(test$N_years)#get a list of unique window lengths
  prop.vec<-c()#create a blank vector to store proportions in
  r.vec<-c()#create vector for storing average rsquare values in
  for(i in 1:length(windows)){#for each window length, compute proportion 'wrong'
    window_length<-windows[i]
    test_subset<-test[which(test$N_years==window_length),]
    number_of_windows<-nrow(test_subset)#how many windows
    #case 1: true slope is not significant
    if (true_p>significance){
      wrong_windows<-test_subset[which(test_subset$p_value<significance),]
    }else{ #true slope is significant
      if(true_slope>0){#true slope is positive
        wrong_windows<-test_subset[which(test_subset$slope<0|test_subset$p_value>significance),]#wrong means the slope is the wrong sign or 0
      }else{#true slope is negative
        wrong_windows<-test_subset[which(test_subset$slope>0|test_subset$p_value>significance),]#wrong means the slope is the wrong sign or 0
      }
    }
    count_wrong<-nrow(wrong_windows)
    proportion<-count_wrong/number_of_windows
    prop.vec<-c(prop.vec, proportion)
    avg.confidence<-mean(test_subset$r_square)
    r.vec<-c(r.vec, avg.confidence)
  }
  
  x_name <- "window_length"
  y_name <- "proportion_wrong"
  z_name <- "avg_r_square"
  
  dt <- data.table(windows, prop.vec, r.vec)
  names(dt) <- c(x_name, y_name, z_name)
  return(dt)
  
}


#test it

#linear model

proportion_wrong_series_full_lm <- proportion_wrong_series(EBS.dissim.simp[domain == "Full"], significance = 0.05)
proportion_wrong_series_full_lm[,domain := "Full"][,linear_model := "lm"]
proportion_wrong_series_inner_lm <- proportion_wrong_series(EBS.dissim.simp[domain == "Inner"], significance = 0.05)
proportion_wrong_series_inner_lm[,domain := "Inner"][,linear_model := "lm"]
proportion_wrong_series_outer_lm <- proportion_wrong_series(EBS.dissim.simp[domain == "Outer"], significance = 0.05)
proportion_wrong_series_outer_lm[,domain := "Outer"][,linear_model := "lm"]
proportion_wrong_series_middle_lm <- proportion_wrong_series(EBS.dissim.simp[domain == "Middle"], significance = 0.05)
proportion_wrong_series_middle_lm[,domain := "Middle"][,linear_model := "lm"]

#theil

proportion_wrong_series_full_theil_sen <- proportion_wrong_series(EBS.dissim.simp[domain == "Full"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_wrong_series_full_theil_sen[,domain := "Full"][,linear_model := "theil_sen"]
proportion_wrong_series_inner_theil_sen <- proportion_wrong_series(EBS.dissim.simp[domain == "Inner"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_wrong_series_inner_theil_sen[,domain := "Inner"][,linear_model := "theil_sen"]
proportion_wrong_series_outer_theil_sen <- proportion_wrong_series(EBS.dissim.simp[domain == "Outer"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_wrong_series_outer_theil_sen[,domain := "Outer"][,linear_model := "theil_sen"]
proportion_wrong_series_middle_theil_sen <- proportion_wrong_series(EBS.dissim.simp[domain == "Middle"], significance = 0.05, linear_model = "theil_sen_regression")
proportion_wrong_series_middle_theil_sen[,domain := "Middle"][,linear_model := "theil_sen"]
#proportion significantly wrong under stability time- we're going to define this as 'directionally wrong'
#where there is a significant relationship that does not match the direction of the true slope

proportion_wrong_series_all<-rbind(proportion_wrong_series_full_lm,
                                   proportion_wrong_series_full_theil_sen, proportion_wrong_series_inner_lm,
                                   proportion_wrong_series_inner_theil_sen, proportion_wrong_series_middle_lm,
                                   proportion_wrong_series_middle_theil_sen, proportion_wrong_series_outer_lm,
                                   proportion_wrong_series_outer_theil_sen)

#proportion wrong before stability
proportion_wrong_before_stability<- function(data, significance=0.05,
                                             min_percent=95, error_multiplyer=1,
                                             linear_model = "lm"){#returns a single value between 0 and 1
  
  test<-multiple_breakups(data, linear_model = linear_model)
  count<-nrow(test)
  true_slope<-test[count,4] #find the slope of the longest series
  true_p<-test[count,6]
  
  #cut out data below threshold
  threshold<-stability_time(data, min_percent, error_multiplyer)#find stability threshold
  test1<-test[which(test$N_years<threshold),]
  count1<-nrow(test1)
  #case 1: true slope is not significant
  if (true_p>significance){
    wrong_windows<-test1[which(test1$p_value<significance),]
  }else{ #true slope is significant
    if(true_slope>0){#true slope is positive
      wrong_windows<-test1[which(test1$slope<0|test1$p_value>significance),]#wrong means the slope is the wrong sign or 0
    }else{#true slope is negative
      wrong_windows<-test1[which(test1$slope>0|test1$p_value>significance),]#wrong means the slope is the wrong sign or 0
    }
  }
  count_wrong<-nrow(wrong_windows)
  proportion<-count_wrong/count1
  return(proportion)
  
}

#linear model
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Full"], significance = 0.05) #27%
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Inner"], significance = 0.05) #14%
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Outer"], significance = 0.05) #39%
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Middle"], significance = 0.05) #33%

#theil
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Full"], significance = 0.05, linear_model = "theil_sen_regression") #26%
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Inner"], significance = 0.05, linear_model = "theil_sen_regression") #14%
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Outer"], significance = 0.05, linear_model = "theil_sen_regression") #35%
proportion_wrong_before_stability(EBS.dissim.simp[domain == "Middle"], significance = 0.05, linear_model = "theil_sen_regression") #32%


#implement another charting function that gives the proportion wrong by window length

wrongness_plot<-function(data, significance=0.05, min_percent=95, error_multiplyer=1, title ="", linear_model = "lm"){
  threshold<-stability_time(data, min_percent, error_multiplyer, linear_model = linear_model)#find stability threshold
  wrongness<-proportion_wrong_series(data, significance)
  maxyears<-max(wrongness$window_length)
  plot<- ggplot(wrongness) +
    theme_classic() +
    geom_vline(xintercept = (threshold-0.1), linetype = 3, color="grey38") +
    geom_smooth(aes(y = proportion_wrong, x = window_length, 
                    linetype="Propwrong", color="Propwrong"), se=FALSE)+
    geom_point(aes(y = proportion_wrong, x = window_length, 
                   shape="Propwrong", fill="Propwrong"), size=3)+
    geom_smooth(aes(y = avg_r_square, x = window_length, 
                    linetype="rsq", color="rsq"), se=FALSE)+
    geom_point(aes(y = avg_r_square, x = window_length, 
                   shape="rsq", fill="rsq"), size=3)+
    scale_fill_manual(name="", values=c(Propwrong="black",rsq="orange"),
                      labels=c("Proportion\n wrong", expression("Average R"^2)))+
    scale_shape_manual(name="", values=c(Propwrong=21, rsq=24), 
                       labels=c("Proportion\n wrong", expression("Average R"^2)))+
    scale_linetype_manual(name="", values=c(Propwrong=1, rsq=2), 
                          labels=c("Proportion\n wrong", expression("Average R"^2)))+
    scale_color_manual(name="", values=c(Propwrong="blue", rsq="red"), 
                       labels=c("Proportion\n wrong", expression("Average R"^2)))+
    ggtitle(title)+
    xlab("Number of years in window")+
    ylab("Average value")+
    ylim(0,1)
  return(plot)
}


#plot
#linear model
wrongness_plot_full_lm <- wrongness_plot(EBS.dissim.simp[domain == "Full"], title="Full EBS, lm()", significance=0.05)
ggsave(wrongness_plot_full_lm, path = file.path("Figures"), filename = "wrongness_plot_full_lm.jpg")
wrongness_plot_inner_lm <- wrongness_plot(EBS.dissim.simp[domain == "Inner"], title="Inner EBS, lm()", significance=0.05)
ggsave(wrongness_plot_inner_lm, path = file.path("Figures"), filename = "wrongness_plot_inner_lm.jpg")
wrongness_plot_outer_lm <- wrongness_plot(EBS.dissim.simp[domain == "Outer"], title="Outer EBS, lm()", significance=0.05)
ggsave(wrongness_plot_outer_lm, path = file.path("Figures"), filename = "wrongness_plot_outer_lm.jpg")
wrongness_plot_middle_lm <- wrongness_plot(EBS.dissim.simp[domain == "Middle"], title="Middle EBS, lm()", significance=0.05)
ggsave(wrongness_plot_middle_lm, path = file.path("Figures"), filename = "wrongness_plot_middle_lm.jpg")

#theil sein
wrongness_plot_full_theil_sen <- wrongness_plot(EBS.dissim.simp[domain == "Full"], linear_model = "theil_sen_regression", title="Full EBS, theil sen", significance=0.05)
ggsave(wrongness_plot_full_theil_sen, path = file.path("Figures"), filename = "wrongness_plot_full_theil_sen.jpg")
wrongness_plot_inner_theil_sen <- wrongness_plot(EBS.dissim.simp[domain == "Inner"], linear_model = "theil_sen_regression", title="Inner EBS, theil sen", significance=0.05)
ggsave(wrongness_plot_inner_theil_sen, path = file.path("Figures"), filename = "wrongness_plot_inner_theil_sen.jpg")
wrongness_plot_outer_theil_sen <- wrongness_plot(EBS.dissim.simp[domain == "Outer"], linear_model = "theil_sen_regression", title="Outer EBS, theil sen", significance=0.05)
ggsave(wrongness_plot_outer_theil_sen, path = file.path("Figures"), filename = "wrongness_plot_outer_theil_sen.jpg")
wrongness_plot_middle_theil_sen <- wrongness_plot(EBS.dissim.simp[domain == "Middle"], linear_model = "theil_sen_regression", title="Middle EBS, theil sen", significance=0.05)
ggsave(wrongness_plot_middle_theil_sen, path = file.path("Figures"), filename = "wrongness_plot_middle_theil_sen.jpg")

#merge into one figure
wrongness_plot_merge <- plot_grid(
  wrongness_plot_full_lm  +  theme(plot.title = element_text(face = "bold", color = "black"), legend.key.size =unit(0.5, unit = "cm"), legend.position = c(0.3, 0.7), axis.title.x = element_blank(), text = element_text(size = 14)), 
  wrongness_plot_inner_lm  + theme(plot.title = element_text(face = "bold", color = "#AA4499"), legend.key.size =unit(0.5, unit = "cm"), legend.position = "null",  text = element_text(size = 14), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()), 
  wrongness_plot_middle_lm  + theme(plot.title = element_text(face = "bold", color = "#44AA99"), legend.key.size =unit(0.5, unit = "cm"), legend.position = "null",  text = element_text(size = 14), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()), 
  wrongness_plot_outer_lm  + theme(plot.title = element_text(face = "bold", color = "#999933"), legend.key.size =unit(0.5, unit = "cm"), legend.position = "null",  text = element_text(size = 14), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()), 
  wrongness_plot_full_theil_sen  + theme(plot.title = element_text(face = "bold", color = "black"), legend.key.size =unit(0.5, unit = "cm"), legend.position = "null", text = element_text(size = 14)), 
  wrongness_plot_inner_theil_sen  + theme(plot.title = element_text(face = "bold", color = "#AA4499"), legend.key.size =unit(0.5, unit = "cm"), legend.position = "null",  text = element_text(size = 14), axis.text.y = element_blank(), axis.title.y = element_blank()), 
  wrongness_plot_middle_theil_sen  + theme(plot.title = element_text(face = "bold", color = "#44AA99"), legend.key.size =unit(0.5, unit = "cm"), legend.position = "null",  text = element_text(size = 14), axis.text.y = element_blank(), axis.title.y = element_blank()),
  wrongness_plot_outer_theil_sen  + theme(plot.title = element_text(face = "bold", color = "#999933"), legend.key.size =unit(0.5, unit = "cm"), legend.position = "null",  text = element_text(size = 14), axis.text.y = element_blank(), axis.title.y = element_blank()), 
  ncol = 4, nrow = 2)

ggsave(wrongness_plot_merge, path = file.path("Figures"),
       filename = "wrongness_plot_merge.jpg", height = 8, width = 15, unit = "in")



#now for a function that plots all the lines by window length

broken_stick_plot<-function(data, title="", significance=0.05, window_length=3, linear_model = "lm"){
  out<-multiple_breakups(data, linear_model = linear_model)
  years<-length(unique(out$start_year))
  maxyears<-max(out$N_years)
  count<-nrow(out)
  #compute mean of longest series
  true_slope<-out[count,4] #find the slope of the longest series
  true_intercept<-(out[count,7]) #find the intercept of the longest series
  out<-out[which(out$N_years==window_length),] #only work with one window length per plot
  #create a separate frame for significant and not results
  out_sig<-out[which(out$p_value<significance),]
  countsig<-nrow(out_sig)#count the number of rows in the set we want to plot
  out_not<-out[which(out$p_value>significance),]
  countnot<-nrow(out_not)#count the number of rows in the set we want to plot
  plot<- ggplot(data, aes(x=year, y=bray_curtis_dissimilarity_balanced_mean)) +
    theme_classic()+geom_smooth(linetype=0, fill="lightblue1", method=lm, formula='y ~ x', 
                                level=0.99)#99% confidence interval around longest series
  if(countnot>0){
    for(i in 1:countnot){ #plot not significant windows
      slopei<-out_not$slope[i]
      intercepti<-out_not$intercept[i]
      plot<-plot+geom_abline(slope=slopei, intercept=intercepti, linetype=3, colour="grey12")
    }
  }
  if(countsig>0){
    for(i in 1:countsig){ #plot significant windows
      slopei<-out_sig$slope[i]
      intercepti<-out_sig$intercept[i]
      plot<-plot+geom_abline(slope=slopei, intercept=intercepti, linetype=2, colour="red")
    }
  }
  
  
  plot<-plot+ ggtitle(title)+
    geom_abline(slope=true_slope, intercept=true_intercept, linetype=1, colour="grey16", size=1)+
    geom_point(size=3, pch=21, fill="grey22")+
    xlab("Year")+ylab("β diversity")
  return(plot)
}
#test it
broken_stick_plot(EBS.dissim.simp, window_length = 32, significance = 0.5)






#linear model
#window length 3
broken_stick_plot_w3_full_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], window_length = 3, title="Full EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w3_full_lm, path = file.path("Figures"), filename = "broken_stick_plot_w3_full_lm.jpg")
broken_stick_plot_w3_inner_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], window_length = 3, title="Inner EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w3_inner_lm, path = file.path("Figures"), filename = "broken_stick_plot_w3_inner_lm.jpg")
broken_stick_plot_w3_outer_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], window_length = 3, title="Outer EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w3_outer_lm, path = file.path("Figures"), filename = "broken_stick_plot_w3_outer_lm.jpg")
broken_stick_plot_w3_middle_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], window_length = 3, title="Middle EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w3_middle_lm, path = file.path("Figures"), filename = "broken_stick_plot_w3_middle_lm.jpg")

#window length 10
broken_stick_plot_w10_full_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], window_length = 10, title="Full EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w10_full_lm, path = file.path("Figures"), filename = "broken_stick_plot_w10_full_lm.jpg")
broken_stick_plot_w10_inner_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], window_length = 10, title="Inner EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w10_inner_lm, path = file.path("Figures"), filename = "broken_stick_plot_w10_inner_lm.jpg")
broken_stick_plot_w10_outer_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], window_length = 10, title="Outer EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w10_outer_lm, path = file.path("Figures"), filename = "broken_stick_plot_w10_outer_lm.jpg")
broken_stick_plot_w10_middle_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], window_length = 10, title="Middle EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w10_middle_lm, path = file.path("Figures"), filename = "broken_stick_plot_w10_middle_lm.jpg")

#window length 20
broken_stick_plot_w20_full_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], window_length = 20, title="Full EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w20_full_lm, path = file.path("Figures"), filename = "broken_stick_plot_w20_full_lm.jpg")
broken_stick_plot_w20_inner_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], window_length = 20, title="Inner EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w20_inner_lm, path = file.path("Figures"), filename = "broken_stick_plot_w20_inner_lm.jpg")
broken_stick_plot_w20_outer_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], window_length = 20, title="Outer EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w20_outer_lm, path = file.path("Figures"), filename = "broken_stick_plot_w20_outer_lm.jpg")
broken_stick_plot_w20_middle_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], window_length = 20, title="Middle EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w20_middle_lm, path = file.path("Figures"), filename = "broken_stick_plot_w20_middle_lm.jpg")

#window length 30
broken_stick_plot_w30_full_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], window_length = 30, title="Full EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w30_full_lm, path = file.path("Figures"), filename = "broken_stick_plot_w30_full_lm.jpg")
broken_stick_plot_w30_inner_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], window_length = 30, title="Inner EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w30_inner_lm, path = file.path("Figures"), filename = "broken_stick_plot_w30_inner_lm.jpg")
broken_stick_plot_w30_outer_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], window_length = 30, title="Outer EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w30_outer_lm, path = file.path("Figures"), filename = "broken_stick_plot_w30_outer_lm.jpg")
broken_stick_plot_w30_middle_lm <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], window_length = 30, title="Middle EBS, lm()", significance=0.05)
ggsave(broken_stick_plot_w30_middle_lm, path = file.path("Figures"), filename = "broken_stick_plot_w30_middle_lm.jpg")

#theil sein
#window length 3
broken_stick_plot_w3_full_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], window_length = 3, linear_model = "theil_sen_regression", title="Full EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w3_full_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w3_full_theil_sen.jpg")
broken_stick_plot_w3_inner_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], window_length = 3, linear_model = "theil_sen_regression", title="Inner EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w3_inner_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w3_inner_theil_sen.jpg")
broken_stick_plot_w3_outer_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], window_length = 3, linear_model = "theil_sen_regression", title="Outer EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w3_outer_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w3_outer_theil_sen.jpg")
broken_stick_plot_w3_middle_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], window_length = 3, linear_model = "theil_sen_regression", title="Middle EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w3_middle_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w3_middle_theil_sen.jpg")

#window length 10
broken_stick_plot_w10_full_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], linear_model = "theil_sen_regression", window_length = 10, title="Full EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w10_full_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w10_full_theil_sen.jpg")
broken_stick_plot_w10_inner_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], linear_model = "theil_sen_regression", window_length = 10, title="Inner EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w10_inner_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w10_inner_theil_sen.jpg")
broken_stick_plot_w10_outer_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], linear_model = "theil_sen_regression", window_length = 10, title="Outer EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w10_outer_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w10_outer_theil_sen.jpg")
broken_stick_plot_w10_middle_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], linear_model = "theil_sen_regression", window_length = 10, title="Middle EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w10_middle_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w10_middle_theil_sen.jpg")

#window length 20
broken_stick_plot_w20_full_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], linear_model = "theil_sen_regression", window_length = 20, title="Full EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w20_full_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w20_full_theil_sen.jpg")
broken_stick_plot_w20_inner_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], linear_model = "theil_sen_regression", window_length = 20, title="Inner EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w20_inner_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w20_inner_theil_sen.jpg")
broken_stick_plot_w20_outer_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], linear_model = "theil_sen_regression", window_length = 20, title="Outer EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w20_outer_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w20_outer_theil_sen.jpg")
broken_stick_plot_w20_middle_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], linear_model = "theil_sen_regression", window_length = 20, title="Middle EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w20_middle_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w20_middle_theil_sen.jpg")

#window length 30
broken_stick_plot_w30_full_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Full"], linear_model = "theil_sen_regression", window_length = 30, title="Full EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w30_full_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w30_full_theil_sen.jpg")
broken_stick_plot_w30_inner_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Inner"], linear_model = "theil_sen_regression", window_length = 30, title="Inner EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w30_inner_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w30_inner_theil_sen.jpg")
broken_stick_plot_w30_outer_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Outer"], linear_model = "theil_sen_regression", window_length = 30, title="Outer EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w30_outer_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w30_outer_theil_sen.jpg")
broken_stick_plot_w30_middle_theil_sen <- broken_stick_plot(EBS.dissim.simp[domain == "Middle"], linear_model = "theil_sen_regression", window_length = 30, title="Middle EBS, theil sen", significance=0.05)
ggsave(broken_stick_plot_w30_middle_theil_sen, path = file.path("Figures"), filename = "broken_stick_plot_w30_middle_theil_sen.jpg")

#merge into one figure
#lm
broken_stick_plot_merge_lm <- plot_grid(
  #top row
  broken_stick_plot_w3_full_lm + ggtitle("Full EBS") + ylab("3-year window\nβ diversity") + theme(plot.title = element_text(face = "bold", color = "black"), text = element_text(size = 20), axis.title.x = element_blank()),
  broken_stick_plot_w3_inner_lm + ggtitle("Inner EBS") + theme(plot.title = element_text(face = "bold", color = "#AA4499"), text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()),
  broken_stick_plot_w3_middle_lm + ggtitle("Middle EBS") + theme(plot.title = element_text(face = "bold", color = "#44AA99"), text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()),
  broken_stick_plot_w3_outer_lm + ggtitle("Outer EBS") + theme(plot.title = element_text(face = "bold", color = "#999933"), text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()),
  #second row
  broken_stick_plot_w10_full_lm + ylab("10-year window\nβ diversity") + theme(text = element_text(size = 20), axis.title.x = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w10_inner_lm + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w10_middle_lm + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w10_outer_lm + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  #third row
  broken_stick_plot_w20_full_lm + ylab("20-year window\nβ diversity") + theme(text = element_text(size = 20), axis.title.x = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w20_inner_lm + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w20_middle_lm + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w20_outer_lm + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  #fourth row
  broken_stick_plot_w30_full_lm + ylab("30-year window\nβ diversity") + theme(text = element_text(size = 20), plot.title = element_blank()),
  broken_stick_plot_w30_inner_lm + theme(text = element_text(size = 20), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w30_middle_lm + theme(text = element_text(size = 20), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w30_outer_lm + theme(text = element_text(size = 20), axis.title.y = element_blank(), plot.title = element_blank()),
  ncol = 4, nrow = 4)

ggsave(broken_stick_plot_merge_lm, path = file.path("Figures"),
       filename = "broken_stick_plot_merge_lm.jpg", height = 15, width = 22, unit = "in")




#theil sen

broken_stick_plot_merge_theil_sen <- plot_grid(
  #top row
  broken_stick_plot_w3_full_theil_sen + ggtitle("Full EBS") + ylab("3-year window\nβ diversity") + theme(plot.title = element_text(face = "bold", color = "black"), text = element_text(size = 20), axis.title.x = element_blank()),
  broken_stick_plot_w3_inner_theil_sen + ggtitle("Inner EBS") + theme(plot.title = element_text(face = "bold", color = "#AA4499"), text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()),
  broken_stick_plot_w3_middle_theil_sen + ggtitle("Middle EBS") + theme(plot.title = element_text(face = "bold", color = "#44AA99"), text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()),
  broken_stick_plot_w3_outer_theil_sen + ggtitle("Outer EBS") + theme(plot.title = element_text(face = "bold", color = "#999933"), text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank()),
  #second row
  broken_stick_plot_w10_full_theil_sen + ylab("10-year window\nβ diversity") + theme(text = element_text(size = 20), axis.title.x = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w10_inner_theil_sen + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w10_middle_theil_sen + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w10_outer_theil_sen + theme(text = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  #third row
  broken_stick_plot_w20_full_theil_sen + ylab("20-year window\nβ diversity") + theme(axis.title.x = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w20_inner_theil_sen + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w20_middle_theil_sen + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w20_outer_theil_sen + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank()),
  #fourth row
  broken_stick_plot_w30_full_theil_sen + ylab("30-year window\nβ diversity") + theme(plot.title = element_blank()),
  broken_stick_plot_w30_inner_theil_sen + theme(axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w30_middle_theil_sen + theme(axis.title.y = element_blank(), plot.title = element_blank()),
  broken_stick_plot_w30_outer_theil_sen + theme(axis.title.y = element_blank(), plot.title = element_blank()),
  ncol = 4, nrow = 4)

ggsave(broken_stick_plot_merge_theil_sen, path = file.path("Figures"),
       filename = "broken_stick_plot_merge_theil_sen.jpg", height = 15, width = 19, unit = "in")
