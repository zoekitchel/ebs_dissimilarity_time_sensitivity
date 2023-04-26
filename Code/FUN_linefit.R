###########################
#Linefit adapted by Zoë Kitchel from Bahlai et al. 2021 for Eastern Bering Sea β diversity (dissimilarity) Through Time

#######################
##VERSIONS##
#R 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
#macOS Big Sur 11.7

#######################
##PACKAGES##
#######################
library(RobustLinearReg)
#######################
##FUNCTION
#######################

# next we need a function that runs a simple linear model of x=year, y=response variable

#linear model can be "lm" or "theil_sen_regression" or other

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
