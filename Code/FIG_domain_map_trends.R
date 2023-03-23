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

#Pull grobs of map and trends (facet)

year_beta_bydomain <- readRDS(file.path("Figures","year_beta_bydomain.Rds"))
Alaska_domains <- readRDS(file.path("Figures","Alaska_domains.Rds"))

#Pull in distances dissimilarities csv to calculate slopes of year~dissimilarity
EBS.distances_dissimilarities_allyears <-  fread(file.path("Output","EBS.distances_dissimilarities_allyears.csv"))

#and then calculate slopes etc.
#slopes
slope_full <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_inner <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_middle <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_outer <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)

#R^2 values summary(model)
R2_full <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_inner <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_middle <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_outer <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)

#p_values summary(model)
p_value_full <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_inner <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_middle <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_outer <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)

#######################
##PLOT
#######################

#place model values on plot
map_year_beta_bydomain_annotate <- ggdraw(xlim = c(0,10), ylim = c(0,20)) +
  draw_plot(year_beta_bydomain, x = 0, y = 0, width = 9, height = 10) +
  draw_plot(Alaska_domains, x = 0.2, y = 10, width = 10, height = 10) +
  geom_text(aes(x = 2.8, y = 1.5, label = paste0("Slope = ",slope_middle, "  R^2 = ",R2_middle, "  p = ",p_value_middle)), size = 3) +
  geom_text(aes(x = 7, y = 1.5, label = paste0("Slope = ",slope_outer, "  R^2 = ",R2_outer, "  p = ",p_value_outer)), size = 3) +
  geom_text(aes(x = 2.8, y = 5.8, label = paste0("Slope = ",slope_full, "  R^2 = ",R2_full, "  p = ",p_value_full)), size = 3) +
  geom_text(aes(x = 7, y = 8.8, label = paste0("Slope = ",slope_inner, "  R^2 = ",R2_inner, "  p = ",p_value_inner)), size = 3)

map_year_beta_bydomain_annotate 

save_plot(map_year_beta_bydomain_annotate, filename = file.path("Figures","map_year_beta_bydomain_annotate.jpg"), base_height = 7, base_width = 6)
   