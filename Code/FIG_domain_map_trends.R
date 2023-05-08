###########################
#This script merges the map figure with the trend figure for year~dissimilarity for each domain
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
library(RobustLinearReg)

#######################
##DATA##
#######################

#Pull grobs of map and trends (facet)

#trends
year_beta_bydomain <- readRDS(
  file.path("Figures","year_beta_bydomain.Rds")) #bray curtis balanced
year_beta_bray_curtis_balanced_theilsenreg_bydomain <- readRDS(
  file.path("Figures","Supplement","Theil_Sen","year_beta_bray_curtis_balanced_theilsenreg_bydomain.Rds")) #theil sen regression
year_beta_bray_curtis_biomassgradient_bydomain <- readRDS(
  file.path("Figures","Supplement","Biomass_gradient","year_beta_bray_curtis_biomassgradient_bydomain.Rds")) #BC biomass gradient
year_beta_jaccard_bydomain <- readRDS(file.path("Figures","Supplement","Jaccard","year_beta_jaccard_bydomain.Rds")) #jaccard turnover

Alaska_domains <- readRDS(file.path("Figures","Alaska_domains.Rds")) #map

#pull in broken stick plots

  #balanced BC
  broken_stick_plot_w3_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w3_full_lm_fig1.Rds"))
  broken_stick_plot_w10_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w10_full_lm_fig1.Rds"))
  broken_stick_plot_w20_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w20_full_lm_fig1.Rds"))
  broken_stick_plot_w35_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w35_full_lm_fig1.Rds"))
  
  #Theil Sen
  broken_stick_plot_w3_full_theilsen_fig1 <- readRDS(file.path("Figures","Supplement", "Theil_Sen","broken_stick_plot_w3_full_theilsen_fig1.Rds"))
  broken_stick_plot_w10_full_theilsen_fig1 <- readRDS(file.path("Figures","Supplement", "Theil_Sen","broken_stick_plot_w10_full_theilsen_fig1.Rds"))
  broken_stick_plot_w20_full_theilsen_fig1 <- readRDS(file.path("Figures","Supplement", "Theil_Sen","broken_stick_plot_w20_full_theilsen_fig1.Rds"))
  broken_stick_plot_w35_full_theilsen_fig1 <- readRDS(file.path("Figures","Supplement", "Theil_Sen","broken_stick_plot_w35_full_theilsen_fig1.Rds"))
  
  #biomass gradient BC
  broken_stick_plot_w3_full_biomass_gradient_fig1 <- readRDS(file.path("Figures", "Supplement", "Biomass_gradient","broken_stick_plot_w3_full_biomass_gradient_fig1.Rds"))
  broken_stick_plot_w10_full_biomass_gradient_fig1 <- readRDS(file.path("Figures", "Supplement", "Biomass_gradient","broken_stick_plot_w10_full_biomass_gradient_fig1.Rds"))
  broken_stick_plot_w20_full_biomass_gradient_fig1 <- readRDS(file.path("Figures", "Supplement", "Biomass_gradient","broken_stick_plot_w20_full_biomass_gradient_fig1.Rds"))
  broken_stick_plot_w35_full_biomass_gradient_fig1 <- readRDS(file.path("Figures", "Supplement", "Biomass_gradient","broken_stick_plot_w35_full_biomass_gradient_fig1.Rds"))
  
  #jaccard turnover
  broken_stick_plot_w3_full_jaccard_turnover_fig1 <- readRDS(file.path("Figures", "Supplement", "Jaccard","broken_stick_plot_w3_full_jaccard_turnover_fig1.Rds"))
  broken_stick_plot_w10_full_jaccard_turnover_fig1 <- readRDS(file.path("Figures", "Supplement", "Jaccard","broken_stick_plot_w10_full_jaccard_turnover_fig1.Rds"))
  broken_stick_plot_w20_full_jaccard_turnover_fig1 <- readRDS(file.path("Figures", "Supplement", "Jaccard","broken_stick_plot_w20_full_jaccard_turnover_fig1.Rds"))
  broken_stick_plot_w35_full_jaccard_turnover_fig1 <- readRDS(file.path("Figures", "Supplement", "Jaccard","broken_stick_plot_w35_full_jaccard_turnover_fig1.Rds"))
  

#Pull in distances dissimilarities csv to calculate slopes of year~dissimilarity
EBS.distances_dissimilarities_allyears <-  fread(file.path("Output","EBS.distances_dissimilarities_allyears.csv"))
  #reorder factors
  EBS.distances_dissimilarities_allyears[,Domain := factor(domain, levels = c("Full","Outer","Middle","Inner"))]


#also take a look at no pollock
EBS.distances_dissimilarities_allyears_nopollock <-  fread(file.path("Output","EBS.distances_dissimilarities_allyears_nopollock.csv"))

#and then calculate slopes etc.
#slopes
slope_full <- signif(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],2)
slope_inner <- signif(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],2)
slope_middle <- signif(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],2)
slope_outer <- signif(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],2)

slope_se_full <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year,
                                   data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,2],2)
slope_se_inner <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,2],2)
slope_se_middle <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,2],2)
slope_se_outer <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,2],2)

#R^2 values summary(model)
R2_full <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_inner <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_middle <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_outer <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)

#p_values summary(model)
p_value_full <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_inner <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_middle <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_outer <- signif(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)


#THEIL SEN
#slopes
slope_full_theil <- signif(coefficients(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],2)

slope_se_full_theil <- signif(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year,
                                   data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,2],2)
#R^2 values summary(model)
R2_full_theil <- signif(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)

#p_values summary(model)
p_value_full_theil <- signif(summary(theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)


#JACCARD (same)

#slopes
slope_jaccard_full <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_jaccard_inner <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_jaccard_middle <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_jaccard_outer <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)

#SE of slope
slope_se_jaccard_full <- signif(summary(lm(jaccard_dissimilarity_turnover_mean~year,
                                   data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,2],2)
slope_se_jaccard_inner <- signif(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,2],2)
slope_se_jaccard_middle <- signif(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,2],2)
slope_se_jaccard_outer <- signif(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,2],2)



#R^2 values summary(model)
R2_jaccard_full <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_jaccard_inner <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_jaccard_middle <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_jaccard_outer <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)


#p_values summary(model)
p_value_jaccard_full <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_jaccard_inner <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_jaccard_middle <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_jaccard_outer <- round(summary(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)

####Biomass gradient
#slopes
slope_bray_curtis_full_biomassgradient <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_bray_curtis_inner_biomassgradient <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_bray_curtis_middle_biomassgradient <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_bray_curtis_outer_biomassgradient <- round(coefficients(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)

#SE of slope
slope_se_bray_curtis_full_biomassgradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year,
                                           data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,2],2)
slope_se_bray_curtis_inner_biomassgradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,2],2)
slope_se_bray_curtis_middle_biomassgradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,2],2)
slope_se_bray_curtis_outer_biomassgradient <- signif(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,2],2)



#R^2 values summary(model)
R2_bray_curtis_full_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$adj.r.squared,2)
R2_bray_curtis_inner_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$adj.r.squared,2)
R2_bray_curtis_middle_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$adj.r.squared,2)
R2_bray_curtis_outer_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$adj.r.squared,2)


#p_values summary(model)
p_value_bray_curtis_full_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))$coefficients[2,4],2)
p_value_bray_curtis_inner_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))$coefficients[2,4],2)
p_value_bray_curtis_middle_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))$coefficients[2,4],2)
p_value_bray_curtis_outer_biomassgradient <- round(summary(lm(bray_curtis_dissimilarity_gradient_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))$coefficients[2,4],2)



#NO POLLOCK (same)
#slopes
slope_full_nopollock <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Full",]))[[2]],5)
slope_inner_nopollock <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Inner",]))[[2]],5)
slope_middle_nopollock <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Middle",]))[[2]],5)
slope_outer_nopollock <- round(coefficients(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Outer",]))[[2]],5)

#R^2 values summary(model)
R2_full_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Full",]))$adj.r.squared,2)
R2_inner_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Inner",]))$adj.r.squared,2)
R2_middle_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Middle",]))$adj.r.squared,2)
R2_outer_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Outer",]))$adj.r.squared,2)

#p_values summary(model)
p_value_full_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Full",]))$coefficients[2,4],2)
p_value_inner_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Inner",]))$coefficients[2,4],2)
p_value_middle_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Middle",]))$coefficients[2,4],2)
p_value_outer_nopollock <- round(summary(lm(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears_nopollock[domain == "Outer",]))$coefficients[2,4],2)

#######################
##PLOT
#######################
#Bray Curtis balaned with Pollock
#outer, middle, inner
BC_dissim_O_M_I <- ggplot(EBS.distances_dissimilarities_allyears[Domain != "Full",]) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean, color = Domain), size = 1) +
  labs(color = "Domain", x = "Year", y = "β diversity") +
  geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean), color = "black", method = "lm", se = F, linetype = "longdash", linewidth = 0.6) +
  scale_color_manual(values = c("#999933", "#44AA99","#AA4499")) +
  facet_wrap(~Domain) +
  theme_classic() +
  theme(legend.position = "null", strip.text = element_blank())

#full
BC_dissim_F <- ggplot(EBS.distances_dissimilarities_allyears[Domain == "Full",]) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean)) +
  labs(x = "Year", y = "β diversity") +
  geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean), color = "black",
              method = "lm", se = F, linetype = "longdash") +
#  geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean), color = "darkgrey",
         #     method = "gam", se = F, linetype = "dotted") +
  theme_classic() +
  theme(legend.position = "null")

#place model values on plot
map_year_beta_bydomain_annotate <- ggdraw(xlim = c(0,30), ylim = c(0,20)) +
  draw_plot(BC_dissim_O_M_I + theme(axis.title.x = element_blank()), x = 0, y = 0, width = 10, height = 8) +
  draw_plot(Alaska_domains, x = 0, y =6.5, width = 10, height = 15) +
  draw_plot(BC_dissim_F + theme(axis.title.x = element_blank()), x = 10, y = 0, width = 13, height = 20) +
  geom_text(aes(x = 13, y = 18, label = paste0("slope = ",slope_full,"+/-",slope_se_full,"\np = ",p_value_full)), size = 3) +
  draw_plot(broken_stick_plot_w3_full_lm_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 15,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w10_full_lm_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 10,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w20_full_lm_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 5,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w35_full_lm_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 0,  width = 7,  height = 5) +
  geom_text(aes(x = 0.5, y = 19.7, label = "a."), size = 4, fontface = "bold") +
  geom_text(aes(x = 1.3, y = 7.3, label = "b."), size = 4, fontface = "bold") +
  geom_text(aes(x = 4.3, y =7.3, label = "c."), size = 4, fontface = "bold") +
  geom_text(aes(x = 7.3, y =7.3, label = "d."), size = 4, fontface = "bold") +
  geom_text(aes(x = 10.5, y = 19.7, label = "e."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 19.7, label = "f."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 14.5, label = "g."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 9.5, label = "h."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 4.5, label = "i."), size = 4, fontface = "bold")
  
 
map_year_beta_bydomain_annotate 

ggsave(map_year_beta_bydomain_annotate, path = file.path("Figures"), filename = "map_year_beta_bydomain_annotate.jpg", height =5, width = 13.5)

#################
##SUPPLEMENT
#################

full_theilsen_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Full",])
inner_theilsen_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Inner",])
middle_theilsen_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Middle",])
outer_theilsen_mod <- theil_sen_regression(bray_curtis_dissimilarity_balanced_mean~year, data = EBS.distances_dissimilarities_allyears[Domain == "Outer",])



EBS.distances_dissimilarities_allyears[Domain == "Full",pred_dissim := predict(full_theilsen_mod)]
EBS.distances_dissimilarities_allyears[Domain == "Inner",pred_dissim := predict(inner_theilsen_mod)]
EBS.distances_dissimilarities_allyears[Domain == "Middle",pred_dissim := predict(middle_theilsen_mod)]
EBS.distances_dissimilarities_allyears[Domain == "Outer",pred_dissim := predict(outer_theilsen_mod)]

#THEILSEN
#outer, middle, inner
BC_dissim_O_M_I_theilsen <- ggplot(EBS.distances_dissimilarities_allyears[Domain != "Full",]) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean, color = Domain), size = 1) +
  labs(color = "Domain", x = "Year", y = "β diversity") +
  geom_line(aes(x = year, y = pred_dissim), color = "black", linetype = "longdash", linewidth = 0.6) +
  scale_color_manual(values = c("#999933", "#44AA99","#AA4499")) +
  facet_wrap(~Domain) +
  theme_classic() +
  theme(legend.position = "null", strip.text = element_blank())

#full
BC_dissim_F_theilsen <- ggplot(EBS.distances_dissimilarities_allyears[Domain == "Full",]) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean)) +
  labs(x = "Year", y = "β diversity") +
  geom_line(aes(x = year, y = pred_dissim), color = "black", linetype = "longdash", linewidth = 1) +
  #  geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_balanced_mean), color = "darkgrey",
  #     method = "gam", se = F, linetype = "dotted") +
  theme_classic() +
  theme(legend.position = "null")

#place model values on plot
map_year_beta_bydomain_annotate_theilsen <- ggdraw(xlim = c(0,30), ylim = c(0,20)) +
  draw_plot(BC_dissim_O_M_I_theilsen + theme(axis.title.x = element_blank()), x = 0, y = 0, width = 10, height = 8) +
  draw_plot(Alaska_domains, x = 0, y =6.5, width = 10, height = 15) +
  draw_plot(BC_dissim_F_theilsen + theme(axis.title.x = element_blank()), x = 10, y = 0, width = 13, height = 20) +
  geom_text(aes(x = 13, y = 18, label = paste0("slope = ",slope_full_theil,"+/-",slope_se_full_theil,"\np = ",p_value_full_theil)), size = 3) +
  draw_plot(broken_stick_plot_w3_full_theilsen_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 15,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w10_full_theilsen_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 10,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w20_full_theilsen_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 5,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w35_full_theilsen_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 0,  width = 7,  height = 5) +
  geom_text(aes(x = 0.5, y = 19.7, label = "a."), size = 4, fontface = "bold") +
  geom_text(aes(x = 1.3, y = 7.3, label = "b."), size = 4, fontface = "bold") +
  geom_text(aes(x = 4.3, y =7.3, label = "c."), size = 4, fontface = "bold") +
  geom_text(aes(x = 7.3, y =7.3, label = "d."), size = 4, fontface = "bold") +
  geom_text(aes(x = 10.5, y = 19.7, label = "e."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 19.7, label = "f."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 14.5, label = "g."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 9.5, label = "h."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 4.5, label = "i."), size = 4, fontface = "bold")


map_year_beta_bydomain_annotate_theilsen 

ggsave(map_year_beta_bydomain_annotate_theilsen, path = file.path("Figures","Supplement","Theil_Sen"), filename = "map_year_beta_bydomain_annotate_theilsen.jpg", height =5, width = 13.5)

#biomassgradient component of BC dissimilarity
#outer, middle, inner
BC_dissim_O_M_I_biomassgradient <- ggplot(EBS.distances_dissimilarities_allyears[Domain != "Full",]) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_gradient_mean, color = Domain), size = 1) +
  labs(color = "Domain", x = "Year", y = "β diversity") +
  geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), color = "black", method = "lm", se = F, linetype = "longdash", linewidth = 0.6) +
  scale_color_manual(values = c("#999933", "#44AA99","#AA4499")) +
  facet_wrap(~Domain) +
  theme_classic() +
  theme(legend.position = "null", strip.text = element_blank())

#full
BC_dissim_F_bray_curtis_biomassgradient <- ggplot(EBS.distances_dissimilarities_allyears[Domain == "Full",]) +
  geom_point(aes(x = year, y = bray_curtis_dissimilarity_gradient_mean)) +
  labs(x = "Year", y = "β diversity\n(Biomass gradient component of Bray Curtis dissimilarity)") +
  geom_smooth(aes(x = year, y = bray_curtis_dissimilarity_gradient_mean), color = "black",
              method = "lm", se = F, linetype = "longdash") +
  theme_classic() +
  theme(legend.position = "null")

#place model values on plot
map_year_beta_bydomain_annotate_bray_curtis_biomassgradient <- ggdraw(xlim = c(0,30), ylim = c(0,20)) +
  draw_plot(BC_dissim_O_M_I_biomassgradient + theme(axis.title.x = element_blank()), x = 0, y = 0, width = 10, height = 8) +
  draw_plot(Alaska_domains, x = 0, y =6.5, width = 10, height = 15) +
  draw_plot(BC_dissim_F_bray_curtis_biomassgradient + theme(axis.title.x = element_blank()), x = 10, y = 0, width = 13, height = 20) +
  geom_text(aes(x = 14, y = 18, label = paste0("slope = ",slope_bray_curtis_full_biomassgradient,"+/-",slope_se_bray_curtis_full_biomassgradient,"\np = ",p_value_bray_curtis_full_biomassgradient)), size = 3) +
  draw_plot(broken_stick_plot_w3_full_biomass_gradient_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 15,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w10_full_biomass_gradient_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 10,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w20_full_biomass_gradient_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 5,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w35_full_biomass_gradient_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 0,  width = 7,  height = 5) +
  geom_text(aes(x = 0.5, y = 19.7, label = "a."), size = 4, fontface = "bold") +
  geom_text(aes(x = 1.5, y = 7.3, label = "b."), size = 4, fontface = "bold") +
  geom_text(aes(x = 4.5, y =7.3, label = "c."), size = 4, fontface = "bold") +
  geom_text(aes(x = 7.5, y =7.3, label = "d."), size = 4, fontface = "bold") +
  geom_text(aes(x = 10.5, y = 19.7, label = "e."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 19.7, label = "f."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 14.5, label = "g."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 9.5, label = "h."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 4.5, label = "i."), size = 4, fontface = "bold")


map_year_beta_bydomain_annotate_bray_curtis_biomassgradient

ggsave(map_year_beta_bydomain_annotate_bray_curtis_biomassgradient, path = file.path("Figures","Supplement","Biomass_gradient"), filename = "map_year_beta_bydomain_annotate_bray_curtis_biomassgradient.jpg", height =5, width = 13.5)


#JACCARD
#outer, middle, inner
BC_dissim_O_M_I_jaccard <- ggplot(EBS.distances_dissimilarities_allyears[Domain != "Full",]) +
  geom_point(aes(x = year, y = jaccard_dissimilarity_turnover_mean, color = Domain), size = 1) +
  labs(color = "Domain", x = "Year", y = "β diversity (Jaccard)") +
  geom_smooth(aes(x = year, y = jaccard_dissimilarity_turnover_mean), color = "black", method = "lm", se = F, linetype = "longdash", linewidth = 0.6) +
  scale_color_manual(values = c("#999933", "#44AA99","#AA4499")) +
  facet_wrap(~Domain) +
  theme_classic() +
  theme(legend.position = "null", strip.text = element_blank())

#full
BC_dissim_F_jaccard <- ggplot(EBS.distances_dissimilarities_allyears[Domain == "Full",]) +
  geom_point(aes(x = year, y = jaccard_dissimilarity_turnover_mean)) +
  labs(x = "Year", y = "β diversity (Jaccard)") +
  geom_smooth(aes(x = year, y = jaccard_dissimilarity_turnover_mean), color = "black",
              method = "lm", se = F, linetype = "longdash") +
  theme_classic() +
  theme(legend.position = "null")

#place model values on plot
map_year_beta_bydomain_annotate_jaccard <- ggdraw(xlim = c(0,30), ylim = c(0,20)) +
  draw_plot(BC_dissim_O_M_I_jaccard + theme(axis.title.x = element_blank()), x = 0, y = 0, width = 10, height = 8) +
  draw_plot(Alaska_domains, x = 0, y =6.5, width = 10, height = 15) +
  draw_plot(BC_dissim_F_jaccard + theme(axis.title.x = element_blank()), x = 10, y = 0, width = 13, height = 20) +
  geom_text(aes(x = 13, y = 18, label = paste0("slope = ",slope_jaccard_full,"+/-",slope_se_jaccard_full,"\np = ",p_value_jaccard_full)), size = 3) +
  draw_plot(broken_stick_plot_w3_full_jaccard_turnover_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 15,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w10_full_jaccard_turnover_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 10,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w20_full_jaccard_turnover_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 5,  width = 7,  height = 5) +
  draw_plot(broken_stick_plot_w35_full_jaccard_turnover_fig1 + theme(plot.title = element_text(size = 8), axis.title = element_blank()), x = 23, y = 0,  width = 7,  height = 5) +
  geom_text(aes(x = 0.5, y = 19.7, label = "a."), size = 4, fontface = "bold") +
  geom_text(aes(x = 1.3, y = 7.3, label = "b."), size = 4, fontface = "bold") +
  geom_text(aes(x = 4.3, y =7.3, label = "c."), size = 4, fontface = "bold") +
  geom_text(aes(x = 7.3, y =7.3, label = "d."), size = 4, fontface = "bold") +
  geom_text(aes(x = 10.5, y = 19.7, label = "e."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 19.7, label = "f."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 14.5, label = "g."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 9.5, label = "h."), size = 4, fontface = "bold") +
  geom_text(aes(x = 23, y = 4.5, label = "i."), size = 4, fontface = "bold")


map_year_beta_bydomain_annotate_jaccard

ggsave(map_year_beta_bydomain_annotate_jaccard, path = file.path("Figures","Supplement","Jaccard"), filename = "map_year_beta_bydomain_annotate_jaccard.jpg", height =5, width = 13.5)


#NO POLLOCK

#place model values on plot
map_year_beta_bydomain_annotate_nopollock <- ggdraw(xlim = c(0,10), ylim = c(0,20)) +
  draw_plot(EBS.domains.r.nopollock, x = 0, y = 0, width = 9, height = 10) +
  draw_plot(Alaska_domains, x = 0.2, y = 10, width = 10, height = 10) +
  geom_text(aes(x = 2.8, y = 1.5, label = paste0("Slope = ",slope_middle_nopollock, "  R^2 = ",R2_middle_nopollock, "  p = ",p_value_middle_nopollock)), size = 3) +
  geom_text(aes(x = 7, y = 1.5, label = paste0("Slope = ",slope_outer_nopollock, "  R^2 = ",R2_outer_nopollock, "  p = ",p_value_outer_nopollock)), size = 3) +
  geom_text(aes(x = 2.8, y = 5.8, label = paste0("Slope = ",slope_full_nopollock, "  R^2 = ",R2_full_nopollock, "  p = ",p_value_full_nopollock)), size = 3) +
  geom_text(aes(x = 7, y = 8.8, label = paste0("Slope = ",slope_inner_nopollock, "  R^2 = ",R2_inner_nopollock, "  p = ",p_value_inner_nopollock)), size = 3)

map_year_beta_bydomain_annotate_nopollock

save_plot(map_year_beta_bydomain_annotate_nopollock,
          filename = file.path("Figures","Supplement","Nopollock","map_year_beta_bydomain_annotate_nopollock.jpg"), base_height = 7, base_width = 6)

