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

#######################
##DATA##
#######################

#Pull grobs of map and trends (facet)

year_beta_bydomain <- readRDS(file.path("Figures","year_beta_bydomain.Rds")) #bray curtis balanced
year_beta_jaccard_bydomain <- readRDS(file.path("Figures","Supplement","Jaccard","year_beta_jaccard_bydomain.Rds")) #jaccard turnover
Alaska_domains <- readRDS(file.path("Figures","Alaska_domains.Rds"))

#pull in broken stick plots
broken_stick_plot_w3_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w3_full_lm_fig1.Rds"))
broken_stick_plot_w10_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w10_full_lm_fig1.Rds"))
broken_stick_plot_w20_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w20_full_lm_fig1.Rds"))
broken_stick_plot_w35_full_lm_fig1 <- readRDS(file.path("Figures","broken_stick_plot_w35_full_lm_fig1.Rds"))

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


#JACCARD (same)

#slopes
slope_jaccard_full <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Full",]))[[2]],5)
slope_jaccard_inner <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Inner",]))[[2]],5)
slope_jaccard_middle <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Middle",]))[[2]],5)
slope_jaccard_outer <- round(coefficients(lm(jaccard_dissimilarity_turnover_mean~year, data = EBS.distances_dissimilarities_allyears[domain == "Outer",]))[[2]],5)



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
#Bray Curtis with Pollock
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

#JACCARD
#place model values on plot
map_year_beta_jaccard_bydomain_annotate <- ggdraw(xlim = c(0,10), ylim = c(0,20)) +
  draw_plot(year_beta_jaccard_bydomain, x = 0, y = 0, width = 9, height = 10) +
  draw_plot(Alaska_domains, x = 0.2, y = 10, width = 10, height = 10) +
  geom_text(aes(x = 2.8, y = 1.5, label = paste0("Slope = ",slope_jaccard_middle, "  R^2 = ",R2_jaccard_middle, "  p = ",p_value_jaccard_middle)), size = 3) +
  geom_text(aes(x = 7, y = 1.5, label = paste0("Slope = ",slope_jaccard_outer, "  R^2 = ",R2_jaccard_outer, "  p = ",p_value_jaccard_outer)), size = 3) +
  geom_text(aes(x = 2.8, y = 5.8, label = paste0("Slope = ",slope_jaccard_full, "  R^2 = ",R2_jaccard_full, "  p = ",p_value_jaccard_full)), size = 3) +
  geom_text(aes(x = 7, y = 8.8, label = paste0("Slope = ",slope_jaccard_inner, "  R^2 = ",R2_jaccard_inner, "  p = ",p_value_jaccard_inner)), size = 3)

map_year_beta_jaccard_bydomain_annotate 

save_plot(map_year_beta_jaccard_bydomain_annotate, filename = file.path("Figures","Supplement","Jaccard","map_year_beta_jaccard_bydomain_annotate.jpg"), base_height = 7, base_width = 6)


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

