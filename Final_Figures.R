#All final figures
setwd()

library('ctv')
library('phytools')
library(ape)
library(phangorn)
library(devtools)
library(plyr)
library(dplyr)
library(devtools)
library(shiny)
devtools::github_release()
library(githubinstall)
devtools::install_github("vqv/ggbiplot")
library("ggbiplot") 
library(ggbiplot)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(tibble)
library(magrittr)
library(stringr)
library('gginnards')
library('ggpubr')
library('viridis')
library(dplyr)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)

#Figure 1 ----------------------------------------------------------------------
#Birth-Death Parameter Values -------------------------------------------------------------------------------------------------------------------------------
#Scenario 1 --------
load(file="tree_parameters.RData")
tree_bd_params <- data.frame(goodtreeparameters$birthparameters,
                             goodtreeparameters$deathparameters,
                             goodtreeparameters$rho)
colnames(tree_bd_params) <- c("birth", "death", "rho")
hist_birth_params <- ggplot(tree_bd_params, aes(x=log(birth))) + geom_histogram(aes(y=..count..), binwidth=0.05, color="#33638D", fill="#33638D") +
                        xlim(-10, 10) + geom_vline(aes(xintercept=log(median(birth))), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
                        labs(y = "Count", x = expression("log("*lambda*")")) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
hist_death_params <- ggplot(tree_bd_params, aes(x=log(death))) + geom_histogram(aes(y=..count..), binwidth=0.05, color="#33638D", fill="#33638D") + xlim(-10, 10) +
  geom_vline(aes(xintercept=log(median(death))), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = "log(µ)") + ylim(0,20) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
hist_rho_params <- ggplot(tree_bd_params, aes(x=rho)) + geom_histogram(aes(y=..count..), binwidth=0.01, color="#33638D", fill="#33638D") + xlim(0, 1) +
  geom_vline(aes(xintercept=median(rho)), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = expression(rho)) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
combined_bdparameters <- ggarrange(hist_birth_params, hist_death_params, hist_rho_params,
                                     ncol = 1, nrow = 3)


#Scenario 2------------------------------------------
load(file="finaltreeparams.RData")
load(file="tree_parameters.RData")
nooutgroup <- filter(finaltreeparams, longestpendant != age, EmpOrSim == "emp")
treenums <- c(nooutgroup$treenumber)
df <- data.frame(goodtreeparameters)
nooutgrouptree_params <- df[treenums,]
nooutgrouptree_params
hist_birth_params <- ggplot(nooutgrouptree_params, aes(x=log(birthparameters))) + geom_histogram(aes(y=..count..), binwidth=0.05, color="#33638D", fill="#33638D") +
  xlim(-10, 10) + geom_vline(aes(xintercept=log(median(birthparameters))), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = expression("log("*lambda*")")) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
hist_death_params <- ggplot(nooutgrouptree_params, aes(x=log(deathparameters))) + geom_histogram(aes(y=..count..), binwidth=0.05, color="#33638D", fill="#33638D") + xlim(-10, 10) +
  geom_vline(aes(xintercept=log(median(deathparameters))), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = "log(µ)") + ylim(0,20) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
hist_rho_params <- ggplot(nooutgrouptree_params, aes(x=rho)) + geom_histogram(aes(y=..count..), binwidth=0.01, color="#33638D", fill="#33638D") + xlim(0, 1) +
  geom_vline(aes(xintercept=median(rho)), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = expression(rho)) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
combined_bdparameters_nooutgroup <- ggarrange(hist_birth_params, hist_death_params, hist_rho_params,
                                   ncol = 1, nrow = 3)

#Scenario 3-------------------------------------------
load(file="tree_parameters_rho1.RData")
tree_bd_params_rho1 <- data.frame(goodtreeparameters_rho1$birthparameters,
                                  goodtreeparameters_rho1$deathparameters,
                                  goodtreeparameters_rho1$rho)
colnames(tree_bd_params_rho1) <- c("birth", "death", "rho")
hist_birth_params <- ggplot(tree_bd_params_rho1, aes(x=log(birth))) + geom_histogram(aes(y=..count..), binwidth=0.05, color="#33638D", fill="#33638D") +
  xlim(-10, 10) + geom_vline(aes(xintercept=log(median(birth))), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = expression("log("*lambda*")")) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
hist_death_params <- ggplot(tree_bd_params_rho1, aes(x=log(death))) + geom_histogram(aes(y=..count..), binwidth=0.05, color="#33638D", fill="#33638D") + xlim(-10, 10) +
  geom_vline(aes(xintercept=log(median(death))), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = "log(µ)") + ylim(0,20) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
combined_bdparameters_rho1 <- ggarrange(hist_birth_params, hist_death_params,
                                              ncol = 1, nrow = 2)

#Scenario 4-----------------------------------------
load(file="tree_parameters_yule.RData")
tree_bd_params_yule <- data.frame(goodtreeparameters_yule$birthparameters,
                                  goodtreeparameters_yule$rho)
colnames(tree_bd_params_yule) <- c("birth", "rho")
hist_birth_params <- ggplot(tree_bd_params_yule, aes(x=log(birth))) + geom_histogram(aes(y=..count..), binwidth=0.05, color="#33638D", fill="#33638D") +
  xlim(-10, 10) + geom_vline(aes(xintercept=log(median(birth))), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = expression("log("*lambda*")")) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
hist_rho_params <- ggplot(tree_bd_params_yule, aes(x=rho)) + geom_histogram(aes(y=..count..), binwidth=0.01, color="#33638D", fill="#33638D") + xlim(0, 1) +
  geom_vline(aes(xintercept=median(rho)), color="#20A486", linetype="dashed", size=1)+ theme_bw() +
  labs(y = "Count", x = expression(rho)) + ylim(0,75) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8))
combined_bdparameters_yule <- ggarrange(hist_birth_params, hist_rho_params,
                                              ncol = 1, nrow = 2)

#Combining into one plot -------------
BDparameter_plot <- ggarrange(combined_bdparameters, combined_bdparameters_nooutgroup, combined_bdparameters_rho1,
                              combined_bdparameters_yule + rremove("x.text"), 
                           labels = c("a)", "b)", "c)", "d)"),
                           heights = c(1.2,0.8),
                           label.y = 0.99,
                           label.x = -0.02,
                           ncol = 2, nrow = 2)
BDparameter_plot
ggsave(file="FinalBDparameterPlot.pdf", width = 6.75, height = 9)


#Finding the percentage of death/extinction parameter values that were zero
#Scen 1
table(tree_bd_params$death) #47 with zero = 3.95%
#Scen 2
table(nooutgrouptree_params$deathparameters) #22 with zero, 3.12% 
#Scen 3
table(tree_bd_params_rho1$death) #778 with zero, 65.43%




















#Figure 2 ----------------------------------------------------------------------------------------------------------------------------------------------
#Grouping each index by type: balance, imbalance, other with branch length, other without branch length-----------------------------------------------------------------
imbalance_list <- c("Total I","Mean I", "Median I", "Norm. Colless", "Colless-like index", 
                       "I2", "Norm. Sackin", "Total coph. index", "Avg. leaf depth",
                       "Leaf depth variance", "J index", "SNI", "S-shape statistic",
                       "Max. depth", "Mean depth", "stairs")

balance_list <- c("B1", "B2", "Max. width", "Max. diff. in widths", "Rooted Quartet", "Furnas rank", "stairs2")

other_no_branch_list <- c("Norm. cherries", "APP", "Norm. IL number", "Norm. avg. ladder length", 
                                   "Max. height", "Norm. pitchforks")
other_with_branch_list <- c("Longest pendant edge", "Shortest pendant edge", "Stemminess")

#Scenario 1 ----------------------------------------------
load(file="wilcoxon_distance_values.RData")
wilcoxon_distance_values <- wilcoxon_distance_values %>%
  mutate(Type = ifelse(name %in% imbalance_list, "Imbalance index",
                       ifelse(name %in% balance_list, "Balance index",
                              ifelse(name %in% other_no_branch_list, "Other index, no branch length",
                                     ifelse(name %in% other_with_branch_list, "Other index, with branch length", NA)))))
#Expected U value is n1*n2/2 = 1189*1189000/2 = 706860500
expectedU <- 706860500
#Main plot
plot1_wil <- ggplot(wilcoxon_distance_values, aes(x= reorder(name, abs(distance_value-expectedU)), y=distance_value-expectedU, fill=Type, color=)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  coord_flip() + scale_fill_manual(values=c('#440154', '#33638D', '#20A486', '#AADC32')) + ylim(-5*(10^8), 5*(10^8)) +
  theme(axis.title.x=element_text(size = 8), legend.position='none', axis.title.y=element_blank()) + labs(y="Difference between U Statistic and Expected U Value")
plot1_wil
#Dot plot for p-values
plot2_wil <- ggplot(wilcoxon_distance_values, aes(x=reorder(name, abs(distance_value-expectedU)), y=log(p_value))) + theme_classic() +
  geom_point(stat='identity', size=1.5, fill="#440154", colour="#440154") + coord_flip() + ylim(-350, 1) + 
  ylab("Corresponding log (p-value)") + geom_hline(yintercept = -2, linetype="dashed", color = "#440154", size=0.5) +
  theme(axis.title.x=element_text(size = 8), axis.title.y=element_blank(), legend.position='none', axis.text.y=element_blank())
plot2_wil
#combining into one figure
combined_wil <- ggarrange(plot1_wil, plot2_wil, widths = c(2, 0.6),
                          ncol = 2, nrow = 1)

#Scenario 2 ---------------------------------------------
load(file="wilcoxon_distance_values_nooutgroup.RData")
wilcoxon_distance_values_nooutgroup <- wilcoxon_distance_values_nooutgroup %>%
  mutate(Type = ifelse(name %in% imbalance_list, "Imbalance index",
                       ifelse(name %in% balance_list, "Balance index",
                              ifelse(name %in% other_no_branch_list, "Other index, no branch length",
                                     ifelse(name %in% other_with_branch_list, "Other index, with branch length", NA)))))
#Expected U value is n1*n2/2 = 706*706000/2 = 249218000
expectedU <- 249218000
#Main plot
plot1_wil_nooutgroup <- ggplot(wilcoxon_distance_values_nooutgroup, aes(x= reorder(name, abs(distance_value-expectedU)), y=distance_value-expectedU, fill=Type, color=)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  coord_flip() + scale_fill_manual(values=c('#440154', '#33638D', '#20A486', '#AADC32')) + ylim(-5*(10^8), 5*(10^8)) +
  theme(axis.title.x=element_text(size = 8), legend.position='none', axis.title.y=element_blank()) + labs(y="Difference between U Statistic and Expected U Value")
plot1_wil_nooutgroup
#Dot plot for p-values
plot2_wil_nooutgroup <- ggplot(wilcoxon_distance_values_nooutgroup, aes(x=reorder(name, abs(distance_value-expectedU)), y=log(p_value))) + theme_classic() +
  geom_point(stat='identity', size=1.5, fill="#440154", colour="#440154") + coord_flip() + ylim(-350, 1) + 
  ylab("Corresponding log (p-value)") + geom_hline(yintercept = -2, linetype="dashed", color = "#440154", size=0.5) +
  theme(axis.title.x=element_text(size = 8), axis.title.y=element_blank(), legend.position='none', axis.text.y=element_blank())
plot2_wil_nooutgroup
#combining into one figure
combined_wil_nooutgroup <- ggarrange(plot1_wil_nooutgroup, plot2_wil_nooutgroup, widths = c(2, 0.6),
                                     ncol = 2, nrow = 1)

#Scenario 3 --------------------------------------------------
load(file="wilcoxon_distance_values_rho1.RData")
wilcoxon_distance_values_rho1 <- wilcoxon_distance_values_rho1 %>%
  mutate(Type = ifelse(name %in% imbalance_list, "Imbalance index",
                       ifelse(name %in% balance_list, "Balance index",
                              ifelse(name %in% other_no_branch_list, "Other index, no branch length",
                                     ifelse(name %in% other_with_branch_list, "Other index, with branch length", NA)))))
#Expected U value is n1*n2/2 = 1189*1189000/2 = 706860500
expectedU <- 706860500
#Main plot
plot1_wil_rho1 <- ggplot(wilcoxon_distance_values_rho1, aes(x= reorder(name, abs(distance_value-expectedU)), y=distance_value-expectedU, fill=Type, color=)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + ylim(-5*(10^8), 5*(10^8)) +
  coord_flip() + scale_fill_manual(values=c('#440154', '#33638D', '#20A486', '#AADC32')) +
  theme(axis.title.x=element_text(size = 8), legend.position='none', axis.title.y=element_blank()) + labs(y="Difference between U Statistic and Expected U Value")
plot1_wil_rho1
#Dot plot for p-values
plot2_wil_rho1 <- ggplot(wilcoxon_distance_values_rho1, aes(x=reorder(name, abs(distance_value-expectedU)), y=log(p_value))) + theme_classic() +
  geom_point(stat='identity', size=1.5, fill="#440154", colour="#440154") + coord_flip() + ylim(-350, 1) + 
  ylab("Corresponding log (p-value)") + geom_hline(yintercept = -2, linetype="dashed", color = "#440154", size=0.5) +
  theme(axis.title.x=element_text(size = 8), axis.title.y=element_blank(), legend.position='none', axis.text.y=element_blank())
plot2_wil_rho1
#combining into one figure
combined_wil_rho1 <- ggarrange(plot1_wil_rho1, plot2_wil_rho1, widths = c(2, 0.6),
          ncol = 2, nrow = 1)

#Scenario 4------------------------------------------------------------
load(file="wilcoxon_distance_values_yule.RData")
wilcoxon_distance_values_yule <- wilcoxon_distance_values_yule %>%
  mutate(Type = ifelse(name %in% imbalance_list, "Imbalance index",
                       ifelse(name %in% balance_list, "Balance index",
                              ifelse(name %in% other_no_branch_list, "Other index, no branch length",
                                     ifelse(name %in% other_with_branch_list, "Other index, with branch length", NA)))))
#Expected U value is n1*n2/2 = 1189*1189000/2 = 706860500
expectedU <- 706860500
#Main plot
plot1_wil_yule <- ggplot(wilcoxon_distance_values_yule, aes(x= reorder(name, abs(distance_value-expectedU)), y=distance_value-expectedU, fill=Type, color=)) +
  theme_bw() +
  geom_bar(stat = "identity") + ylim(-5*(10^8), 5*(10^8)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  coord_flip() + scale_fill_manual(values=c('#440154', '#33638D', '#20A486', '#AADC32')) +
  theme(axis.title.x=element_text(size = 8), legend.position='none', axis.title.y=element_blank()) + labs(y="Difference between U Statistic and Expected U Value")
plot1_wil_yule
#Dot plot for p-values
plot2_wil_yule <- ggplot(wilcoxon_distance_values_yule, aes(x=reorder(name, abs(distance_value-expectedU)), y=log(p_value))) + theme_classic() +
  geom_point(stat='identity', size=1.5, fill="#440154", colour="#440154") + coord_flip() + ylim(-350, 1) + 
  ylab("Corresponding log (p-value)") + geom_hline(yintercept = -2, linetype="dashed", color = "#440154", size=0.5) +
  theme(axis.title.x=element_text(size = 8), axis.title.y=element_blank(), legend.position='none', axis.text.y=element_blank())
plot2_wil_yule
#combining into one figure
combined_wil_yule <- ggarrange(plot1_wil_yule, plot2_wil_yule, widths = c(2, 0.6),
          ncol = 2, nrow = 1)

#Getting the legend
# Extract the legend. Returns a gtable
legend_plot <- plot1_wil_yule <- ggplot(wilcoxon_distance_values_yule, aes(x= reorder(name, abs(distance_value-expectedU)), y=distance_value-expectedU, fill=Type, color=)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  coord_flip() + scale_fill_manual(values=c('#440154', '#33638D', '#20A486', '#AADC32')) +
  theme(axis.title.x=element_text(size = 8), axis.title.y=element_blank()) + labs(y="Difference between U Statistic and Expected U Value") +
  theme(legend.direction ="horizontal")
leg <- get_legend(legend_plot)
# Convert to a ggplot and print
legend <- as_ggplot(leg)


#Combining into one figure-------------------------------------
wilcoxon_plot <- ggarrange(combined_wil, combined_wil_nooutgroup, combined_wil_rho1,
          combined_wil_yule, legend + rremove("x.text"), 
          labels = c("a)", "b)", "c)", "d)"),
          label.y = 0.90,
          label.x = 0.05,
          heights = c(1,1,0.1),
          ncol = 2, nrow = 3)
wilcoxon_plot 
ggsave(file="FinalWilcoxonPlot.pdf", width = 15, height = 11.25)







#Figure 3 -----------------------------------------------------------------------

#Principal component analysis ----------------------------------------------------------------------------------------
load(file="finaltreeparams.RData")

#Scenario 1
#Only includes indices that are normalized for all trees, i.e., where they don't depend on tree size --------------------------------------
pca_variables <- data.frame(finaltreeparams[,c("MeanI", "NormColless", "I2", "NormSackin", "normcherry", "normILnumber", "normaverageladder",
                                               "normpitchforks", "stairs", "stairs2", "stemminess")])
colnames(pca_variables) <- c("Mean I", "Norm. Colless", 
         "I2", "Norm. Sackin", "Norm. # cherries", "Norm. IL number", "Norm. average ladder length", 
         "Norm. # pitchforks",
         "stairs", "stairs2", "Stemminess (non-cumulative)")
tree.pca <- prcomp(pca_variables, center = TRUE, scale. = TRUE)

summary(tree.pca)
print(tree.pca)
plot(tree.pca)
require(ggbiplot)
g <- ggbiplot(tree.pca, ellipse = TRUE, groups = finaltreeparams$EmpOrSim, alpha = 0, varname.size = 2)
g1 <- g + theme_bw() + 
  scale_color_manual(name="Tree Type", labels=c("Empirical", "Simulated"), values=c('#440154', '#20A486')) + ylim(-2, 3.5) + xlim(-3.5, 3) +
  xlab("Standardized PC 1 (57.5% explained var.)") +
  ylab("Standardized PC 2 (24.0% explained var.)") +
  scale_fill_discrete(labels=c('Empirical', 'Simulated')) +
  theme(legend.direction ="horizontal", 
        legend.position = "bottom") +
        guides(colour = guide_legend(override.aes = list(size=10, linewidth=5))) +
        labs(title="PCA, Scenario 1") +
  theme(legend.position='none')

#Scenario 2
#Only includes parameters that are normalized for all trees, i.e., where they don't depend on tree size --------------------------------------
load(file="finaltreeparams.RData")
#filtering to get the right trees, fitting the criteria for scenario 2
nooutgroup <- filter(finaltreeparams, longestpendant != age, EmpOrSim == "emp")
treenums <- nooutgroup$treenumber
nooutgrouptrees <- finaltreeparams %>% filter(treenumber %in% treenums)
nooutgrouptrees
pca_variables_nooutgroup <- data.frame(nooutgrouptrees[,c("MeanI", "NormColless", "I2", "NormSackin", "normcherry", "normILnumber", "normaverageladder",
                                                          "normpitchforks", "stairs", "stairs2", "stemminess")])
colnames(pca_variables_nooutgroup) <- c("Mean I", "Norm. Colless", 
                                        "I2", "Norm. Sackin", "Norm. # cherries", "Norm. IL number", "Norm. average ladder length", 
                                        "Norm. # pitchforks",
                                        "stairs", "stairs2", "Stemminess (non-cumulative)")

tree.pca.nooutgroup <- prcomp(pca_variables_nooutgroup, center = TRUE, scale. = TRUE)
require(ggbiplot)
g <- ggbiplot(tree.pca.nooutgroup, ellipse = TRUE, groups = nooutgrouptrees$EmpOrSim, alpha = 0, varname.size=2)
g4 <- g + theme_bw() + 
  scale_color_manual(name="Tree Type", labels=c("Empirical", "Simulated"), values=c('#440154', '#20A486')) + ylim(-2, 3.5) + xlim(-3.5, 3) +
  xlab("Standardized PC 1 (53.6% explained var.)") +
  ylab("Standardized PC 2 (20.6% explained var.)") +
  scale_fill_discrete(labels=c('Empirical', 'Simulated')) +
  theme(legend.direction ="horizontal", 
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=10, linewidth=5))) +
  labs(title="PCA, Scenario 2") +
  theme(legend.position='none')

#Scenario 3
#Only includes parameters that are normalized for all trees, i.e., where they don't depend on tree size --------------------------------------
load(file="finaltreeparams_rho1.RData")
pca_variables_rho1 <- data.frame(finaltreeparams_rho1[,c("MeanI", "NormColless", "I2", "NormSackin", "normcherry", "normILnumber", "normaverageladder",
                                               "normpitchforks", "stairs", "stairs2", "stemminess")])

colnames(pca_variables_rho1) <- c("Mean I", "Norm. Colless", 
                                  "I2", "Norm. Sackin", "Norm. # cherries", "Norm. IL number", "Norm. average ladder length", 
                                  "Norm. # pitchforks",
                                  "stairs", "stairs2", "Stemminess (non-cumulative)")

#Running PCA
tree.pca.rho1 <- prcomp(pca_variables_rho1, center = TRUE, scale. = TRUE)

#Plotting
require(ggbiplot)
g <- ggbiplot(tree.pca.rho1, ellipse = TRUE, groups = finaltreeparams_rho1$EmpOrSim, alpha = 0, varname.size = 2)
g5 <- g + theme_bw() + 
  scale_color_manual(name="Tree Type", labels=c("Empirical", "Simulated"), values=c('#440154', '#20A486')) + ylim(-2, 3.5) + xlim(-3.5, 3) +
  xlab("Standardized PC 1 (56.4% explained var.)") +
  ylab("Standardized PC 2 (21.7% explained var.)") +
  scale_fill_discrete(labels=c('Empirical', 'Simulated')) +
  theme(legend.direction ="horizontal", 
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=10, linewidth=5))) +
  labs(title="PCA, Scenario 3") +
  theme(legend.position='none')

#Scenario 4
#Only includes parameters that are normalized for all trees, i.e., where they don't depend on tree size --------------------------------------
load(file="finaltreeparams_yule.RData")
pca_variables_yule <- data.frame(finaltreeparams_yule[,c("MeanI", "NormColless", "I2", "NormSackin", "normcherry", "normILnumber", "normaverageladder",
                                                         "normpitchforks", "stairs", "stairs2", "stemminess")])

colnames(pca_variables_yule) <- c("Mean I", "Norm. Colless", 
                                  "I2", "Norm. Sackin", "Norm. # cherries", "Norm. IL number", "Norm. average ladder length", 
                                  "Norm. # pitchforks",
                                  "stairs", "stairs2", "Stemminess (non-cumulative)")

#Running PCA
tree.pca.yule <- prcomp(pca_variables_yule, center = TRUE, scale. = TRUE)

#Plotting
require(ggbiplot)
g <- ggbiplot(tree.pca.yule, ellipse = TRUE, groups = finaltreeparams_yule$EmpOrSim, alpha = 0, varname.size = 2)
g6 <- g + theme_bw() + 
  scale_color_manual(name="Tree Type", labels=c("Empirical", "Simulated"), values=c('#440154', '#20A486')) + ylim(-2, 3.5) + xlim(-3.5, 3) +
  xlab("Standardized PC 1 (56.3% explained var.)") +
  ylab("Standardized PC 2 (20.3% explained var.)") +
  scale_fill_discrete(labels=c('Empirical', 'Simulated')) +
  theme(legend.direction ="horizontal", 
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=10, linewidth=5))) +
  labs(title="PCA, Scenario 4") +
  theme(legend.position='none')

#Getting the legend
# Extract the legend. Returns a gtable
g7 <- g + theme_bw() + 
  scale_color_manual(name="Tree Type", labels=c("Empirical", "Simulated"), values=c('#440154', '#20A486')) + ylim(-2, 3.5) + xlim(-3.5, 3) +
  xlab("Standardized PC 1 (56.4% explained var.)") +
  ylab("Standardized PC 2 (21.7% explained var.)") +
  scale_fill_discrete(labels=c('Empirical', 'Simulated')) +
  theme(legend.direction ="horizontal", 
        legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=10, linewidth=5))) +
  labs(title="PCA, Scenario 4")
leg <- get_legend(g7)
# Convert to a ggplot and print
legend <- as_ggplot(leg)

#Combining into one big figure for all PCA plots
g8 <- ggarrange(g1, g4, g5, g6, legend + rremove("x.text"), 
                labels = c("a)", "b)", "c)", "d)", ""),
                label.y = 0.90,
                label.x = 0.05,
                heights = c(1,1,0.1), ncol = 2, nrow = 3)
g8
ggsave(file="PCAwithoutPoints.pdf", g8,  width = 12, height = 9)














#Figure 4 ----------------------------------------------------------------------------------------------

#Z-Score Plots --------------------------------------------------------------------------------------------------
# Five from each scenario------------------------
load(file="zscoresdf.RData")
load(file="zscoresdf_nooutgroup.RData")
load(file="zscoresdf_rho1.RData")
load(file="zscoresdf_yule.RData")


#Scenario 1: Normalized Colless, Leaf depth variance, Stemminess, Mean I, B2-------------------
hist_MeanI <- ggplot(zscoresdf, aes(x=MeanI)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-5, 7.5) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Mean I Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Mean I")
hist_NormColless <- ggplot(zscoresdf, aes(x=NormColless)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Normalized Colless Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Normalized Colless")
hist_leafdepthvariance <- ggplot(zscoresdf, aes(x=leafdepthvariance)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Leaf depth variance") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Leaf Depth Variance")
hist_stemminess <- ggplot(zscoresdf, aes(x=stemminess)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 3) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Non-Cumulative Stemminess Index Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Non-Cumulative Stemminess")
hist_B2 <- ggplot(zscoresdf, aes(x=B2)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 2.5) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "B2 Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="B2")
#Combining into one big figure for Scenario 1
combined_Z_S1 <- ggarrange(hist_MeanI + rremove("x.title"), hist_NormColless + rremove("x.title"),
                           hist_leafdepthvariance + rremove("x.title"), hist_stemminess + rremove("x.title"), hist_B2 + rremove("x.title"),
                ncol = 5, nrow = 1)

#Scenario 2: No Outgroup
hist_MeanI <- ggplot(zscoresdf_nooutgroupdf, aes(x=MeanI)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-5, 7.5) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Mean I Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Mean I")
hist_NormColless <- ggplot(zscoresdf_nooutgroupdf, aes(x=NormColless)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Normalized Colless Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Normalized Colless")
hist_leafdepthvariance <- ggplot(zscoresdf_nooutgroupdf, aes(x=leafdepthvariance)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Leaf depth variance") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Leaf Depth Variance")
hist_stemminess <- ggplot(zscoresdf_nooutgroupdf, aes(x=stemminess)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 3) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Non-Cumulative Stemminess Index Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Non-Cumulative Stemminess")
hist_B2 <- ggplot(zscoresdf_nooutgroupdf, aes(x=B2)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 2.5) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "B2 Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="B2")
#Combining into one big figure for Scenario 2
combined_Z_S2 <- ggarrange(hist_MeanI + rremove("x.title"), hist_NormColless + rremove("x.title"),
                           hist_leafdepthvariance + rremove("x.title"), hist_stemminess + rremove("x.title"), hist_B2 + rremove("x.title"),
                           ncol = 5, nrow = 1)


#Scenario 3: Rho = 1
hist_MeanI <- ggplot(zscoresdf_rho1, aes(x=MeanI)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-5, 7.5) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Mean I Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Mean I")
hist_NormColless <- ggplot(zscoresdf_rho1, aes(x=NormColless)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Normalized Colless Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Normalized Colless")
hist_leafdepthvariance <- ggplot(zscoresdf_rho1, aes(x=leafdepthvariance)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Leaf depth variance") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Leaf Depth Variance")
hist_stemminess <- ggplot(zscoresdf_rho1, aes(x=stemminess)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 3) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Non-Cumulative Stemminess Index Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Non-Cumulative Stemminess")
hist_B2 <- ggplot(zscoresdf_rho1, aes(x=B2)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 2.5) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "B2 Value") + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="B2")
#Combining into one big figure for Scenario 3
combined_Z_S3 <- ggarrange(hist_MeanI + rremove("x.title"), hist_NormColless + rremove("x.title"),
                           hist_leafdepthvariance + rremove("x.title"), hist_stemminess + rremove("x.title"), hist_B2 + rremove("x.title"),
                           ncol = 5, nrow = 1)


#Scenario 4: Yule
hist_MeanI <- ggplot(zscoresdf_yule, aes(x=MeanI)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-5, 7.5) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Z-score") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Mean I")
hist_NormColless <- ggplot(zscoresdf_yule, aes(x=NormColless)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Z-score") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Normalized Colless")
hist_leafdepthvariance <- ggplot(zscoresdf_yule, aes(x=leafdepthvariance)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-2.5, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Z-score") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Leaf Depth Variance")
hist_stemminess <- ggplot(zscoresdf_yule, aes(x=stemminess)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 3) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Z-score") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="Non-Cumulative Stemminess")
hist_B2 <- ggplot(zscoresdf_yule, aes(x=B2)) + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 2.5) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)+
  theme_bw() + labs(y = "Count", x = "Z-score") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 8), plot.title = element_text(size=10)) +
  labs(title="B2")
#Combining into one big figure for Scenario 4
combined_Z_S4 <- ggarrange(hist_MeanI, hist_NormColless,
                           hist_leafdepthvariance, hist_stemminess, hist_B2,
                           ncol = 5, nrow = 1)


#Combining into one figure for all scenarios ----------------------
combined_Z_ALL <- ggarrange(combined_Z_S1, combined_Z_S2 + rremove("x.title"), combined_Z_S3 + rremove("x.title"),
                            combined_Z_S4 + rremove("x.title"), labels = c("a)", "b)", "c)", "d)"),
                            label.y = 0.95,
                            label.x = 0.01,
                            heights = c(1,1,1,1.1),
                ncol = 1, nrow = 4)
combined_Z_ALL
ggsave(file="Z_Score_Distribution_5_Indices.pdf", combined_Z_ALL,  width = 12, height = 9)







#Figure 5 ----------------------------------------------------------------------------------------------------------------------------------
#Average Z-score per study - trying to find whether specific studies are outliers ---------------------------
load(file="zscoresdf.RData")
load(file="zscoresdf_nooutgroup.RData")
load(file="zscoresdf_rho1.RData")
load(file="zscoresdf_yule.RData")

#Scenario 1 --------------------------------------------------------------------------------------------------------------------
d <- c(3:8, 10:15,17:36) #Removing the CP Rank index, since it was no computationally tractable
e <- zscoresdf[,d]
scen1_study_z <- data.frame((rowMeans(abs(e))))
colnames(scen1_study_z) <- c("Z_score_mean")

hist_studyZscore <- ggplot(scen1_study_z, aes(x=Z_score_mean)) + geom_histogram(aes(y=..count..),binwidth=0.05, color="#33638D", fill="#33638D") + xlim(0, 10.5) + ylim(0, 55) +
  theme_bw() + labs(y = "Count", x = "Mean of Absolute Z-Score Values") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 10), plot.title = element_text(size=10))
hist_studyZscore

#Scenario 2 --------------------------------------------------------------------------------------------------------------------
d <- c(4:9, 11:16,18:36) #Removing the CP Rank index, since it was no computationally tractable
e <- zscoresdf_nooutgroupdf[,d]
scen1_study_z_nooutgroup <- data.frame((rowMeans(abs(e))))
colnames(scen1_study_z_nooutgroup) <- c("Z_score_mean")

hist_studyZscore_nooutgroup <- ggplot(scen1_study_z_nooutgroup, aes(x=Z_score_mean)) + geom_histogram(aes(y=..count..),binwidth=0.05, color="#33638D", fill="#33638D") + xlim(0, 10.5) + ylim(0, 32) +
  theme_bw() + labs(y = "Count", x = "Mean of Absolute Z-Score Values") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 10), plot.title = element_text(size=10))
hist_studyZscore_nooutgroup

#Scenario 3 --------------------------------------------------------------------------------------------------------------------
d <- c(2:7, 9:14,16:35) #Removing the CP Rank index, since it was no computationally tractable
e <- zscoresdf_rho1[,d]
scen1_study_z_rho1 <- data.frame(rowMeans(abs(e)))
colnames(scen1_study_z_rho1) <- c("Z_score_mean")

hist_studyZscore_rho1 <- ggplot(scen1_study_z_rho1, aes(x=Z_score_mean)) + geom_histogram(aes(y=..count..),binwidth=0.05, color="#33638D", fill="#33638D") + xlim(0, 10.5) + ylim(0, 40) +
  theme_bw() + labs(y = "Count", x = "Mean of Absolute Z-Score Values") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 10), plot.title = element_text(size=10))
hist_studyZscore_rho1

#Scenario 4 --------------------------------------------------------------------------------------------------------------------
d <- c(2:7, 9:14,16:35) #Removing the CP Rank index, since it was no computationally tractable
e <- zscoresdf_yule[,d]
scen1_study_z_yule <- data.frame((rowMeans(abs(e))))
colnames(scen1_study_z_yule) <- c("Z_score_mean")

hist_studyZscore_yule <- ggplot(scen1_study_z_yule, aes(x=Z_score_mean)) + geom_histogram(aes(y=..count..),binwidth=0.05, color="#33638D", fill="#33638D") + xlim(0, 10.5) + ylim(0, 50) +
  theme_bw() + labs(y = "Count", x = "Mean of absolute z-score values") + theme(axis.title.x=element_text(size = 10), axis.title.y=element_text(size = 10), plot.title = element_text(size=10))
hist_studyZscore_yule

#combining into one figure----------------------------------------------------------
combined_hist_Z_study <- ggarrange(hist_studyZscore + rremove("x.title"), hist_studyZscore_nooutgroup + rremove("x.title"), hist_studyZscore_rho1, hist_studyZscore_yule,
                            labels = c("a)", "b)", "c)", "d)"),
                            label.y = 0.95,
                            label.x = 0.01,
                            heights = c(1,1,1,1.1),
                            ncol = 2, nrow = 2)
combined_hist_Z_study
ggsave(file="Mean_Absolute_Z_Score_Distribution_ByStudy.pdf", combined_hist_Z_study, width = 9, height = 9)


#Finding which studies were outliers (mean z-score >5)
outlier <- which(scen1_study_z > 5) #Scenario 1
which(scen1_study_z_yule > 5) #Scenario 4, same list as for Scenario 1
load(file="final_timetrees.RData")
outliertrees <- trfinallist[outlier]
finaltreeparams$leaves[outlier]
#Finding tree size of outliers
mean(finaltreeparams$leaves[outlier])





#SUPPLEMENTARY FIGURES -----------------------------------------------------------------------------

#Supplementary 4-7 - All Z-scores for all indices for each scenario

#Scenario 1 --------------------
load(file="zscoresdf.RData")
library('ggpubr')
library(gridExtra)
#Histograms of Z-scores -----------------------------------------------------------------------------------------
hist_TotalI <- ggplot(zscoresdf, aes(x=TotalI)) + labs(y = "Count", x = "Z-score", title="Total I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalI)), color="#20A486", linetype="dashed", size=1)
hist_MeanI <- ggplot(zscoresdf, aes(x=MeanI)) + labs(y = "Count", x = "Z-score", title="Mean I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)
hist_MedianI <- ggplot(zscoresdf, aes(x=MedianI)) + labs(y = "Count", x = "Z-score", title="Median I") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MedianI)), color="#20A486", linetype="dashed", size=1)
hist_NormColless <- ggplot(zscoresdf, aes(x=NormColless)) + labs(y = "Count", x = "Z-score", title="Norm. Colless") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)
hist_CollessLikeExpMDM <- ggplot(zscoresdf, aes(x=CollessLikeExpMDM)) +  labs(y = "Count", x = "Z-score", title="Colless Like Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(CollessLikeExpMDM)), color="#20A486", linetype="dashed", size=1)
hist_I2 <- ggplot(zscoresdf, aes(x=I2)) + labs(y = "Count", x = "Z-score", title="I2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(I2)), color="#20A486", linetype="dashed", size=1)
hist_NormSackin <- ggplot(zscoresdf, aes(x=NormSackin)) + labs(y = "Count", x = "Z-score", title="Norm. Sackin") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormSackin)), color="#20A486", linetype="dashed", size=1)
hist_TotalCophenetic <- ggplot(zscoresdf, aes(x=TotalCophenetic)) + labs(y = "Count", x = "Z-score", title="Total Cophenetic Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalCophenetic)), color="#20A486", linetype="dashed", size=1)
hist_B1 <- ggplot(zscoresdf, aes(x=B1)) + labs(y = "Count", x = "Z-score", title="B1") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B1)), color="#20A486", linetype="dashed", size=1)
hist_B2 <- ggplot(zscoresdf, aes(x=B2)) + labs(y = "Count", x = "Z-score", title="B2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)
hist_AvgLeafDepth <- ggplot(zscoresdf, aes(x=AvgLeafDepth)) + labs(y = "Count", x = "Z-score", title="Average Leaf Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(AvgLeafDepth)), color="#20A486", linetype="dashed", size=1)
hist_leafdepthvariance <- ggplot(zscoresdf, aes(x=leafdepthvariance)) + labs(y = "Count", x = "Z-score", title="Leaf Depth Variance") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)
hist_normcherry <- ggplot(zscoresdf, aes(x=normcherry)) + labs(y = "Count", x = "Z-score", title="Norm. Cherry Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normcherry)), color="#20A486", linetype="dashed", size=1)
hist_Jindex <- ggplot(zscoresdf, aes(x=Jindex)) + labs(y = "Count", x = "Z-score", title="J Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Jindex)), color="#20A486", linetype="dashed", size=1)
hist_SNI <- ggplot(zscoresdf, aes(x=SNI)) +  labs(y = "Count", x = "Z-score", title="SNI") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(SNI)), color="#20A486", linetype="dashed", size=1)
hist_Sshape <- ggplot(zscoresdf, aes(x=Sshape)) +  labs(y = "Count", x = "Z-score", title="S-shape Statistic") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Sshape)), color="#20A486", linetype="dashed", size=1)
hist_APP <- ggplot(zscoresdf, aes(x=APP)) + labs(y = "Count", x = "Z-score", title="APP") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(APP)), color="#20A486", linetype="dashed", size=1)
hist_normILnumber <- ggplot(zscoresdf, aes(x=normILnumber)) + labs(y = "Count", x = "Z-score", title="Norm. IL number") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normILnumber)), color="#20A486", linetype="dashed", size=1)
hist_normaverageladder <- ggplot(zscoresdf, aes(x=normaverageladder)) + labs(y = "Count", x = "Z-score", title="Norm. Avg. Ladder Length") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normaverageladder)), color="#20A486", linetype="dashed", size=1)
hist_maxdepth <- ggplot(zscoresdf, aes(x=maxdepth)) + labs(y = "Count", x = "Z-score", title="Max Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdepth)), color="#20A486", linetype="dashed", size=1)
hist_maxwidth <- ggplot(zscoresdf, aes(x=maxwidth)) + labs(y = "Count", x = "Z-score", title="Max Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxwidth)), color="#20A486", linetype="dashed", size=1)
hist_maxdiffwidth <- ggplot(zscoresdf, aes(x=maxdiffwidth)) + labs(y = "Count", x = "Z-score", title="Max. Diff. in Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdiffwidth)), color="#20A486", linetype="dashed", size=1)
hist_meandepth <- ggplot(zscoresdf, aes(x=meandepth)) + labs(y = "Count", x = "Z-score", title="Mean Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(meandepth)), color="#20A486", linetype="dashed", size=1)
hist_maxheight <- ggplot(zscoresdf, aes(x=maxheight)) + labs(y = "Count", x = "Z-score", title="Max Height") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxheight)), color="#20A486", linetype="dashed", size=1)
hist_normpitchforks <- ggplot(zscoresdf, aes(x=normpitchforks)) + labs(y = "Count", x = "Z-score", title="Norm. # of Pitchforks") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normpitchforks)), color="#20A486", linetype="dashed", size=1)
hist_rootedquartet <- ggplot(zscoresdf, aes(x=rootedquartet)) + labs(y = "Count", x = "Z-score", title="Rooted Quartet Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(rootedquartet)), color="#20A486", linetype="dashed", size=1)
hist_furnas <- ggplot(zscoresdf, aes(x=furnas)) + labs(y = "Count", x = "Z-score", title="Furnas Rank") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(furnas)), color="#20A486", linetype="dashed", size=1)
hist_stairs <- ggplot(zscoresdf, aes(x=stairs)) + labs(y = "Count", x = "Z-score", title="stairs") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs)), color="#20A486", linetype="dashed", size=1)
hist_stairs2 <- ggplot(zscoresdf, aes(x=stairs2)) + labs(y = "Count", x = "Z-score", title="stairs2") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs2)), color="#20A486", linetype="dashed", size=1)
hist_longestpendant <- ggplot(zscoresdf, aes(x=longestpendant)) + labs(y = "Count", x = "Z-score", title="Longest Pendant Edge") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(longestpendant)), color="#20A486", linetype="dashed", size=1)
hist_shortestpendant <- ggplot(zscoresdf, aes(x=shortestpendant)) + labs(y = "Count", x = "Z-score", title="Shortest Pendant Edge") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(shortestpendant)), color="#20A486", linetype="dashed", size=1)
hist_stemminess <- ggplot(zscoresdf, aes(x=stemminess)) + labs(y = "Count", x = "Z-score", title="Non-Cumulative Stemminess") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)


zhist_all <- list(hist_TotalI, hist_MeanI, hist_MedianI, hist_NormColless, hist_CollessLikeExpMDM, hist_I2, hist_NormSackin, 
                  hist_TotalCophenetic, hist_B1, hist_B2, hist_AvgLeafDepth, hist_leafdepthvariance, hist_normcherry, hist_Jindex, hist_SNI, 
                  hist_Sshape, hist_APP, hist_normILnumber, hist_normaverageladder,  hist_maxdepth, 
                  hist_maxwidth, hist_maxdiffwidth, hist_meandepth, hist_maxheight,
                  hist_normpitchforks, hist_rootedquartet, hist_stairs, hist_stairs2, hist_longestpendant,
                  hist_shortestpendant, hist_stemminess)
save(zhist_all, file="zscore_histograms.RData")
load(file="zscore_histograms.RData")
ggsave(file="AllZscorePlots_Scenario1.pdf", marrangeGrob(grobs = zhist_all, nrow=4, ncol=3, top = "Z-scores for Scenario 1 Indices; Dotted Lines Denote the Sample Mean"), width = 9, height = 12,
       device = "pdf")

#Scenario 2 --------------------
load(file="zscoresdf_nooutgroup.RData")
#Histograms of Z-scores -----------------------------------------------------------------------------------------
hist_TotalI <- ggplot(zscoresdf_nooutgroupdf, aes(x=TotalI)) + labs(y = "Count", x = "Z-score", title="Total I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalI)), color="#20A486", linetype="dashed", size=1)
hist_MeanI <- ggplot(zscoresdf_nooutgroupdf, aes(x=MeanI)) + labs(y = "Count", x = "Z-score", title="Mean I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)
hist_MedianI <- ggplot(zscoresdf_nooutgroupdf, aes(x=MedianI)) + labs(y = "Count", x = "Z-score", title="Median I") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MedianI)), color="#20A486", linetype="dashed", size=1)
hist_NormColless <- ggplot(zscoresdf_nooutgroupdf, aes(x=NormColless)) + labs(y = "Count", x = "Z-score", title="Norm. Colless") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)
hist_CollessLikeExpMDM <- ggplot(zscoresdf_nooutgroupdf, aes(x=CollessLikeExpMDM)) +  labs(y = "Count", x = "Z-score", title="Colless Like Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(CollessLikeExpMDM)), color="#20A486", linetype="dashed", size=1)
hist_I2 <- ggplot(zscoresdf_nooutgroupdf, aes(x=I2)) + labs(y = "Count", x = "Z-score", title="I2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(I2)), color="#20A486", linetype="dashed", size=1)
hist_NormSackin <- ggplot(zscoresdf_nooutgroupdf, aes(x=NormSackin)) + labs(y = "Count", x = "Z-score", title="Norm. Sackin") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormSackin)), color="#20A486", linetype="dashed", size=1)
hist_TotalCophenetic <- ggplot(zscoresdf_nooutgroupdf, aes(x=TotalCophenetic)) + labs(y = "Count", x = "Z-score", title="Total Cophenetic Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalCophenetic)), color="#20A486", linetype="dashed", size=1)
hist_B1 <- ggplot(zscoresdf_nooutgroupdf, aes(x=B1)) + labs(y = "Count", x = "Z-score", title="B1") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B1)), color="#20A486", linetype="dashed", size=1)
hist_B2 <- ggplot(zscoresdf_nooutgroupdf, aes(x=B2)) + labs(y = "Count", x = "Z-score", title="B2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)
hist_AvgLeafDepth <- ggplot(zscoresdf_nooutgroupdf, aes(x=AvgLeafDepth)) + labs(y = "Count", x = "Z-score", title="Average Leaf Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(AvgLeafDepth)), color="#20A486", linetype="dashed", size=1)
hist_leafdepthvariance <- ggplot(zscoresdf_nooutgroupdf, aes(x=leafdepthvariance)) + labs(y = "Count", x = "Z-score", title="Leaf Depth Variance") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)
hist_normcherry <- ggplot(zscoresdf_nooutgroupdf, aes(x=normcherry)) + labs(y = "Count", x = "Z-score", title="Norm. Cherry Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normcherry)), color="#20A486", linetype="dashed", size=1)
hist_Jindex <- ggplot(zscoresdf_nooutgroupdf, aes(x=Jindex)) + labs(y = "Count", x = "Z-score", title="J Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Jindex)), color="#20A486", linetype="dashed", size=1)
hist_SNI <- ggplot(zscoresdf_nooutgroupdf, aes(x=SNI)) +  labs(y = "Count", x = "Z-score", title="SNI") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(SNI)), color="#20A486", linetype="dashed", size=1)
hist_Sshape <- ggplot(zscoresdf_nooutgroupdf, aes(x=Sshape)) +  labs(y = "Count", x = "Z-score", title="S-shape Statistic") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Sshape)), color="#20A486", linetype="dashed", size=1)
hist_APP <- ggplot(zscoresdf_nooutgroupdf, aes(x=APP)) + labs(y = "Count", x = "Z-score", title="APP") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(APP)), color="#20A486", linetype="dashed", size=1)
hist_normILnumber <- ggplot(zscoresdf_nooutgroupdf, aes(x=normILnumber)) + labs(y = "Count", x = "Z-score", title="Norm. IL number") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normILnumber)), color="#20A486", linetype="dashed", size=1)
hist_normaverageladder <- ggplot(zscoresdf_nooutgroupdf, aes(x=normaverageladder)) + labs(y = "Count", x = "Z-score", title="Norm. Avg. Ladder Length") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normaverageladder)), color="#20A486", linetype="dashed", size=1)
hist_maxdepth <- ggplot(zscoresdf_nooutgroupdf, aes(x=maxdepth)) + labs(y = "Count", x = "Z-score", title="Max Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdepth)), color="#20A486", linetype="dashed", size=1)
hist_maxwidth <- ggplot(zscoresdf_nooutgroupdf, aes(x=maxwidth)) + labs(y = "Count", x = "Z-score", title="Max Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxwidth)), color="#20A486", linetype="dashed", size=1)
hist_maxdiffwidth <- ggplot(zscoresdf_nooutgroupdf, aes(x=maxdiffwidth)) + labs(y = "Count", x = "Z-score", title="Max. Diff. in Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdiffwidth)), color="#20A486", linetype="dashed", size=1)
hist_meandepth <- ggplot(zscoresdf_nooutgroupdf, aes(x=meandepth)) + labs(y = "Count", x = "Z-score", title="Mean Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(meandepth)), color="#20A486", linetype="dashed", size=1)
hist_maxheight <- ggplot(zscoresdf_nooutgroupdf, aes(x=maxheight)) + labs(y = "Count", x = "Z-score", title="Max Height") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxheight)), color="#20A486", linetype="dashed", size=1)
hist_normpitchforks <- ggplot(zscoresdf_nooutgroupdf, aes(x=normpitchforks)) + labs(y = "Count", x = "Z-score", title="Norm. # of Pitchforks") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normpitchforks)), color="#20A486", linetype="dashed", size=1)
hist_rootedquartet <- ggplot(zscoresdf_nooutgroupdf, aes(x=rootedquartet)) + labs(y = "Count", x = "Z-score", title="Rooted Quartet Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(rootedquartet)), color="#20A486", linetype="dashed", size=1)
hist_furnas <- ggplot(zscoresdf_nooutgroupdf, aes(x=furnas)) + labs(y = "Count", x = "Z-score", title="Furnas Rank") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(furnas)), color="#20A486", linetype="dashed", size=1)
hist_stairs <- ggplot(zscoresdf_nooutgroupdf, aes(x=stairs)) + labs(y = "Count", x = "Z-score", title="stairs") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs)), color="#20A486", linetype="dashed", size=1)
hist_stairs2 <- ggplot(zscoresdf_nooutgroupdf, aes(x=stairs2)) + labs(y = "Count", x = "Z-score", title="stairs2") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs2)), color="#20A486", linetype="dashed", size=1)
hist_shortestpendant <- ggplot(zscoresdf_nooutgroupdf, aes(x=shortestpendant)) + labs(y = "Count", x = "Z-score", title="Shortest Pendant Edge") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(shortestpendant)), color="#20A486", linetype="dashed", size=1)
hist_stemminess <- ggplot(zscoresdf_nooutgroupdf, aes(x=stemminess)) + labs(y = "Count", x = "Z-score", title="Non-Cumulative Stemminess") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)


zhist_all <- list(hist_TotalI, hist_MeanI, hist_MedianI, hist_NormColless, hist_CollessLikeExpMDM, hist_I2, hist_NormSackin, 
                  hist_TotalCophenetic, hist_B1, hist_B2, hist_AvgLeafDepth, hist_leafdepthvariance, hist_normcherry, hist_Jindex, hist_SNI, 
                  hist_Sshape, hist_APP, hist_normILnumber, hist_normaverageladder,  hist_maxdepth, 
                  hist_maxwidth, hist_maxdiffwidth, hist_meandepth, hist_maxheight,
                  hist_normpitchforks, hist_rootedquartet, hist_stairs, hist_stairs2,
                  hist_shortestpendant, hist_stemminess)
ggsave(file="AllZscorePlots_Scenario2.pdf", marrangeGrob(grobs = zhist_all, nrow=4, ncol=3, top = "Z-scores for Scenario 2 Indices; Dotted Lines Denote the Sample Mean"), width = 9, height = 12,
       device = "pdf")


#Scenario 3 --------------------
load(file="zscoresdf_rho1.RData")
#Histograms of Z-scores -----------------------------------------------------------------------------------------
hist_TotalI <- ggplot(zscoresdf_rho1, aes(x=TotalI)) + labs(y = "Count", x = "Z-score", title="Total I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalI)), color="#20A486", linetype="dashed", size=1)
hist_MeanI <- ggplot(zscoresdf_rho1, aes(x=MeanI)) + labs(y = "Count", x = "Z-score", title="Mean I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)
hist_MedianI <- ggplot(zscoresdf_rho1, aes(x=MedianI)) + labs(y = "Count", x = "Z-score", title="Median I") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MedianI)), color="#20A486", linetype="dashed", size=1)
hist_NormColless <- ggplot(zscoresdf_rho1, aes(x=NormColless)) + labs(y = "Count", x = "Z-score", title="Norm. Colless") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)
hist_CollessLikeExpMDM <- ggplot(zscoresdf_rho1, aes(x=CollessLikeExpMDM)) +  labs(y = "Count", x = "Z-score", title="Colless Like Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(CollessLikeExpMDM)), color="#20A486", linetype="dashed", size=1)
hist_I2 <- ggplot(zscoresdf_rho1, aes(x=I2)) + labs(y = "Count", x = "Z-score", title="I2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(I2)), color="#20A486", linetype="dashed", size=1)
hist_NormSackin <- ggplot(zscoresdf_rho1, aes(x=NormSackin)) + labs(y = "Count", x = "Z-score", title="Norm. Sackin") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormSackin)), color="#20A486", linetype="dashed", size=1)
hist_TotalCophenetic <- ggplot(zscoresdf_rho1, aes(x=TotalCophenetic)) + labs(y = "Count", x = "Z-score", title="Total Cophenetic Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalCophenetic)), color="#20A486", linetype="dashed", size=1)
hist_B1 <- ggplot(zscoresdf_rho1, aes(x=B1)) + labs(y = "Count", x = "Z-score", title="B1") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B1)), color="#20A486", linetype="dashed", size=1)
hist_B2 <- ggplot(zscoresdf_rho1, aes(x=B2)) + labs(y = "Count", x = "Z-score", title="B2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)
hist_AvgLeafDepth <- ggplot(zscoresdf_rho1, aes(x=AvgLeafDepth)) + labs(y = "Count", x = "Z-score", title="Average Leaf Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(AvgLeafDepth)), color="#20A486", linetype="dashed", size=1)
hist_leafdepthvariance <- ggplot(zscoresdf_rho1, aes(x=leafdepthvariance)) + labs(y = "Count", x = "Z-score", title="Leaf Depth Variance") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)
hist_normcherry <- ggplot(zscoresdf_rho1, aes(x=normcherry)) + labs(y = "Count", x = "Z-score", title="Norm. Cherry Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normcherry)), color="#20A486", linetype="dashed", size=1)
hist_Jindex <- ggplot(zscoresdf_rho1, aes(x=Jindex)) + labs(y = "Count", x = "Z-score", title="J Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Jindex)), color="#20A486", linetype="dashed", size=1)
hist_SNI <- ggplot(zscoresdf_rho1, aes(x=SNI)) +  labs(y = "Count", x = "Z-score", title="SNI") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(SNI)), color="#20A486", linetype="dashed", size=1)
hist_Sshape <- ggplot(zscoresdf_rho1, aes(x=Sshape)) +  labs(y = "Count", x = "Z-score", title="S-shape Statistic") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Sshape)), color="#20A486", linetype="dashed", size=1)
hist_APP <- ggplot(zscoresdf_rho1, aes(x=APP)) + labs(y = "Count", x = "Z-score", title="APP") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(APP)), color="#20A486", linetype="dashed", size=1)
hist_normILnumber <- ggplot(zscoresdf_rho1, aes(x=normILnumber)) + labs(y = "Count", x = "Z-score", title="Norm. IL number") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normILnumber)), color="#20A486", linetype="dashed", size=1)
hist_normaverageladder <- ggplot(zscoresdf_rho1, aes(x=normaverageladder)) + labs(y = "Count", x = "Z-score", title="Norm. Avg. Ladder Length") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normaverageladder)), color="#20A486", linetype="dashed", size=1)
hist_maxdepth <- ggplot(zscoresdf_rho1, aes(x=maxdepth)) + labs(y = "Count", x = "Z-score", title="Max Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdepth)), color="#20A486", linetype="dashed", size=1)
hist_maxwidth <- ggplot(zscoresdf_rho1, aes(x=maxwidth)) + labs(y = "Count", x = "Z-score", title="Max Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxwidth)), color="#20A486", linetype="dashed", size=1)
hist_maxdiffwidth <- ggplot(zscoresdf_rho1, aes(x=maxdiffwidth)) + labs(y = "Count", x = "Z-score", title="Max. Diff. in Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdiffwidth)), color="#20A486", linetype="dashed", size=1)
hist_meandepth <- ggplot(zscoresdf_rho1, aes(x=meandepth)) + labs(y = "Count", x = "Z-score", title="Mean Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(meandepth)), color="#20A486", linetype="dashed", size=1)
hist_maxheight <- ggplot(zscoresdf_rho1, aes(x=maxheight)) + labs(y = "Count", x = "Z-score", title="Max Height") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxheight)), color="#20A486", linetype="dashed", size=1)
hist_normpitchforks <- ggplot(zscoresdf_rho1, aes(x=normpitchforks)) + labs(y = "Count", x = "Z-score", title="Norm. # of Pitchforks") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normpitchforks)), color="#20A486", linetype="dashed", size=1)
hist_rootedquartet <- ggplot(zscoresdf_rho1, aes(x=rootedquartet)) + labs(y = "Count", x = "Z-score", title="Rooted Quartet Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(rootedquartet)), color="#20A486", linetype="dashed", size=1)
hist_furnas <- ggplot(zscoresdf_rho1, aes(x=furnas)) + labs(y = "Count", x = "Z-score", title="Furnas Rank") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(furnas)), color="#20A486", linetype="dashed", size=1)
hist_stairs <- ggplot(zscoresdf_rho1, aes(x=stairs)) + labs(y = "Count", x = "Z-score", title="stairs") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs)), color="#20A486", linetype="dashed", size=1)
hist_stairs2 <- ggplot(zscoresdf_rho1, aes(x=stairs2)) + labs(y = "Count", x = "Z-score", title="stairs2") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs2)), color="#20A486", linetype="dashed", size=1)
hist_longestpendant <- ggplot(zscoresdf_rho1, aes(x=longestpendant)) + labs(y = "Count", x = "Z-score", title="Longest Pendant Edge") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(longestpendant)), color="#20A486", linetype="dashed", size=1)
hist_shortestpendant <- ggplot(zscoresdf_rho1, aes(x=shortestpendant)) + labs(y = "Count", x = "Z-score", title="Shortest Pendant Edge") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(shortestpendant)), color="#20A486", linetype="dashed", size=1)
hist_stemminess <- ggplot(zscoresdf_rho1, aes(x=stemminess)) + labs(y = "Count", x = "Z-score", title="Non-Cumulative Stemminess") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)


zhist_all <- list(hist_TotalI, hist_MeanI, hist_MedianI, hist_NormColless, hist_CollessLikeExpMDM, hist_I2, hist_NormSackin, 
                  hist_TotalCophenetic, hist_B1, hist_B2, hist_AvgLeafDepth, hist_leafdepthvariance, hist_normcherry, hist_Jindex, hist_SNI, 
                  hist_Sshape, hist_APP, hist_normILnumber, hist_normaverageladder,  hist_maxdepth, 
                  hist_maxwidth, hist_maxdiffwidth, hist_meandepth, hist_maxheight,
                  hist_normpitchforks, hist_rootedquartet, hist_stairs, hist_stairs2, hist_longestpendant,
                  hist_shortestpendant, hist_stemminess)
ggsave(file="AllZscorePlots_Scenario3.pdf", marrangeGrob(grobs = zhist_all, nrow=4, ncol=3, top = "Z-scores for Scenario 3 Indices; Dotted Lines Denote the Sample Mean"), width = 9, height = 12,
       device = "pdf")


#Scenario 4 --------------------
load(file="zscoresdf_yule.RData")
#Histograms of Z-scores -----------------------------------------------------------------------------------------
hist_TotalI <- ggplot(zscoresdf_yule, aes(x=TotalI)) + labs(y = "Count", x = "Z-score", title="Total I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalI)), color="#20A486", linetype="dashed", size=1)
hist_MeanI <- ggplot(zscoresdf_yule, aes(x=MeanI)) + labs(y = "Count", x = "Z-score", title="Mean I") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MeanI)), color="#20A486", linetype="dashed", size=1)
hist_MedianI <- ggplot(zscoresdf_yule, aes(x=MedianI)) + labs(y = "Count", x = "Z-score", title="Median I") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(MedianI)), color="#20A486", linetype="dashed", size=1)
hist_NormColless <- ggplot(zscoresdf_yule, aes(x=NormColless)) + labs(y = "Count", x = "Z-score", title="Norm. Colless") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormColless)), color="#20A486", linetype="dashed", size=1)
hist_CollessLikeExpMDM <- ggplot(zscoresdf_yule, aes(x=CollessLikeExpMDM)) +  labs(y = "Count", x = "Z-score", title="Colless Like Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(CollessLikeExpMDM)), color="#20A486", linetype="dashed", size=1)
hist_I2 <- ggplot(zscoresdf_yule, aes(x=I2)) + labs(y = "Count", x = "Z-score", title="I2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(I2)), color="#20A486", linetype="dashed", size=1)
hist_NormSackin <- ggplot(zscoresdf_yule, aes(x=NormSackin)) + labs(y = "Count", x = "Z-score", title="Norm. Sackin") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(NormSackin)), color="#20A486", linetype="dashed", size=1)
hist_TotalCophenetic <- ggplot(zscoresdf_yule, aes(x=TotalCophenetic)) + labs(y = "Count", x = "Z-score", title="Total Cophenetic Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(TotalCophenetic)), color="#20A486", linetype="dashed", size=1)
hist_B1 <- ggplot(zscoresdf_yule, aes(x=B1)) + labs(y = "Count", x = "Z-score", title="B1") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B1)), color="#20A486", linetype="dashed", size=1)
hist_B2 <- ggplot(zscoresdf_yule, aes(x=B2)) + labs(y = "Count", x = "Z-score", title="B2") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(B2)), color="#20A486", linetype="dashed", size=1)
hist_AvgLeafDepth <- ggplot(zscoresdf_yule, aes(x=AvgLeafDepth)) + labs(y = "Count", x = "Z-score", title="Average Leaf Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(AvgLeafDepth)), color="#20A486", linetype="dashed", size=1)
hist_leafdepthvariance <- ggplot(zscoresdf_yule, aes(x=leafdepthvariance)) + labs(y = "Count", x = "Z-score", title="Leaf Depth Variance") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(leafdepthvariance)), color="#20A486", linetype="dashed", size=1)
hist_normcherry <- ggplot(zscoresdf_yule, aes(x=normcherry)) + labs(y = "Count", x = "Z-score", title="Norm. Cherry Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normcherry)), color="#20A486", linetype="dashed", size=1)
hist_Jindex <- ggplot(zscoresdf_yule, aes(x=Jindex)) + labs(y = "Count", x = "Z-score", title="J Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Jindex)), color="#20A486", linetype="dashed", size=1)
hist_SNI <- ggplot(zscoresdf_yule, aes(x=SNI)) +  labs(y = "Count", x = "Z-score", title="SNI") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(SNI)), color="#20A486", linetype="dashed", size=1)
hist_Sshape <- ggplot(zscoresdf_yule, aes(x=Sshape)) +  labs(y = "Count", x = "Z-score", title="S-shape Statistic") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(Sshape)), color="#20A486", linetype="dashed", size=1)
hist_APP <- ggplot(zscoresdf_yule, aes(x=APP)) + labs(y = "Count", x = "Z-score", title="APP") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(APP)), color="#20A486", linetype="dashed", size=1)
hist_normILnumber <- ggplot(zscoresdf_yule, aes(x=normILnumber)) + labs(y = "Count", x = "Z-score", title="Norm. IL number") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normILnumber)), color="#20A486", linetype="dashed", size=1)
hist_normaverageladder <- ggplot(zscoresdf_yule, aes(x=normaverageladder)) + labs(y = "Count", x = "Z-score", title="Norm. Avg. Ladder Length") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normaverageladder)), color="#20A486", linetype="dashed", size=1)
hist_maxdepth <- ggplot(zscoresdf_yule, aes(x=maxdepth)) + labs(y = "Count", x = "Z-score", title="Max Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdepth)), color="#20A486", linetype="dashed", size=1)
hist_maxwidth <- ggplot(zscoresdf_yule, aes(x=maxwidth)) + labs(y = "Count", x = "Z-score", title="Max Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxwidth)), color="#20A486", linetype="dashed", size=1)
hist_maxdiffwidth <- ggplot(zscoresdf_yule, aes(x=maxdiffwidth)) + labs(y = "Count", x = "Z-score", title="Max. Diff. in Width") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxdiffwidth)), color="#20A486", linetype="dashed", size=1)
hist_meandepth <- ggplot(zscoresdf_yule, aes(x=meandepth)) + labs(y = "Count", x = "Z-score", title="Mean Depth") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(meandepth)), color="#20A486", linetype="dashed", size=1)
hist_maxheight <- ggplot(zscoresdf_yule, aes(x=maxheight)) + labs(y = "Count", x = "Z-score", title="Max Height") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(maxheight)), color="#20A486", linetype="dashed", size=1)
hist_normpitchforks <- ggplot(zscoresdf_yule, aes(x=normpitchforks)) + labs(y = "Count", x = "Z-score", title="Norm. # of Pitchforks") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(normpitchforks)), color="#20A486", linetype="dashed", size=1)
hist_rootedquartet <- ggplot(zscoresdf_yule, aes(x=rootedquartet)) + labs(y = "Count", x = "Z-score", title="Rooted Quartet Index") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(rootedquartet)), color="#20A486", linetype="dashed", size=1)
hist_furnas <- ggplot(zscoresdf_yule, aes(x=furnas)) + labs(y = "Count", x = "Z-score", title="Furnas Rank") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(furnas)), color="#20A486", linetype="dashed", size=1)
hist_stairs <- ggplot(zscoresdf_yule, aes(x=stairs)) + labs(y = "Count", x = "Z-score", title="stairs") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs)), color="#20A486", linetype="dashed", size=1)
hist_stairs2 <- ggplot(zscoresdf_yule, aes(x=stairs2)) + labs(y = "Count", x = "Z-score", title="stairs2") + geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stairs2)), color="#20A486", linetype="dashed", size=1)
hist_longestpendant <- ggplot(zscoresdf_yule, aes(x=longestpendant)) + labs(y = "Count", x = "Z-score", title="Longest Pendant Edge") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(longestpendant)), color="#20A486", linetype="dashed", size=1)
hist_shortestpendant <- ggplot(zscoresdf_yule, aes(x=shortestpendant)) + labs(y = "Count", x = "Z-score", title="Shortest Pendant Edge") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(shortestpendant)), color="#20A486", linetype="dashed", size=1)
hist_stemminess <- ggplot(zscoresdf_yule, aes(x=stemminess)) + labs(y = "Count", x = "Z-score", title="Non-Cumulative Stemminess") +  geom_histogram(aes(y=..count..),binwidth=0.1, color="#33638D", fill="#33638D") + xlim(-10, 10) + geom_vline(aes(xintercept=mean(stemminess)), color="#20A486", linetype="dashed", size=1)


zhist_all <- list(hist_TotalI, hist_MeanI, hist_MedianI, hist_NormColless, hist_CollessLikeExpMDM, hist_I2, hist_NormSackin, 
                  hist_TotalCophenetic, hist_B1, hist_B2, hist_AvgLeafDepth, hist_leafdepthvariance, hist_normcherry, hist_Jindex, hist_SNI, 
                  hist_Sshape, hist_APP, hist_normILnumber, hist_normaverageladder,  hist_maxdepth, 
                  hist_maxwidth, hist_maxdiffwidth, hist_meandepth, hist_maxheight,
                  hist_normpitchforks, hist_rootedquartet, hist_stairs, hist_stairs2, hist_longestpendant,
                  hist_shortestpendant, hist_stemminess)
ggsave(file="AllZscorePlots_Scenario4.pdf", marrangeGrob(grobs = zhist_all, nrow=4, ncol=3, top = "Z-scores for Scenario 4 Indices; Dotted Lines Denote the Sample Mean"), width = 9, height = 12,
       device = "pdf")







#Supplementary Figure 8, Scree Plots  --------------------------------------------------------
#Scenario 1
tree.pca
scree_plot <- ggscreeplot(tree.pca, type = c("pev")) + labs(y = "Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 1") +
  xlim(1, 5) + theme_bw() + ylim(0,0.6) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))
scree_plot_cumulative <- ggscreeplot(tree.pca, type = "cev") + labs(y = "Cumulative Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 1") +
  xlim(1, 5) + theme_bw() + ylim(0,1) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))

#Scenario 2
tree.pca.nooutgroup
scree_plot_nooutgroup <- ggscreeplot(tree.pca.nooutgroup, type = c("pev")) + labs(y = "Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 2") +
  xlim(1, 5) + theme_bw() + ylim(0,0.6) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))
scree_plot_cumulative_nooutgroup <- ggscreeplot(tree.pca.nooutgroup, type = "cev") + labs(y = "Cumulative Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 2") +
  xlim(1, 5) + theme_bw() + ylim(0,1) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))

#Scenario 3
tree.pca.rho1
scree_plot_rho1 <- ggscreeplot(tree.pca.rho1, type = c("pev")) + labs(y = "Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 3") +
  xlim(1, 5) + theme_bw() + ylim(0,0.6) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))
scree_plot_cumulative_rho1 <- ggscreeplot(tree.pca.rho1, type = "cev") + labs(y = "Cumulative Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 3") +
  xlim(1, 5) + theme_bw() + ylim(0,1) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))

#Scenario 4
tree.pca.yule
scree_plot_yule <- ggscreeplot(tree.pca.yule, type = c("pev")) + labs(y = "Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 4") +
  xlim(1, 5) + theme_bw() + ylim(0,0.6) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))
scree_plot_cumulative_yule <- ggscreeplot(tree.pca.yule, type = "cev") + labs(y = "Cumulative Proportion of Explained Variance", x = "Principal Component Number", title="Scree Plot, Scenario 4") +
  xlim(1, 5) + theme_bw() + ylim(0,1) + theme(axis.title.x=element_text(size = 8), axis.title.y=element_text(size = 8), plot.title = element_text(size=10))

#Arranging into one plot
combined_scree_plots <- ggarrange(scree_plot, scree_plot_cumulative, 
                                  scree_plot_nooutgroup, scree_plot_cumulative_nooutgroup,
                                  scree_plot_rho1, scree_plot_cumulative_rho1, 
                                  scree_plot_yule, scree_plot_cumulative_yule,
                                   heights = c(1,1,1,1),
                                   ncol = 2, nrow = 4)
combined_scree_plots
ggsave(file="Combined_Scree_Plots.pdf", combined_scree_plots, width = 9, height = 12)














