#Data analysis for Scenario 2
#To see the code for the final figures and PCA, see the Final_Figures file

setwd()

#Loading packages
library('ctv')
library('phytools')
library(ape)
library(phangorn)
library(seqinr)
library(adegenet)
library(tidyverse)
library(TreeSim)
library(apTreeshape)
library(phytools)
library('treebalance')
library('phyloTop')
devtools::install_github('Leonardini/treeCentrality')
library('treeCentrality')
library(devtools)
library(plyr)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(tibble)
library(magrittr)
library(stringr)
library('gginnards')

#Loading empirical trees
load(file="finaltreeparams.RData")

#Excluding trees with pendants that extend all the way (i.e., those with potential outgroups)--------------------------------
#Use a which statement to exclude trees with potential outgroups
nooutgroup <- filter(finaltreeparams, longestpendant != age, EmpOrSim == "emp")
treenums <- nooutgroup$treenumber
nooutgrouptrees <- finaltreeparams %>% filter(treenumber %in% treenums)
nooutgrouptrees

#Wilcoxon Rank Sum Test to compare empirical and simulated trees; no need to simulate new trees before analysis
library('dplyr')
empirical <- filter(nooutgrouptrees, EmpOrSim == "emp")
simulated <- filter(nooutgrouptrees, EmpOrSim == "sim")
wil_TotalI <- wilcox.test(empirical$TotalI, simulated$TotalI, paired = FALSE, alternative = "two.sided")
wil_MeanI <- wilcox.test(empirical$MeanI, simulated$MeanI, paired = FALSE, alternative = "two.sided")
wil_MedianI <- wilcox.test(empirical$MedianI, simulated$MedianI, paired = FALSE, alternative = "two.sided")
wil_NormColless <- wilcox.test(empirical$NormColless, simulated$NormColless, paired = FALSE, alternative = "two.sided")
wil_CollessLikeExpMDM <- wilcox.test(empirical$CollessLikeExpMDM, simulated$CollessLikeExpMDM, paired = FALSE, alternative = "two.sided")
wil_I2 <- wilcox.test(empirical$I2, simulated$I2, paired = FALSE, alternative = "two.sided")
wil_NormSackin <- wilcox.test(empirical$NormSackin, simulated$NormSackin, paired = FALSE, alternative = "two.sided")
wil_TotalCophenetic <- wilcox.test(empirical$TotalCophenetic, simulated$TotalCophenetic, paired = FALSE, alternative = "two.sided")
wil_B1 <- wilcox.test(empirical$B1, simulated$B1, paired = FALSE, alternative = "two.sided")
wil_B2 <- wilcox.test(empirical$B2, simulated$B2, paired = FALSE, alternative = "two.sided")
wil_AvgLeafDepth <- wilcox.test(empirical$AvgLeafDepth, simulated$AvgLeafDepth, paired = FALSE, alternative = "two.sided")
wil_leafdepthvariance <- wilcox.test(empirical$leafdepthvariance, simulated$leafdepthvariance, paired = FALSE, alternative = "two.sided")
wil_CPrank <- wilcox.test(empirical$CPrank, simulated$CPrank, paired = FALSE, alternative = "two.sided")
wil_normcherry <- wilcox.test(empirical$normcherry, simulated$normcherry, paired = FALSE, alternative = "two.sided")
wil_Jindex <- wilcox.test(empirical$Jindex, simulated$Jindex, paired = FALSE, alternative = "two.sided")
wil_SNI <- wilcox.test(empirical$SNI, simulated$SNI, paired = FALSE, alternative = "two.sided")
wil_Sshape <- wilcox.test(empirical$Sshape, simulated$Sshape, paired = FALSE, alternative = "two.sided")
wil_APP <- wilcox.test(empirical$APP, simulated$APP, paired = FALSE, alternative = "two.sided")
wil_normILnumber <- wilcox.test(empirical$normILnumber, simulated$normILnumber, paired = FALSE, alternative = "two.sided")
wil_normaverageladder <- wilcox.test(empirical$normaverageladder, simulated$normaverageladder, paired = FALSE, alternative = "two.sided")
wil_maxdepth <- wilcox.test(empirical$maxdepth, simulated$maxdepth, paired = FALSE, alternative = "two.sided")
wil_maxwidth <- wilcox.test(empirical$maxwidth, simulated$maxwidth, paired = FALSE, alternative = "two.sided")
wil_maxwidthtodepth <- wilcox.test((empirical$maxwidth/empirical$maxdepth), (simulated$maxwidth/simulated$maxdepth), paired = FALSE, alternative = "two.sided")
wil_maxdiffwidth <- wilcox.test(empirical$maxdiffwidth, simulated$maxdiffwidth, paired = FALSE, alternative = "two.sided")
wil_meandepth <- wilcox.test(empirical$meandepth, simulated$meandepth, paired = FALSE, alternative = "two.sided")
wil_maxheight <- wilcox.test(empirical$maxheight, simulated$maxheight, paired = FALSE, alternative = "two.sided")
wil_diameter <- wilcox.test(empirical$diameter, simulated$diameter, paired = FALSE, alternative = "two.sided")
wil_normpitchforks <- wilcox.test(empirical$normpitchforks, simulated$normpitchforks, paired = FALSE, alternative = "two.sided")
wil_rootedquartet <- wilcox.test(empirical$rootedquartet, simulated$rootedquartet, paired = FALSE, alternative = "two.sided")
wil_furnas <- wilcox.test(empirical$furnas, simulated$furnas, paired = FALSE, alternative = "two.sided")
wil_stairs <- wilcox.test(empirical$stairs, simulated$stairs, paired = FALSE, alternative = "two.sided")
wil_stairs2 <- wilcox.test(empirical$stairs2, simulated$stairs2, paired = FALSE, alternative = "two.sided")
wil_longestpendant <- wilcox.test(empirical$longestpendant, simulated$longestpendant, paired = FALSE, alternative = "two.sided")
wil_shortestpendant <- wilcox.test(empirical$shortestpendant, simulated$shortestpendant, paired = FALSE, alternative = "two.sided")
wil_stemminess <- wilcox.test(empirical$stemminess, simulated$stemminess, paired = FALSE, alternative = "two.sided")


wil_all_nooutgroup <- list(wil_TotalI, wil_MeanI, wil_MedianI, wil_NormColless, wil_CollessLikeExpMDM, wil_I2, wil_NormSackin, 
                wil_TotalCophenetic, wil_B1, wil_B2, wil_AvgLeafDepth, wil_leafdepthvariance, wil_CPrank, wil_normcherry, wil_Jindex, wil_SNI, 
                wil_Sshape, wil_APP, wil_normILnumber, wil_normaverageladder,  wil_maxdepth, 
                wil_maxwidth, wil_maxdiffwidth, wil_meandepth, wil_maxheight,
                wil_furnas, wil_normpitchforks, wil_rootedquartet, wil_stairs, wil_stairs2, wil_longestpendant,
                wil_shortestpendant, wil_stemminess)

save(wil_all_nooutgroup, file="wilcoxon_values_nooutgroup.RData")
load(file="wilcoxon_values_nooutgroup.RData")

#Creating a data frame with all the Wilcoxon distance values and p-values
wilcoxon_distance_values_nooutgroup <- data.frame(
  name=c("Total I","Mean I", "Median I", "Norm. Colless", "Colless-like index", 
         "I2", "Norm. Sackin", "Total coph. index", "B1", "B2", "Avg. leaf depth",
         "Leaf depth variance", "Norm. cherries", "J index", "SNI", "S-shape statistic",
         "APP", "Norm. IL number", "Norm. avg. ladder length", "Max. depth", "Max. width", 
         "Max. diff. in widths", "Mean depth", "Max. height",  
         "Norm. pitchforks", "Rooted Quartet", "Furnas rank", 
         "stairs", "stairs2", "Longest pendant edge", "Shortest pendant edge", "Stemminess"),
  distance_value=c(wil_TotalI$statistic, wil_MeanI$statistic, wil_MedianI$statistic, wil_NormColless$statistic, wil_CollessLikeExpMDM$statistic, 
                   wil_I2$statistic, wil_NormSackin$statistic, 
                   wil_TotalCophenetic$statistic, wil_B1$statistic, wil_B2$statistic, wil_AvgLeafDepth$statistic, wil_leafdepthvariance$statistic, 
                   wil_normcherry$statistic, wil_Jindex$statistic, wil_SNI$statistic, 
                   wil_Sshape$statistic, wil_APP$statistic, wil_normILnumber$statistic, wil_normaverageladder$statistic,  wil_maxdepth$statistic, 
                   wil_maxwidth$statistic, wil_maxdiffwidth$statistic, wil_meandepth$statistic, wil_maxheight$statistic,
                   wil_normpitchforks$statistic, wil_rootedquartet$statistic, wil_furnas$statistic, wil_stairs$statistic, wil_stairs2$statistic, 
                   wil_longestpendant$statistic,
                   wil_shortestpendant$statistic, wil_stemminess$statistic),
  p_value=c(wil_TotalI$p.value, wil_MeanI$p.value, wil_MedianI$p.value, wil_NormColless$p.value, wil_CollessLikeExpMDM$p.value, 
            wil_I2$p.value, wil_NormSackin$p.value, 
            wil_TotalCophenetic$p.value, wil_B1$p.value, wil_B2$p.value, wil_AvgLeafDepth$p.value, wil_leafdepthvariance$p.value, 
            wil_normcherry$p.value, wil_Jindex$p.value, wil_SNI$p.value, 
            wil_Sshape$p.value, wil_APP$p.value, wil_normILnumber$p.value, wil_normaverageladder$p.value,  wil_maxdepth$p.value, 
            wil_maxwidth$p.value, wil_maxdiffwidth$p.value, wil_meandepth$p.value, wil_maxheight$p.value,
            wil_normpitchforks$p.value, wil_rootedquartet$p.value, wil_furnas$p.value, wil_stairs$p.value, wil_stairs2$p.value, wil_longestpendant$p.value,
            wil_shortestpendant$p.value, wil_stemminess$p.value))

save(wilcoxon_distance_values_nooutgroup, file="wilcoxon_distance_values_nooutgroup.RData")
load(file="wilcoxon_distance_values_nooutgroup.RData")

#Calculating the z-score for each index, where the mean and SD are taken from the simulated distribution but the z-score is calculated for the empirical tree
zscores <- list("age", "gamma", "TotalI","MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", 
                "I2", "Sackin", "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth",
                "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape",
                "APP", "normILnumber", "normaverageladder", "maxdepth", "maxwidth", 
                "maxdiffwidth", "meandepth", "maxheight",  
                "diameter", "normpitchforks", "rootedquartet", "furnas", 
                "stairs", "stairs2", "longestpendant", "shortestpendant", "longestpendanttt", "stemminess")

#For loop to get the z-score values
i <- 1
for (i in i:1189) {
  print(i)
  {
    tree2 <- NULL
    tree2 <- filter(finaltreeparams, EmpOrSim == "emp", treenumber == i)
    tree <- NULL
    tree <- filter(finaltreeparams, EmpOrSim == "sim", treenumber == i)
    
    zscores$age <- c(zscores$age,tree2$age)
    zscores$longestpendanttt <- c(zscores$longestpendanttt,tree2$longestpendant)
    
    zscores$gamma <- c(zscores$gamma,(tree2$gamma-(mean(tree$gamma)))/sd(tree$gamma))
    zscores$TotalI <- c(zscores$TotalI,(tree2$TotalI-(mean(tree$TotalI)))/sd(tree$TotalI))
    zscores$MeanI <- c(zscores$MeanI,(tree2$MeanI-(mean(tree$MeanI)))/sd(tree$MeanI))
    zscores$MedianI <- c(zscores$MedianI,(tree2$MedianI-(mean(tree$MedianI)))/sd(tree$MedianI))
    zscores$NormColless <- c(zscores$NormColless,(tree2$NormColless-(mean(tree$NormColless)))/sd(tree$NormColless))
    zscores$CollessLikeExpMDM <- c(zscores$CollessLikeExpMDM,(tree2$CollessLikeExpMDM-(mean(tree$CollessLikeExpMDM)))/sd(tree$CollessLikeExpMDM))
    zscores$I2 <- c(zscores$I2,(tree2$I2-(mean(tree$I2)))/sd(tree$I2))
    zscores$Sackin <- c(zscores$Sackin,(tree2$Sackin-(mean(tree$Sackin)))/sd(tree$Sackin))
    zscores$NormSackin <- c(zscores$NormSackin,(tree2$NormSackin-(mean(tree$NormSackin)))/sd(tree$NormSackin))
    zscores$TotalCophenetic <- c(zscores$TotalCophenetic,(tree2$TotalCophenetic-(mean(tree$TotalCophenetic)))/sd(tree$TotalCophenetic))
    zscores$B1 <- c(zscores$B1,(tree2$B1-(mean(tree$B1)))/sd(tree$B1))
    zscores$B2 <- c(zscores$B2,(tree2$B2-(mean(tree$B2)))/sd(tree$B2))
    zscores$AvgLeafDepth <- c(zscores$AvgLeafDepth,(tree2$AvgLeafDepth-(mean(tree$AvgLeafDepth)))/sd(tree$AvgLeafDepth))
    zscores$leafdepthvariance <- c(zscores$leafdepthvariance,(tree2$leafdepthvariance-(mean(tree$leafdepthvariance)))/sd(tree$leafdepthvariance))
    zscores$CPrank <- c(zscores$CPrank,(tree2$CPrank-(mean(tree$CPrank)))/sd(tree$CPrank))
    zscores$normcherry <- c(zscores$normcherry,(tree2$normcherry-(mean(tree$normcherry)))/sd(tree$normcherry))
    zscores$Jindex <- c(zscores$Jindex,(tree2$Jindex-(mean(tree$Jindex)))/sd(tree$Jindex))
    zscores$SNI <- c(zscores$SNI,(tree2$SNI-(mean(tree$SNI)))/sd(tree$SNI))
    zscores$Sshape <- c(zscores$Sshape,(tree2$Sshape-(mean(tree$Sshape)))/sd(tree$Sshape))
    zscores$APP <- c(zscores$APP,(tree2$APP-(mean(tree$APP)))/sd(tree$APP))
    zscores$normILnumber <- c(zscores$normILnumber,(tree2$normILnumber-(mean(tree$normILnumber)))/sd(tree$normILnumber))
    zscores$normaverageladder <- c(zscores$normaverageladder,(tree2$normaverageladder-(mean(tree$normaverageladder)))/sd(tree$normaverageladder))
    zscores$maxdepth <- c(zscores$maxdepth,(tree2$maxdepth-(mean(tree$maxdepth)))/sd(tree$maxdepth))
    zscores$maxwidth <- c(zscores$maxwidth,(tree2$maxwidth-(mean(tree$maxwidth)))/sd(tree$maxwidth))
    zscores$maxdiffwidth <- c(zscores$maxdiffwidth,(tree2$maxdiffwidth-(mean(tree$maxdiffwidth)))/sd(tree$maxdiffwidth))
    zscores$meandepth <- c(zscores$meandepth,(tree2$meandepth-(mean(tree$meandepth)))/sd(tree$meandepth))
    zscores$maxheight <- c(zscores$maxheight,(tree2$maxheight-(mean(tree$maxheight)))/sd(tree$maxheight))
    zscores$averagepath <- c(zscores$averagepath,(tree2$averagepath-(mean(tree$averagepath)))/sd(tree$averagepath))
    zscores$normpitchforks <- c(zscores$normpitchforks,(tree2$normpitchforks-(mean(tree$normpitchforks)))/sd(tree$normpitchforks))
    zscores$rootedquartet <- c(zscores$rootedquartet,(tree2$rootedquartet-(mean(tree$rootedquartet)))/sd(tree$rootedquartet))
    zscores$furnas <- c(zscores$furnas,(tree2$furnas-(mean(tree$furnas)))/sd(tree$furnas))
    zscores$stairs <- c(zscores$stairs,(tree2$stairs-(mean(tree$stairs)))/sd(tree$stairs))
    zscores$stairs2 <- c(zscores$stairs2,(tree2$stairs2-(mean(tree$stairs2)))/sd(tree$stairs2))
    zscores$shortestpendant <- c(zscores$shortestpendant,(tree2$shortestpendant-(mean(tree$shortestpendant)))/sd(tree$shortestpendant))
    zscores$stemminess <- c(zscores$stemminess,(tree2$stemminess-(mean(tree$stemminess)))/sd(tree$stemminess))
    
     }
}


#Putting all the z-scores in a data frame
zscoresdf <- data.frame(zscores$age, zscores$longestpendanttt, zscores$gamma, zscores$TotalI, zscores$MeanI, zscores$MedianI, zscores$NormColless, 
                        zscores$CollessLikeExpMDM, zscores$I2, zscores$Sackin, zscores$NormSackin, zscores$TotalCophenetic, zscores$B1, 
                        zscores$B2, zscores$AvgLeafDepth,
                        zscores$leafdepthvariance, zscores$CPrank, zscores$normcherry, zscores$Jindex, zscores$SNI, zscores$Sshape,
                        zscores$APP, zscores$normILnumber, 
                        zscores$normaverageladder, zscores$maxdepth, zscores$maxwidth, 
                        zscores$maxdiffwidth, zscores$meandepth, zscores$maxheight,  
                        zscores$normpitchforks, zscores$rootedquartet, zscores$furnas, 
                        zscores$stairs, zscores$stairs2, zscores$shortestpendant, zscores$stemminess)
colnames(zscoresdf) <- c("age", "longestpendanttt", "gamma", "TotalI","MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", 
                         "I2", "Sackin", "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth",
                         "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape",
                         "APP", "normILnumber", "normaverageladder", "maxdepth", "maxwidth", 
                         "maxdiffwidth", "meandepth", "maxheight",  
                         "normpitchforks", "rootedquartet", "furnas", 
                         "stairs", "stairs2", "shortestpendant", "stemminess")
zscoresdf_nooutgroupdf <- filter(zscoresdf, longestpendanttt != age)

save(zscoresdf_nooutgroupdf, file="zscoresdf_nooutgroup.RData")
load(file="zscoresdf_nooutgroup.RData")


#All the final figures used in the paper are coded in the file: Final_Figures


































