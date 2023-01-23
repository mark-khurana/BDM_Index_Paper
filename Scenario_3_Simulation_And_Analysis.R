#Data analysis and simulations for Scenario 3
#To see the code for the final figures and PCA, see the Final_Figures file

setwd("~/Desktop/PhD/TimeTreeData/NwkFiles")

library('ctv')
library('phytools')
library(ape)
library(phangorn)
library(seqinr)
library(adegenet)
library(tidyverse)
library(TreeSim)
library(apTreeshape)
library('treebalance')
library('phyloTop')
library('treeCentrality')
library(Claddis)
library(dispRity)
library(phylobase)
library(parallel)
library(plyr)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(tibble)
library(magrittr)
library(stringr)
library('gginnards')

#Getting the 1189 final trees included in the study
load(file="final_timetrees.RData")

#Finding Parameters for Each Empirical Tree, where Rho = 1 (i.e., where we set the sampling fraction equal to 1)
numbgoodtrees <- length(trfinallist)
goodtreeparameters_rho1 <- list("birthparameters", "deathparameters", "likelihoodparameters", "rho")
bestbd <- NULL
i <- 1
for (i in i:numbgoodtrees) {
  print(i)
  {
    fitbdlist <- NULL
    bestbd[[i]] <- fit.bd(trfinallist[[i]], rho=1)
    goodtreeparameters_rho1$birthparameters <- c(goodtreeparameters_rho1$birthparameters, (bestbd[[i]])$b)
    goodtreeparameters_rho1$deathparameters <- c(goodtreeparameters_rho1$deathparameters, (bestbd[[i]])$d)
    goodtreeparameters_rho1$likelihoodparameters <- c(goodtreeparameters_rho1$likelihoodparameters, (bestbd[[i]])$logL)
    goodtreeparameters_rho1$rho <- c(goodtreeparameters_rho1$rho, (bestbd[[i]])$rho)
  }
}

save(goodtreeparameters_rho1, file="tree_parameters_rho1.RData")

#Simulating 1000 trees for each empirical tree based on the parameters we have inferred above
library('geiger')
library('TreeSim')
load(file="trfinal.RData")

sim.trees_rho1 <- list()

numbgoodtrees <- length(trfinallist)
numbsim <- 1000

i <- 1
for (i in i:numbgoodtrees) {
  print(i)
  {
    lambda <- goodtreeparameters_rho1$birthparameters[i]
    mu <- goodtreeparameters_rho1$deathparameters[i]
    n <- trfinal$realleavesforlength[i]
    age <- trfinal$age[i]
    fract <- goodtreeparameters_rho1$rho[i]
    
    sim.trees_rho1[[i]] <- sim.bd.taxa.age(n=n,numbsim=numbsim,lambda=lambda,mu=mu,frac=fract, age=age, mrca=TRUE) #MRCA true to start from root
  }
}

#Saving the simulated trees for Scenario 3
sim.trees_rho1
save(sim.trees_rho1, file="simulated_trees_rho1.RData")
load(file = "simulated_trees_rho1.RData")

################################################################################


#Changing/Updating the function PhyloCheck in the phylo.top package so it uses is.binary instead of is.binary.tree, since is.binary.tree is outdated and slows the process because of 'error' messages------------
trace(phyloTop::phyloCheck, tracer = quote(is.binary.tree <- is.binary))
trace(ape::is.binary.tree, edit=T) #Removing warning message about is.binary.tree being phased out

#Set up for list with indices --------------------------------------------------------------------------
numbertrees <- length(sim.trees_rho1)
simtreeparams_rho1 <- list("age", "maxbt", "minbt", "realleavesforlength", "gammaStat", "TotalI",
                           "MeanI", "MedianI", "NormColless", "Collesslike", "I2", "Sackin", "NormalizedSackin", "Totalcophenetic", "B1", "B2", "avgLeafDepth",
                           "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape",
                           "APP", "normILnumber", "normaverageladder", "maxdepth", "maxwidth", 
                           "maxdiffwidth", "meandepth", "maxheight", "maxbetween", "maxclose", "averagepath", "diameter",
                           "normpitchforks", "wiener","normwiener", "rootedquartet", "furnas", "stairs", "stairs2", "cladesizes", "longestpendantedge", "shortestpendantedge",
                           "stemminess")

#Finding the index values for the Scenario 3 simulated trees ------------------------------------------------------------------------------
i <- 1
for (i in i:numbertrees) {
  print(i)
  {
    n <- 1
    for (n in n:1000) {
      simtreeparams_rho1$realleavesforlength <- c(simtreeparams_rho1$realleavesforlength, length(sim.trees_rho1[[i]][[n]]$tip.label))
      simtreeparams_rho1$TotalI <- c(simtreeparams_rho1$TotalI, treebalance::IbasedI((sim.trees_rho1[[i]][[n]]), method="total"))
      simtreeparams_rho1$MeanI <- c(simtreeparams_rho1$MeanI, treebalance::IbasedI((sim.trees_rho1[[i]][[n]]), method="mean"))
      simtreeparams_rho1$MedianI <- c(simtreeparams_rho1$MedianI, treebalance::IbasedI((sim.trees_rho1[[i]][[n]]), method="median"))
      simtreeparams_rho1$NormColless <- c(simtreeparams_rho1$NormColless, treebalance::collessI((sim.trees_rho1[[i]][[n]]), method="corrected"))
      simtreeparams_rho1$Collesslike <- c(simtreeparams_rho1$Collesslike, treebalance::collesslikeI((sim.trees_rho1[[i]][[n]]), f.size="exp", dissim="mdm")) #Note settings, exp and median (MDM)
      simtreeparams_rho1$I2 <- c(simtreeparams_rho1$I2, treebalance::ewCollessI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$Sackin <- c(simtreeparams_rho1$Sackin, treebalance::sackinI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$NormalizedSackin <- c(simtreeparams_rho1$NormalizedSackin, phyloTop::sackin.phylo((sim.trees_rho1[[i]][[n]]), normalise = TRUE))
      simtreeparams_rho1$Totalcophenetic <- c(simtreeparams_rho1$Totalcophenetic, treebalance::totCophI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$B1 <- c(simtreeparams_rho1$B1, treebalance::B1I((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$B2 <- c(simtreeparams_rho1$B2, treebalance::B2I((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$avgLeafDepth <- c(simtreeparams_rho1$avgLeafDepth, treebalance::avgLeafDepI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$leafdepthvariance <- c(simtreeparams_rho1$leafdepthvariance, treebalance::varLeafDepI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$CPrank <- c(simtreeparams_rho1$CPrank, treebalance::colPlaLab((sim.trees_rho1[[i]][[n]]), method="binary"))
      simtreeparams_rho1$normcherry <- c(simtreeparams_rho1$normcherry, phyloTop::cherries((sim.trees_rho1[[i]][[n]]), normalise = TRUE))
      simtreeparams_rho1$Jindex <- c(simtreeparams_rho1$Jindex, treebalance::rogersI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$SNI <- c(simtreeparams_rho1$SNI, treebalance::symNodesI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$Sshape <- c(simtreeparams_rho1$Sshape, treebalance::sShapeI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$APP <- c(simtreeparams_rho1$APP, treebalance::areaPerPairI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$normILnumber <- c(simtreeparams_rho1$normILnumber, phyloTop::ILnumber((sim.trees_rho1[[i]][[n]]), normalise = TRUE))
      simtreeparams_rho1$normaverageladder <- c(simtreeparams_rho1$normaverageladder, phyloTop::avgLadder((sim.trees_rho1[[i]][[n]]), normalise = TRUE))
      simtreeparams_rho1$maxdepth <- c(simtreeparams_rho1$maxdepth, treebalance::maxDepth((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$maxwidth <- c(simtreeparams_rho1$maxwidth, treebalance::maxWidth((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$maxdiffwidth <- c(simtreeparams_rho1$maxdiffwidth, treebalance::maxDelW((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$meandepth <- c(simtreeparams_rho1$meandepth, mean(c(phyloTop::getDepths((sim.trees_rho1[[i]][[n]]))$tipDepths, phyloTop::getDepths((sim.trees_rho1[[i]][[n]]))$nodeDepths)))
      simtreeparams_rho1$maxheight <- c(simtreeparams_rho1$maxheight, phyloTop::maxHeight((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$diameter <- c(simtreeparams_rho1$diameter, treeCentrality::computeDiameter((sim.trees_rho1[[i]][[n]]), weight = TRUE))
      simtreeparams_rho1$normpitchforks <- c(simtreeparams_rho1$normpitchforks, phyloTop::pitchforks((sim.trees_rho1[[i]][[n]]), normalise = TRUE))
      simtreeparams_rho1$rootedquartet <- c(simtreeparams_rho1$rootedquartet, treebalance::rQuartetI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$furnas <- c(simtreeparams_rho1$furnas, treebalance::furnasI((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$stairs <- c(simtreeparams_rho1$stairs, treebalance::stairs1((sim.trees_rho1[[i]][[n]])))
      simtreeparams_rho1$stairs2 <- c(simtreeparams_rho1$stairs2, treebalance::stairs2((sim.trees_rho1[[i]][[n]])))

      #Pendant edges (edges connected to leaves)
      tips <- NULL
      pendantlengths <- NULL
      tips <- which(sim.trees_rho1[[i]][[n]]$edge[,2] %in% 1:Ntip(sim.trees_rho1[[i]][[n]]))
      pendantlengths <- (sim.trees_rho1[[i]][[n]]$edge.length[tips])
      simtreeparams_rho1$longestpendantedge <- c(simtreeparams_rho1$longestpendantedge, max(pendantlengths))
      simtreeparams_rho1$shortestpendantedge <- c(simtreeparams_rho1$shortestpendantedge, min(pendantlengths))
    }
  }
}

## Function to compute stemminess index, from Smith et al., DOI: 10.1111/pala.12569 ------------------------------------------------------
stemminess_function <- function(tree){
  Nodes <- Nnode(tree)
  Tips <- Ntip(tree)
  Edges <- tree$edge.length
  tree2 <- phylo4(tree)
  sim_df <- head(tree2, (Tips+Nodes+1))
  stemminess <- data.frame(node = 0, values = 0)
  for(i in (Tips + 2):(Tips + Nodes)){
    dec <- descendants(tree2, i, "all")
    values <- sim_df$edge.length[[i]]
    for(j in dec){
      values <- c(values, sim_df$edge.length[[j]])
    }
    newrow <- data.frame(node = i, values = values[[1]]/sum(values))
    stemminess <- rbind(stemminess, newrow)}
  stemminess <- stemminess[-c(1), ]
  return(mean(stemminess$values))}

#Appending stemminess------------------------------------------------------------------------
stemminess <- c()
i <- 1
for (i in i:numbertrees) {
  print(i)
  {
    n <- 1
    for (n in n:1000) {
      
      stemminess <- c(stemminess, stemminess_function((sim.trees_rho1[[i]][[n]])))
      
    }
  }
}

save(stemminess, file="stemminess_sim_tree_rho1.RData")

simtreeparams <- append(simtreeparams_rho1, "stemminess")
simtreeparams_rho1$stemminess <- stemminess

#saving the final list of index values for Scenario 3
save(simtreeparams_rho1, file="sim_tree_indices_rho1.RData")
load(file="sim_tree_indices_rho1.RData")


#Creating a data frame with all the data, including both empirical and simulated tree index values
#Simulated trees
x1<-rep(c(1:1189),each=100) #Because parameter calculations were done 100 at a time (i.e., 10 intervals of n)
x2<-rep(c(x1),times=10)
treenum <- x2
EmpOrSim <- rep("sim", 1189000)
simparamdf_rho1 <- NULL
simparamdf_rho1 <- data.frame(treenum, EmpOrSim,
                              simtreeparams_rho1$realleavesforlength, simtreeparams_rho1$TotalI, simtreeparams_rho1$MeanI, simtreeparams_rho1$MedianI,
                         simtreeparams_rho1$NormColless, simtreeparams_rho1$Collesslike, simtreeparams_rho1$I2,simtreeparams_rho1$Sackin,
                         simtreeparams_rho1$NormalizedSackin, simtreeparams_rho1$Totalcophenetic, simtreeparams_rho1$B1, simtreeparams_rho1$B2, 
                         simtreeparams_rho1$avgLeafDepth, simtreeparams_rho1$leafdepthvariance, simtreeparams_rho1$CPrank, simtreeparams_rho1$normcherry,
                         simtreeparams_rho1$Jindex, simtreeparams_rho1$SNI, simtreeparams_rho1$Sshape, simtreeparams_rho1$APP,
                         simtreeparams_rho1$normILnumber, simtreeparams_rho1$normaverageladder, simtreeparams_rho1$maxdepth,
                         simtreeparams_rho1$maxwidth, simtreeparams_rho1$maxdiffwidth, simtreeparams_rho1$meandepth, 
                         simtreeparams_rho1$maxheight, simtreeparams_rho1$diameter, simtreeparams_rho1$normpitchforks, 
                         simtreeparams_rho1$rootedquartet, simtreeparams_rho1$furnas, simtreeparams_rho1$stairs, simtreeparams_rho1$stairs2, 
                         simtreeparams_rho1$longestpendantedge, simtreeparams_rho1$shortestpendantedge, simtreeparams_rho1$stemminess)

colnames(simparamdf_rho1) <- c("treenumber", "EmpOrSim", "leaves", "TotalI", "MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", "I2", "Sackin", 
                          "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth", "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape", "APP", 
                          "normILnumber", "normaverageladder", "maxdepth", "maxwidth", "maxdiffwidth", "meandepth", "maxheight", "diameter",
                          "normpitchforks", "rootedquartet", "furnas", "stairs", "stairs2", "longestpendant", "shortestpendant", "stemminess")

#Empirical trees
load(file="emp_tree_indices.RData")
treenum <- c(1:1189)
EmpOrSim <- rep("emp", 1189)
empparamdf_rho1 <- data.frame(treenum, EmpOrSim, 
                         emptreeparams$realleavesforlength, emptreeparams$TotalI, emptreeparams$MeanI, emptreeparams$MedianI,
                         emptreeparams$NormColless, emptreeparams$Collesslike, emptreeparams$I2, emptreeparams$Sackin,
                         emptreeparams$NormalizedSackin, emptreeparams$Totalcophenetic, emptreeparams$B1, emptreeparams$B2, 
                         emptreeparams$avgLeafDepth, emptreeparams$leafdepthvariance, emptreeparams$CPrank, emptreeparams$normcherry,
                         emptreeparams$Jindex, emptreeparams$SNI, emptreeparams$Sshape, emptreeparams$APP,
                         emptreeparams$normILnumber, emptreeparams$normaverageladder, emptreeparams$maxdepth,
                         emptreeparams$maxwidth, emptreeparams$maxdiffwidth, emptreeparams$meandepth,
                         emptreeparams$maxheight, emptreeparams$diameter, emptreeparams$normpitchforks, 
                         emptreeparams$rootedquartet, emptreeparams$furnas, emptreeparams$stairs, emptreeparams$stairs2, emptreeparams$longestpendantedge, 
                         emptreeparams$shortestpendantedge, emptreeparams$stemminess)
colnames(empparamdf_rho1) <- c("treenumber", "EmpOrSim", "leaves", "TotalI", "MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", "I2", "Sackin", 
                          "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth", "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape", "APP",
                          "normILnumber", "normaverageladder", "maxdepth", "maxwidth", "maxdiffwidth", "meandepth", "maxheight", "diameter",
                          "normpitchforks", "rootedquartet", "furnas", "stairs", "stairs2", "longestpendant", "shortestpendant", "stemminess")

#One big dataframe for empirical and simulated tree index values
finaltreeparams_rho1 <- rbind(empparamdf_rho1, simparamdf_rho1)
save(finaltreeparams_rho1, file="finaltreeparams_rho1.RData")
load(file="finaltreeparams_rho1.RData")



#ANALYIS---------------------------------------------------------------------------------------------------------------
#Wilcoxon Rank Sum Test to compare empirical and simulated trees; ------------------------------------------------------------------------------------------------------------------
library('dplyr')
load(file="finaltreeparams_rho1.RData")

empirical <- filter(finaltreeparams_rho1, EmpOrSim == "emp")
simulated <- filter(finaltreeparams_rho1, EmpOrSim == "sim")
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


wil_all_rho1 <- list(wil_TotalI, wil_MeanI, wil_MedianI, wil_NormColless, wil_CollessLikeExpMDM, wil_I2, wil_NormSackin, 
                wil_TotalCophenetic, wil_B1, wil_B2, wil_AvgLeafDepth, wil_leafdepthvariance, wil_CPrank, wil_normcherry, wil_Jindex, wil_SNI, 
                wil_Sshape, wil_APP, wil_normILnumber, wil_normaverageladder,  wil_maxdepth, 
                wil_maxwidth, wil_maxdiffwidth, wil_meandepth, wil_maxheight,
                wil_furnas, wil_normpitchforks, wil_rootedquartet, wil_stairs, wil_stairs2, wil_longestpendant,
                wil_shortestpendant, wil_stemminess)

save(wil_all_rho1, file="wilcoxon_values_rho1.RData")
load(file="wilcoxon_values_rho1.RData")

#Creating a data frame with all the Wilcoxon distance values and p-values
wilcoxon_distance_values_rho1 <- data.frame(
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

save(wilcoxon_distance_values_rho1, file="wilcoxon_distance_values_rho1.RData")
load(file="wilcoxon_distance_values_rho1.RData")

#Z score comparisons-----------------------------------------------------------------------------------------------------------
zscores <- list("treenumber", "TotalI","MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", 
                "I2", "Sackin", "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth",
                "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape",
                "APP", "normILnumber", "normaverageladder", "maxdepth", "maxwidth", 
                "maxdiffwidth", "meandepth", "maxheight",  
                "diameter", "normpitchforks", "rootedquartet", "furnas", 
                "stairs", "stairs2", "longestpendant", "shortestpendant", "stemminess")


i <- 1
for (i in i:1189) {
  print(i)
  {
    tree2 <- NULL
    tree2 <- filter(finaltreeparams_rho1, EmpOrSim == "emp", treenumber == i)
    tree <- NULL
    tree <- filter(finaltreeparams_rho1, EmpOrSim == "sim", treenumber == i)
    
    zscores$treenumber <- c(zscores$treenumber, i)
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
    zscores$longestpendant <- c(zscores$longestpendant,(tree2$longestpendant-(mean(tree$longestpendant)))/sd(tree$longestpendant))
    zscores$shortestpendant <- c(zscores$shortestpendant,(tree2$shortestpendant-(mean(tree$shortestpendant)))/sd(tree$shortestpendant))
    zscores$stemminess <- c(zscores$stemminess,(tree2$stemminess-(mean(tree$stemminess)))/sd(tree$stemminess))
    
  }
}

#Putting all the z-scores in a data frame
zscoresdf_rho1 <- data.frame(zscores$treenumber, zscores$TotalI, zscores$MeanI, zscores$MedianI, zscores$NormColless, 
                        zscores$CollessLikeExpMDM, zscores$I2, zscores$Sackin, zscores$NormSackin, zscores$TotalCophenetic, zscores$B1, 
                        zscores$B2, zscores$AvgLeafDepth,
                        zscores$leafdepthvariance, zscores$CPrank, zscores$normcherry, zscores$Jindex, zscores$SNI, zscores$Sshape,
                        zscores$APP, zscores$normILnumber, 
                        zscores$normaverageladder, zscores$maxdepth, zscores$maxwidth, 
                        zscores$maxdiffwidth, zscores$meandepth, zscores$maxheight,  
                        zscores$normpitchforks, zscores$rootedquartet, zscores$furnas, 
                        zscores$stairs, zscores$stairs2, zscores$longestpendant, zscores$shortestpendant, zscores$stemminess)
colnames(zscoresdf_rho1) <- c("treenumber", "TotalI","MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", 
                         "I2", "Sackin", "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth",
                         "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape",
                         "APP", "normILnumber", "normaverageladder", "maxdepth", "maxwidth", 
                         "maxdiffwidth", "meandepth", "maxheight",  
                         "normpitchforks", "rootedquartet", "furnas", 
                         "stairs", "stairs2", "longestpendant", "shortestpendant", "stemminess")
save(zscoresdf_rho1, file="zscoresdf_rho1.RData")
load(file="zscoresdf_rho1.RData")

#All the final figures used in the paper are coded in the file: Final_Figures






