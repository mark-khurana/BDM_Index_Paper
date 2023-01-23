#This file contains code to select the empirical phylogenetic trees to be included from TimeTree, 
#simulation of Scenario 1 simulated trees, and the calculation of index values for the Scenario 1 simulated trees

#Flow of the files: This file -> Scenario_1_Analysis -> Scenario_2_Analysis -> Scenario_3_Simulation_And_Analysis -> Scenario_4_Simulation_And_Analysis -> Final_Figures

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

trs <- list()
#seeing which directory each item belongs to
for(i in 1:length(dir())){
  try(trs[[i]] <- read.tree(dir()[i]))
}

#Finding the class of every attempted reading of a tree, in which the “phylo” classes were successful:
trclass <- sapply(trs, class)
trclass
#Finding which exactly are the problematic ones:
which(trclass != 'phylo')
#Then create a new list with only the good ones:
trgood <- trs[which(trclass == 'phylo')]
numbtrees <- length(trgood)
save(trgood, file="trgood.RData")
load(file="trgood.RData")

#Preparing for the for loop
numbtrees <- length(trgood)
#To get basic characteristics of the trees, such as age, ultrametricity, number of leaves
trfinal <- list("age","ultrametrictrue", "maxbt", "minbt", "realleavesforlength", "gammaStat")
#Final list of trees to be included in the study
trfinallist <- NULL
#Trees that are removed
trout <- NULL

#for loop to find trees that fit the criteria (see article), such as being ultrametric
i <- 1
for (i in i:numbtrees) {
  print(i)
  trgood[[i]] <- collapse.singles(trgood[[i]], root.edge=TRUE) #Collapsing singletons
  tr <- trgood[[i]]
  #Save the rows in the edge object of branches that correspond to tips.
  tips <- which(tr$edge[,2] %in% 1:Ntip(tr))
  #Test whether the number of unique tip lengths is lower than half the number of tips, a bad scenario
  uniquetips <- (length(unique(tr$edge.length[tips])) <= (Ntip(tr) / 2))
  {
    if ((length(trgood[[i]]$edge.length)>0) & (length(trgood[[i]]$tip.label)>3) & (is.ultrametric(trgood[[i]]) == TRUE) & (length(is.binary(trgood[[i]]))==1))
      {
      if (is.binary(trgood[[i]])==TRUE)
        {
        if (uniquetips == F)
        {
      trfinal$ultrametrictrue <- c(trfinal$ultrametrictrue, i)
      trfinal$age <- c(trfinal$maxbt , max(branching.times(trgood[[i]])))
      trfinal$maxbt <- c(trfinal$maxbt , max(branching.times(trgood[[i]])))
      trfinal$minbt <- c(trfinal$minbt, min(branching.times(trgood[[i]])))
      trfinal$gammaStat <- c(trfinal$gammaStat, gammaStat(trgood[[i]]))
      trfinal$realleavesforlength <- c(trfinal$realleavesforlength, length(trgood[[i]]$tip.label))
      
      trfinallist[[length(trfinallist)+ 1]] <- trgood[[i]]
      }
      else
      {
        trout[[length(trout)+ 1]] <- trgood[[i]]
      }
    }
    else
    {
      trout[[length(trout)+ 1]] <- trgood[[i]]
    }
    }
      else
      {
        trout[[length(trout)+ 1]] <- trgood[[i]]
      }
              }
            }

save(trfinal, file="trfinal.RData")
save(trfinallist, file="final_timetrees.RData")
load(file="final_timetrees.RData")

#Finding BD Parameters (Birth, Death, Rho, and corresponding likelihood values) for Each Empirical Tree
numbgoodtrees <- length(trfinallist)
goodtreeparameters <- list("birthparameters", "deathparameters", "likelihoodparameters", "rho")
bestbd <- NULL
i <- 1
for (i in i:numbgoodtrees) {
  print(i)
  {
    fitbdlist <- NULL
    logliklist <- NULL
    n <- 1 
    for (n in n:100) {
      fitbdlist[[n]] <- fit.bd(trfinallist[[i]], rho=((101-n)/100)) #Because we want to select the highest rho value first, assuming multiple maxima; provides a sequence from 1 to 0.01
      logliklist[[n]] <- fitbdlist[[n]]$logL
      }
    max <- NULL
    max <- which.max(logliklist)
    bestbd[[i]] <- fitbdlist[[max]]
    goodtreeparameters$birthparameters <- c(goodtreeparameters$birthparameters, (bestbd[[i]])$b)
    goodtreeparameters$deathparameters <- c(goodtreeparameters$deathparameters, (bestbd[[i]])$d)
    goodtreeparameters$likelihoodparameters <- c(goodtreeparameters$likelihoodparameters, (bestbd[[i]])$logL)
    goodtreeparameters$rho <- c(goodtreeparameters$rho, (bestbd[[i]])$rho)
    }
}

#Birth-death parameters to be used for Scenario 1 (crBDM)
save(goodtreeparameters, file="tree_parameters.RData")
load(file="tree_parameters.RData")

#Simulating Trees from the Parameters -------------------------------------------
library('geiger')
library('TreeSim')

sim.trees <- list()
numbgoodtrees <- length(trfinallist)
numbsim <- 1000

#Creating 1000 simulated trees for each empirical tree
i <- 1
for (i in i:numbgoodtrees) {
  print(i)
  {
    lambda <- goodtreeparameters$birthparameters[i]
    mu <- goodtreeparameters$deathparameters[i]
    n <- trfinal$realleavesforlength[i]
    age <- trfinal$age[i]
    fract <- goodtreeparameters$rho[i]
    
    sim.trees[[i]] <- sim.bd.taxa.age(n=n,numbsim=numbsim,lambda=lambda,mu=mu,frac=fract, age=age, mrca=TRUE) #MRCA true to start from root
  }
}

#Saving simulated tree list--------------------------------------------------------
sim.trees
save(sim.trees, file="simulated_trees.RData")
load(file="simulated_trees.RData")
sim.trees[[1]][1] #to pick out specific simulation
#The list below is the empirical tree list
trfinallist

#Parameters of interest ------------------------------------------------------
install.packages('treebalance')
library('treebalance')
install.packages('phyloTop')
library('phyloTop')
devtools::install_github('Leonardini/treeCentrality')
library('treeCentrality')
library(Claddis)
library(dispRity)
library(phangorn)
library(phytools)
library(apTreeshape)
library(phylobase)
library(parallel)

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


#Indices that could be included in the study, empirical trees ------------------------------------------------------
ntrees <- length(trfinallist)
emptreeparams <- list("age", "maxbt", "minbt", "realleavesforlength", "gammaStat", "TotalI",
                "MeanI", "MedianI", "NormColless", "Collesslike", "I2", "Sackin", "NormalizedSackin", "Totalcophenetic", "B1", "B2", "avgLeafDepth",
                "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape",
                "APP", "ILnumber", "normILnumber", "averageladder", "normaverageladder", "maxdepth", "maxwidth", 
                "maxdiffwidth", "meandepth", "meantopdist", "maxheight", "averagepath", "diameter", "pitchforks",
                "normpitchforks", "wiener", "normwiener", "rootedquartet", "furnas", "stairs", "stairs2", "cladesizes", "longestpendantedge", "shortestpendantedge")

#for loop to find the index values for the empirical trees
i <- 1
for (i in i:ntrees) {
  print(i)
  {
  
    emptreeparams$age <- c(emptreeparams$maxbt, max(branching.times(trfinallist[[i]])))
    emptreeparams$maxbt <- c(emptreeparams$maxbt , max(branching.times(trfinallist[[i]])))
    emptreeparams$minbt <- c(emptreeparams$minbt, min(branching.times(trfinallist[[i]])))
    emptreeparams$gammaStat <- c(emptreeparams$gammaStat, gammaStat(trfinallist[[i]]))
    emptreeparams$realleavesforlength <- c(emptreeparams$realleavesforlength, length(trfinallist[[i]]$tip.label))
    emptreeparams$TotalI <- c(emptreeparams$TotalI, treebalance::IbasedI((trfinallist[[i]]), method="total"))
    emptreeparams$MeanI <- c(emptreeparams$MeanI, treebalance::IbasedI((trfinallist[[i]]), method="mean"))
    emptreeparams$MedianI <- c(emptreeparams$MedianI, treebalance::IbasedI((trfinallist[[i]]), method="median"))
    emptreeparams$NormColless <- c(emptreeparams$NormColless, treebalance::collessI((trfinallist[[i]]), method="corrected"))
    emptreeparams$Collesslike <- c(emptreeparams$Collesslike, treebalance::collesslikeI((trfinallist[[i]]), f.size="exp", dissim="mdm")) #Note settings, exp and median (MDM)
    emptreeparams$I2 <- c(emptreeparams$I2, treebalance::ewCollessI((trfinallist[[i]])))
    emptreeparams$Sackin <- c(emptreeparams$Sackin, treebalance::sackinI((trfinallist[[i]])))
    emptreeparams$NormalizedSackin <- c(emptreeparams$NormalizedSackin, phyloTop::sackin.phylo((trfinallist[[i]]), normalise = TRUE))
    emptreeparams$Totalcophenetic <- c(emptreeparams$Totalcophenetic, treebalance::totCophI((trfinallist[[i]])))
    emptreeparams$B1 <- c(emptreeparams$B1, treebalance::B1I((trfinallist[[i]])))
    emptreeparams$B2 <- c(emptreeparams$B2, treebalance::B2I((trfinallist[[i]])))
    emptreeparams$avgLeafDepth <- c(emptreeparams$avgLeafDepth, treebalance::avgLeafDepI((trfinallist[[i]])))
    emptreeparams$leafdepthvariance <- c(emptreeparams$leafdepthvariance, treebalance::varLeafDepI((trfinallist[[i]])))
    emptreeparams$CPrank <- c(emptreeparams$CPrank, treebalance::colPlaLab((trfinallist[[i]]), method="binary"))
    emptreeparams$normcherry <- c(emptreeparams$normcherry, phyloTop::cherries((trfinallist[[i]]), normalise = TRUE))
    emptreeparams$Jindex <- c(emptreeparams$Jindex, treebalance::rogersI((trfinallist[[i]])))
    emptreeparams$SNI <- c(emptreeparams$SNI, treebalance::symNodesI((trfinallist[[i]])))
    emptreeparams$Sshape <- c(emptreeparams$Sshape, treebalance::sShapeI((trfinallist[[i]])))
    emptreeparams$APP <- c(emptreeparams$APP, treebalance::areaPerPairI((trfinallist[[i]])))
    emptreeparams$ILnumber <- c(emptreeparams$ILnumber, phyloTop::ILnumber((trfinallist[[i]]), normalise = FALSE))
    emptreeparams$normILnumber <- c(emptreeparams$normILnumber, phyloTop::ILnumber((trfinallist[[i]]), normalise = TRUE)) #Normalized for number of tips
    emptreeparams$averageladder <- c(emptreeparams$averageladder, phyloTop::avgLadder((trfinallist[[i]]), normalise = FALSE))
    emptreeparams$normaverageladder <- c(emptreeparams$normaverageladder, phyloTop::avgLadder((trfinallist[[i]]), normalise = TRUE)) #Normalized for number of tips
    emptreeparams$maxdepth <- c(emptreeparams$maxdepth, treebalance::maxDepth((trfinallist[[i]]))) 
    emptreeparams$maxwidth <- c(emptreeparams$maxwidth, treebalance::maxWidth((trfinallist[[i]])))
    emptreeparams$maxdiffwidth <- c(emptreeparams$maxdiffwidth, treebalance::maxDelW((trfinallist[[i]])))
    emptreeparams$meandepth <- c(emptreeparams$meandepth, mean(c(phyloTop::getDepths((trfinallist[[i]]))$tipDepths, phyloTop::getDepths((trfinallist[[i]]))$nodeDepths)))
    emptreeparams$meantopdist <- c(emptreeparams$meandepth, (phyloTop::getDepths((trfinallist[[i]])))$tipdepths)
    emptreeparams$maxheight <- c(emptreeparams$maxheight, phyloTop::maxHeight((trfinallist[[i]])))
    emptreeparams$diameter <- c(emptreeparams$diameter, treeCentrality::computeDiameter((trfinallist[[i]]), weight = TRUE))
    emptreeparams$pitchforks <- c(emptreeparams$pitchforks, phyloTop::pitchforks((trfinallist[[i]]), normalise = FALSE))
    emptreeparams$normpitchforks <- c(emptreeparams$normpitchforks, phyloTop::pitchforks((trfinallist[[i]]), normalise = TRUE)) #Normalized for number of tips
    emptreeparams$rootedquartet <- c(emptreeparams$rootedquartet, treebalance::rQuartetI((trfinallist[[i]])))
    emptreeparams$furnas <- c(emptreeparams$furnas, treebalance::furnasI((trfinallist[[i]])))
    emptreeparams$stairs <- c(emptreeparams$stairs, treebalance::stairs1((trfinallist[[i]])))
    emptreeparams$stairs2 <- c(emptreeparams$stairs2, treebalance::stairs2((trfinallist[[i]])))
    emptreeparams$cladesizes <- c(emptreeparams$cladesizes, treeCentrality::computeNum4to8((trfinallist[[i]])))
    
    #Pendant edges (edges connected to leaves)
    tips <- NULL
    pendantlengths <- NULL
    tips <- which(trfinallist[[i]]$edge[,2] %in% 1:Ntip(trfinallist[[i]]))
    pendantlengths <- (trfinallist[[i]]$edge.length[tips])
    emptreeparams$longestpendantedge <- c(emptreeparams$longestpendantedge, max(pendantlengths))
    emptreeparams$shortestpendantedge <- c(emptreeparams$shortestpendantedge, min(pendantlengths))
  }
}

#Appending stemminess------------------------------------------------------------------------
stemminess_emp <- c()
i <- 1
for (i in i:numbertrees) {
  print(i)
  {
      try(stemminess_emp <- c(stemminess_emp, stemminess_function((trfinallist[[i]]))))
    }
  }

emptreeparams <- append(emptreeparams, "stemminess")
emptreeparams$stemminess <- stemminess_emp

#Saving the index values for the empirical trees
save(emptreeparams, file="emp_tree_indices.RData")
load(file="emp_tree_indices.RData")
    
#Repeating the process to find index values, but this time for the simulated trees for Scenario 1
#Simulated trees
numbertrees <- length(sim.trees)
simtreeparams <- list("age", "maxbt", "minbt", "realleavesforlength", "gammaStat", "TotalI",
                      "MeanI", "MedianI", "NormColless", "Collesslike", "I2", "Sackin", "NormalizedSackin", "Totalcophenetic", "B1", "B2", "avgLeafDepth",
                      "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape",
                      "APP", "ILnumber", "normILnumber", "averageladder", "normaverageladder", "maxdepth", "maxwidth", 
                      "maxdiffwidth", "meandepth", "meantopdist", "maxheight", "maxbetween", "maxclose", "averagepath", "diameter", "pitchforks",
                      "normpitchforks", "wiener","normwiener", "rootedquartet", "furnas", "stairs", "stairs2", "cladesizes", "longestpendantedge", "shortestpendantedge")

#finding index values for the simulated trees, Scenario 1
i <- 1
for (i in i:numbertrees) {
  print(i)
  {
    n <- 1
    for (n in n:1000) {
      simtreeparams$age <- c(simtreeparams$maxbt , max(branching.times(sim.trees[[i]][[n]])))
      simtreeparams$maxbt <- c(simtreeparams$maxbt , max(branching.times(sim.trees[[i]][[n]])))
      simtreeparams$minbt <- c(simtreeparams$minbt, min(branching.times(sim.trees[[i]][[n]])))
      simtreeparams$gammaStat <- c(simtreeparams$gammaStat, gammaStat(sim.trees[[i]][[n]]))
      simtreeparams$realleavesforlength <- c(simtreeparams$realleavesforlength, length(sim.trees[[i]][[n]]$tip.label))
      simtreeparams$TotalI <- c(simtreeparams$TotalI, treebalance::IbasedI((sim.trees[[i]][[n]]), method="total"))
      simtreeparams$MeanI <- c(simtreeparams$MeanI, treebalance::IbasedI((sim.trees[[i]][[n]]), method="mean"))
      simtreeparams$MedianI <- c(simtreeparams$MedianI, treebalance::IbasedI((sim.trees[[i]][[n]]), method="median"))
      simtreeparams$NormColless <- c(simtreeparams$NormColless, treebalance::collessI((sim.trees[[i]][[n]]), method="corrected"))
      simtreeparams$Collesslike <- c(simtreeparams$Collesslike, treebalance::collesslikeI((sim.trees[[i]][[n]]), f.size="exp", dissim="mdm")) #Note settings, exp and median (MDM)
      simtreeparams$I2 <- c(simtreeparams$I2, treebalance::ewCollessI((sim.trees[[i]][[n]])))
      simtreeparams$Sackin <- c(simtreeparams$Sackin, treebalance::sackinI((sim.trees[[i]][[n]])))
      simtreeparams$NormalizedSackin <- c(simtreeparams$NormalizedSackin, phyloTop::sackin.phylo((sim.trees[[i]][[n]]), normalise = TRUE))
      simtreeparams$Totalcophenetic <- c(simtreeparams$Totalcophenetic, treebalance::totCophI((sim.trees[[i]][[n]])))
      simtreeparams$B1 <- c(simtreeparams$B1, treebalance::B1I((sim.trees[[i]][[n]])))
      simtreeparams$B2 <- c(simtreeparams$B2, treebalance::B2I((sim.trees[[i]][[n]])))
      simtreeparams$avgLeafDepth <- c(simtreeparams$avgLeafDepth, treebalance::avgLeafDepI((sim.trees[[i]][[n]])))
      simtreeparams$leafdepthvariance <- c(simtreeparams$leafdepthvariance, treebalance::varLeafDepI((sim.trees[[i]][[n]])))
      simtreeparams$CPrank <- c(simtreeparams$CPrank, treebalance::colPlaLab((sim.trees[[i]][[n]]), method="binary"))
      simtreeparams$normcherry <- c(simtreeparams$normcherry, phyloTop::cherries((sim.trees[[i]][[n]]), normalise = TRUE))
      simtreeparams$Jindex <- c(simtreeparams$Jindex, treebalance::rogersI((sim.trees[[i]][[n]])))
      simtreeparams$SNI <- c(simtreeparams$SNI, treebalance::symNodesI((sim.trees[[i]][[n]])))
      simtreeparams$Sshape <- c(simtreeparams$Sshape, treebalance::sShapeI((sim.trees[[i]][[n]])))
      simtreeparams$APP <- c(simtreeparams$APP, treebalance::areaPerPairI((sim.trees[[i]][[n]])))
      simtreeparams$ILnumber <- c(simtreeparams$ILnumber, phyloTop::ILnumber((sim.trees[[i]][[n]]), normalise = FALSE))
      simtreeparams$normILnumber <- c(simtreeparams$normILnumber, phyloTop::ILnumber((sim.trees[[i]][[n]]), normalise = TRUE))
      simtreeparams$averageladder <- c(simtreeparams$averageladder, phyloTop::avgLadder((sim.trees[[i]][[n]]), normalise = FALSE))
      simtreeparams$normaverageladder <- c(simtreeparams$normaverageladder, phyloTop::avgLadder((sim.trees[[i]][[n]]), normalise = TRUE))
      simtreeparams$maxdepth <- c(simtreeparams$maxdepth, treebalance::maxDepth((sim.trees[[i]][[n]])))
      simtreeparams$maxwidth <- c(simtreeparams$maxwidth, treebalance::maxWidth((sim.trees[[i]][[n]])))
      simtreeparams$maxdiffwidth <- c(simtreeparams$maxdiffwidth, treebalance::maxDelW((sim.trees[[i]][[n]])))
      simtreeparams$meandepth <- c(simtreeparams$meandepth, mean(c(phyloTop::getDepths((sim.trees[[i]][[n]]))$tipDepths, phyloTop::getDepths((sim.trees[[i]][[n]]))$nodeDepths)))
      simtreeparams$meantopdist <- c(simtreeparams$meandepth, (phyloTop::getDepths((sim.trees[[i]][[n]])))$tipdepths)
      simtreeparams$maxheight <- c(simtreeparams$maxheight, phyloTop::maxHeight((sim.trees[[i]][[n]])))
      simtreeparams$diameter <- c(simtreeparams$diameter, treeCentrality::computeDiameter((sim.trees[[i]][[n]]), weight = TRUE))
      simtreeparams$pitchforks <- c(simtreeparams$pitchforks, phyloTop::pitchforks((sim.trees[[i]][[n]]), normalise = FALSE))
      simtreeparams$normpitchforks <- c(simtreeparams$normpitchforks, phyloTop::pitchforks((sim.trees[[i]][[n]]), normalise = TRUE))
      simtreeparams$rootedquartet <- c(simtreeparams$rootedquartet, treebalance::rQuartetI((sim.trees[[i]][[n]])))
      simtreeparams$furnas <- c(simtreeparams$furnas, treebalance::furnasI((sim.trees[[i]][[n]])))
      simtreeparams$stairs <- c(simtreeparams$stairs, treebalance::stairs1((sim.trees[[i]][[n]])))
      simtreeparams$stairs2 <- c(simtreeparams$stairs2, treebalance::stairs2((sim.trees[[i]][[n]])))
      simtreeparams$cladesizes <- c(simtreeparams$cladesizes, treeCentrality::computeNum4to8((sim.trees[[i]][[n]])))
      
    #Pendant edges (edges connected to leaves)
    tips <- NULL
    pendantlengths <- NULL
    tips <- which(sim.trees[[i]][[n]]$edge[,2] %in% 1:Ntip(sim.trees[[i]][[n]]))
    pendantlengths <- (sim.trees[[i]][[n]]$edge.length[tips])
    simtreeparams$longestpendantedge <- c(simtreeparams$longestpendantedge, max(pendantlengths))
    simtreeparams$shortestpendantedge <- c(simtreeparams$shortestpendantedge, min(pendantlengths))
    }
  }
}

#Appending stemminess------------------------------------------------------------------------
i <- 1 #So that all have tip names and all branch lengths incorrectly labeled as negative are positive
for (i in i:numbertrees) {
  print(i)
  {
    trfinallist[[i]]$tip.label[is.na(trfinallist[[i]]$tip.label)] <- "Unknown"
    trfinallist[[i]]$edge.length[0>(trfinallist[[i]]$edge.length)] <- abs(trfinallist[[i]]$edge.length[0>(trfinallist[[i]]$edge.length)])
  }
}

stemminess <- c()
i <- 1
for (i in i:numbertrees) {
  print(i)
  {
    n <- 1
    for (n in n:1000) {
      
      stemminess <- c(stemminess, stemminess_function((sim.trees[[i]][[n]])))
      
    }
  }
}
save(stemminess, file="stemminess_sim.RData")

simtreeparams <- append(simtreeparams, "stemminess")
simtreeparams$stemminess <- stemminess

#Saving the index values for Scenario 1 simulated trees
save(simtreeparams, file="sim_tree_indices.RData")
load(file="sim_tree_indices.RData")

#Creating a data frame with all the data
#Simulated trees
x1<-rep(c(1:1189),each=10) #Because parameter calculations were done 10 at a time (i.e., 100 intervals of n)
x2<-rep(c(x1),times=100)
treenum <- x2
string1 <- rep("sim", 1189000)
simparamdf <- data.frame(treenum, string1, simtreeparams$age, simtreeparams$maxbt, simtreeparams$minbt, simtreeparams$gammaStat, 
                         simtreeparams$realleavesforlength, simtreeparams$TotalI, simtreeparams$MeanI, simtreeparams$MedianI,
                         simtreeparams$NormColless, simtreeparams$Collesslike, simtreeparams$I2, simtreeparams$Sackin,
                         simtreeparams$NormalizedSackin, simtreeparams$Totalcophenetic, simtreeparams$B1, simtreeparams$B2, 
                         simtreeparams$avgLeafDepth, simtreeparams$leafdepthvariance, simtreeparams$CPrank, simtreeparams$normcherry,
                         simtreeparams$Jindex, simtreeparams$SNI, simtreeparams$Sshape, simtreeparams$APP, simtreeparams$ILnumber,
                         simtreeparams$normILnumber, simtreeparams$averageladder, simtreeparams$normaverageladder, simtreeparams$maxdepth,
                         simtreeparams$maxwidth, simtreeparams$maxdiffwidth, simtreeparams$meandepth, simtreeparams$meantopdist,
                         simtreeparams$maxheight, simtreeparams$diameter, simtreeparams$pitchforks, simtreeparams$normpitchforks, 
                         simtreeparams$rootedquartet, simtreeparams$furnas, simtreeparams$stairs, simtreeparams$stairs2, 
                         simtreeparams$longestpendantedge, simtreeparams$shortestpendantedge, simtreeparams$stemminess)
colnames(simparamdf) <- c("treenumber", "EmpOrSim", "age", "maxbt", "minbt", "gamma", "leaves", "TotalI", "MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", "I2", "Sackin", 
                          "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth", "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape", "APP", "ILnumber",
                          "normILnumber", "averageladder", "normaverageladder", "maxdepth", "maxwidth", "maxdiffwidth", "meandepth", "meantopdist", "maxheight", "diameter", "pitchforks",
                          "normpitchforks", "rootedquartet", "furnas", "stairs", "stairs2", "longestpendant", "shortestpendant", "stemminess")

#Empirical trees - not including clade sizes
treenumb <- c(1:1189)
string2 <- rep("emp", 1189)
empparamdf <- data.frame(treenumb, string2, emptreeparams$age, emptreeparams$maxbt, emptreeparams$minbt, emptreeparams$gammaStat, 
                         emptreeparams$realleavesforlength, emptreeparams$TotalI, emptreeparams$MeanI, emptreeparams$MedianI,
                         emptreeparams$NormColless, emptreeparams$Collesslike, emptreeparams$I2, emptreeparams$Sackin,
                         emptreeparams$NormalizedSackin, emptreeparams$Totalcophenetic, emptreeparams$B1, emptreeparams$B2, 
                         emptreeparams$avgLeafDepth, emptreeparams$leafdepthvariance, emptreeparams$CPrank, emptreeparams$normcherry,
                         emptreeparams$Jindex, emptreeparams$SNI, emptreeparams$Sshape, emptreeparams$APP, emptreeparams$ILnumber,
                         emptreeparams$normILnumber, emptreeparams$averageladder, emptreeparams$normaverageladder, emptreeparams$maxdepth,
                         emptreeparams$maxwidth, emptreeparams$maxdiffwidth, emptreeparams$meandepth, emptreeparams$meantopdist,
                         emptreeparams$maxheight, emptreeparams$diameter, emptreeparams$pitchforks, emptreeparams$normpitchforks, 
                         emptreeparams$rootedquartet, emptreeparams$furnas, emptreeparams$stairs, emptreeparams$stairs2, emptreeparams$longestpendantedge, 
                         emptreeparams$shortestpendantedge, emptreeparams$stemminess)
colnames(empparamdf) <- c("treenumber", "EmpOrSim", "age", "maxbt", "minbt", "gamma", "leaves", "TotalI", "MeanI", "MedianI", "NormColless", "CollessLikeExpMDM", "I2", "Sackin", 
                          "NormSackin", "TotalCophenetic", "B1", "B2", "AvgLeafDepth", "leafdepthvariance", "CPrank", "normcherry", "Jindex", "SNI", "Sshape", "APP", "ILnumber",
                          "normILnumber", "averageladder", "normaverageladder", "maxdepth", "maxwidth", "maxdiffwidth", "meandepth", "meantopdist", "maxheight", "diameter", "pitchforks",
                          "normpitchforks", "rootedquartet", "furnas", "stairs", "stairs2", "longestpendant", "shortestpendant", "stemminess")

#One big dataframe with the empirical index values and the index values for the Scenario 1 simulated trees
finaltreeparams <- rbind(empparamdf, simparamdf)
save(finaltreeparams, file="finaltreeparams.RData")
load(file="finaltreeparams.RData")

#The analysis for Scenario 1 can be found in the file: Scenario_1_Analysis 




