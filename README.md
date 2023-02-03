## The departure between constant-rate birth-death and empirically inferred diversification processes

This repo contains code to reproduce the simulations and figures of the "The departure between constant-rate birth-death and empirically inferred diversification processes".

## Environment setup:
R version: xxx

The following libraries are used:
```R
library(adegenet)
library(ape)
library(apTreeshape)
library(Claddis)
library(ctv)
library(devtools)
library(dispRity)
library(dplyr)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(githubinstall)
library(ggbiplot)
library(gginnards)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(gridExtra)
library(magrittr)
library(parallel)
library(phangorn)
library(phylobase)
library(phyloTop)
library(phytools)
library(plyr)
library(seqinr)
library(stringr)
library(tibble)
library(tidyverse)
library(treebalance)
library(treeCentrality)
library(TreeSim)
library(viridis)
```

## To reproduce the results:
1. Run ```Tree_Selection_And_1_Trees.R``` to select the empirical phylogenetic trees to be included from TimeTree
2. Run the analysis and simulation scripts for each scenario
3. Run ```Final_Figures.R``` to visualize the results

## Citation

TBD
