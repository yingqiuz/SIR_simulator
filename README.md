## Reference
Zheng, Ying-Qiu, et al. "Local vulnerability and global connectivity jointly shape neurodegenerative disease propagation." PLoS biology 17.11 (2019): e3000495. [(webpage)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000495)

## Overview
This project contains the source code for an agent-based Suspectible-Infectious-Recovered/Removed (SIR) model to predict spreading of pathogenous proteins and atrophy progression in Parkinson's Disease (PD). It contains:

 - SIRsimulator.m: a function to simulate the spread of misfolded alpha-synuclein

 - main.m: a demo script to simulate neuronal loss
 - [data](https://github.com/yingqiuz/SIR_simulator/tree/master/data) folder contains the connectivity and gene expression data required by the model
 - [results](https://github.com/yingqiuz/SIR_simulator/tree/master/results) contains the scripts to generate the figures.

## Getting Started
 - Specify number of regions `N_regions`, velocity `v`, time increment `dt`, total number of steps `T_total`, transmission rate `trans_rate`, seed region `seed`, mobility parameter `prob_stay`. For example, the parameters used in [Zheng et al., 2019](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000495)
```matlab
N_regions = 42;
v = 1;
dt = 0.01;
T_total = 20000;
init_number = 1;
syn_control = ROIsize;
prob_stay = 0.5;
trans_rate = 1;
seed = N_regions;
```
 - Load gene expressions, atrophy, region size, functional and structural connectivity
