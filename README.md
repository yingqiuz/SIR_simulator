## Reference
Rahayel, Shady, et al. "Differentially targeted seeding reveals unique pathological alpha-synuclein propagation patterns." Brain (2021). [[webpage](https://academic.oup.com/brain/advance-article/doi/10.1093/brain/awab440/6461986?login=true)]

Zheng, Ying-Qiu, et al. "Local vulnerability and global connectivity jointly shape neurodegenerative disease propagation." PLoS biology 17.11 (2019): e3000495. [[webpage](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000495)]

## Overview
This project contains the MATLAB source code for an agent-based Suspectible-Infectious-Recovered/Removed (SIR) model to predict spreading of pathogenous proteins and atrophy progression in Parkinson's Disease (PD). It contains:
 - [SIRsimulator.m](https://github.com/yingqiuz/SIR_simulator/blob/master/SIRsimulator.m): a function to simulate the spread of misfolded alpha-synuclein

 - [main.m](https://github.com/yingqiuz/SIR_simulator/blob/master/main.m): a demo script to simulate neuronal loss
 - [data](https://github.com/yingqiuz/SIR_simulator/tree/master/data) folder contains the connectivity and gene expression data required by the model
 - [results](https://github.com/yingqiuz/SIR_simulator/tree/master/results) contains the scripts to generate the figures.

## Getting Started
The code does not need installation; users have to install MATLAB to run the code. Typically, the model is not computational expensive, but this depends on the size of networks (e.g., number of nodes and edges). For a network of ~100 nodes with 30% connection density, it typically takes less than 5s.
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
 - Load gene expressions, atrophy, region size, functional and structural connectivity on your local system. To use the data in this repository
```matlab
% load gene expressions, real atrophy, ROIsize, functional connectivity...
load('data/42regions/workspace.mat');
% load structural connectivity
load('data/42regions/sc35.mat');
```
This will load GBA and SNCA expressions (`GBA` and `SNCA`), connectivity length `sconnLen`, connectivity strength `sconnDen`, region size `ROIsize` to your workspace
 - Run the model 
```matlab
[Rnor_all, Rmis_all, Rnor0] = SIRsimulator(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate);
```
The region-wise numbers of normal and misfolded proteins are stored in `Rnor_all` and `Rmis_all`, with `Rnor0` being the number of normal proteins at healthy state.
 - Simulate regional atrophy growth. This involves another two parameters to control contribution of deafferentation and endogenous neuronal death. For example, if the two factors have equal contribution:
```matlab
k1 = 0.5;
k2 = 1 - k1;
% input weigths of deafferentation (scaled by structrual connectivity)
weights = sconnDen ./ repmat(sum(sconnDen, 2), 1, N_regions);

% neuronal loss caused by lack of input from neighbouring regions
ratio_cum = weights * (1-exp(-ratio * dt));
% one time step back
ratio_cum = [zeros(N_regions, 1), ratio_cum(:, 1:end-1)];
ratio_cum = k2 * ratio_cum + k1 * (1-exp(-ratio * dt));

% add all the increments across t
simulated_atrophy = cumsum(ratio_cum, 2);
```

## Bugs and Questions
Please contain Ying-Qiu Zheng at [yingqiu.zheng@mail.mcgill.ca](yingqiu.zheng@mail.mcgill.ca)
