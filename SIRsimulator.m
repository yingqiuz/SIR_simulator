% SIRsimulator.m
function [Rnor_all, Rmis_all, Rnor0, Pnor0, Pnor_all, Pmis_all] = SIRsimulator(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate)
% A function to simulate the spread of misfolded alpha-syn

%% input parameters (inside parenthesis are values used in the paper)
% N_regions: number of regions (42)
% v: speed (1)
% dt: time step (0.01)
% T_total: total time steps (10000)
% GBA: GBA gene expression (zscore, N_regions * 1 vector) (empirical GBA expression)
% SNCA: SNCA gene expression after normalization (zscore, N_regions * 1 vector) (empirical SNCA expression)
% sconnLen: structural connectivity matrix (length) (estimated from HCP data)
% sconnDen: structural connectivity matrix (strength) (estimated from HCP data)
% ROIsize: region sizes (voxel counts)
% seed: seed region of misfolded alpha-syn injection (choose as you like? (^?^)= here substantia nigra)
% syn_control: a parameter to control the number of voxels in which
% alpha-syn may get synthesized (region size, i.e., ROIsize)
% init_number: number of injected misfolded alpha-syn (1)
% prob_stay: the probability of staying in the same region per unit time (0.5)
% trans_rate: a scalar value, controlling the baseline infectivity

%% output parameters
% Rnor_all: A N_regions * T_total matrix, recording the number of normal
% alpha-syn in regions
% Rmis_all: A N_regions * T_total matrix, recording the number of
% misfolded alph-syn in regions
% Pnor_all: a N_regions * N_regions * T_total matrix, recording the number of normal alpha-syn in paths
% could be memory-consuming)
% Pmis_all: a N_regions * N_regions * T_total matrix, recording the number of misfolded alpha-syn in paths
% could be memory-consuming)
% Rnor0: a N_Regions * 1 vector, the population of normal agents in regions before pathogenic spreading 
% Pnor0: a N_Regions * 1 vecotr, the population of normal agents in edges before pathogenic spreading 

% make sure the diag is zero
sconnDen(eye(N_regions)==1) = 0;
sconnLen(eye(N_regions)==1) = 0;

% set the mobility pattern
weights = sconnDen;
weights = (1 - prob_stay) .* weights + prob_stay .* diag(sum(weights, 2)) ;

% multinomial distribution
% element (i,j) is the probability of moving from region i to edge (i,j)
weights = weights ./ repmat(sum(weights, 2), 1, N_regions);

% convert gene expression scores to probabilities
clearance_rate = normcdf(zscore(GBA)); 
synthesis_rate = normcdf(zscore(SNCA));

% store the number of normal/misfoled alpha-syn at each time step
[Rnor_all, Rmis_all] = deal( zeros([N_regions, T_total]) );
[Pnor_all, Pmis_all] = deal( zeros([N_regions, N_regions, T_total]) );

% Rnor, Rmis, Pnor, Pmis store results of single simulation at each time
[Rnor, Rmis] = deal(zeros(N_regions, 1)); % number of normal/misfolded alpha-syn in regions
[Pnor, Pmis] = deal(zeros(N_regions)); % number of normal/misfolded alpha-syn in paths

%% normal alpha-syn growth 
% fill the network with normal proteins
iter_max = 1000000000;
display('normal alpha synuclein growth')
for t = 1:iter_max   
    %%% moving process
    % regions towards paths
    % movDrt stores the number of proteins towards each region. i.e.
    % element in kth row lth col denotes the number of proteins in region k
    % moving towards l
    movDrt = repmat(Rnor, 1, N_regions) .* weights;
    movDrt = movDrt .* dt; 
    movDrt(eye(N_regions)==1) = 0;
    
    % paths towards regions
    % update moving
    movOut = Pnor .* v ./ sconnLen ;  % longer path & smaller v = lower probability of moving out of paths
    movOut(sconnLen == 0) = 0;
    
    Pnor = Pnor - movOut .* dt + movDrt;
    Pnor(eye(N_regions)==1) = 0;
    
    Rtmp = Rnor;
    Rnor = Rnor + sum(movOut, 1)' .* dt -  sum(movDrt, 2);
    
    %%% growth process
    Rnor = Rnor -  Rnor.* (1-exp(-clearance_rate.*dt))   + (synthesis_rate .* syn_control).*dt ;
    
    if abs(Rnor - Rtmp) < (1e-7* Rtmp)
        break;
    end
end

Pnor0 = Pnor;
Rnor0 = Rnor;

%% misfolded protein spreading process

% inject misfolded alpha-syn
Rmis(seed) = init_number;
display('misfolded alpha synuclein spreading')
for t = 1:T_total
    %%% moving process
    % normal proteins: region -->> paths
    movDrt_nor = repmat(Rnor, 1, N_regions) .* weights .* dt;
    movDrt_nor(eye(N_regions) == 1) = 0;
    
    % normal proteins: paths -->> regions
    movOut_nor = Pnor .* v  ./ sconnLen; 
    movOut_nor(sconnLen == 0) = 0;
    
    
    % misfolded proteins: region -->> paths
    movDrt_mis = repmat(Rmis, 1, N_regions) .* weights .* dt;
    movDrt_mis(eye(N_regions) == 1) = 0;
    
    % misfolded proteins: paths -->> regions
    movOut_mis = Pmis .* v ./ sconnLen; 
    movOut_mis(sconnLen == 0) = 0;
    
    % update regions and paths
    Pnor = Pnor - movOut_nor.* dt + movDrt_nor; 
    Pnor(eye(N_regions)==1) = 0;
    Rnor = Rnor + sum(movOut_nor, 1)'.* dt - sum(movDrt_nor, 2);
    
    Pmis = Pmis - movOut_mis.*dt + movDrt_mis; 
    Pmis(eye(N_regions)==1) = 0;
    Rmis = Rmis + sum(movOut_mis, 1)'.*dt - sum(movDrt_mis, 2);    
            
    Rnor_cleared = Rnor .* (1-exp(-clearance_rate.* dt)) ;
    Rmis_cleared = Rmis .* (1-exp(-clearance_rate.* dt))  ;
    % the probability of getting misfolded
    gamma0 = 1 .* trans_rate ./ROIsize;
    misProb = 1 - exp( -Rmis .* gamma0 .* dt ) ; % trans_rate: default
    % number of newly infected
    N_misfolded = Rnor .* exp(-clearance_rate) .* misProb ;
    % update
    Rnor = Rnor - Rnor_cleared - N_misfolded + (synthesis_rate .* syn_control).*dt;
    Rmis = Rmis - Rmis_cleared + N_misfolded;

    Rnor_all(:, t) = Rnor ;
    Rmis_all(:, t) = Rmis ;
    
    % uncomment the following lines if you want outputs of alpha-syn in
    % paths
    %Pnor_ave(:, :, t) = Pnor;
    %Pmis_ave(:, :, t) = Pmis;
   
end
end




