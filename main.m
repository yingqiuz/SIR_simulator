% main.m
% a script to simulate atrophy accrual due to the accumulation of misfolded alpha-syn aggregates
% default
N_regions = 42;
v = 1;
dt = 0.01;
T_total = 10000;
init_number = 1;
syn_control = ROIsize;
prob_stay = 0.5;
trans_rate = 1;
% load your GBA, SNCA, sconnDen, sconnLen, ROISize ....

%%%%% simulation ------ >>>
[Rnor_all, Rmis_all, Rnor0] = SIRsimulator(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate);
ratio = Rmis_all ./(Rnor_all + Rmis_all) ;
ratio(isnan(ratio)) = 0; % remove possible NaNs...

% atrophy growth
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
