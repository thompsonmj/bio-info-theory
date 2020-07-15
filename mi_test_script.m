%% Setup



% load('gap_data_raw_dorsal_wt_XY-Aligned.mat')

% load('psmad_sox2_2020-07-14_13-53.mat')
gene1 = Y_align(:,:,1);
[nX,nEmbryos] = size(gene1);

x = 1/nX:1/nX:1;

nTrials = 1;
nBoots = 1e2;
binCounts = [10:2:50];
subSamps = [0.50, 0.75, 0.80, 0.85, 0.90, 0.95];
% subSamps = 0.95
% binCounts = 50

% genes = {'Hb','Kr','Gt','Kni'};
genes = {'Sox2','P-Smad'};

[MI_est] = ...
    midirectestimate2(Y_align(:,:,[1,2]), nTrials, nBoots, binCounts, subSamps);

extrapolatepi

% while true
%     beep
%     pause(2)
% end

%% Track results
%Hb: 2.1841, 2.1844, 2.1864
%Kr: 1.7353, 1.7366
%Gt: 1.8836, 1.8830
%Kni: 1.6586, 1.6602

%Hb+Kr: 3.4079
%Hb+Gt: 3.4037
%Hb+Kni: 3.1336
%Kr+Gt: 3.0974
%Kr+Kni: 2.8909, 2.8929
%Gt+Kni: 3.2207

%% Plot comparisons (1D)
% Fig. 8 | Table 1 | my calc.
kni = [2.183, 1.83, 1.66]; 
kr = [2.256, 1.96, 1.74]; 
gt = [2.079, 1.84, 1.883]; 
hb = [2.396, 2.24, 2.185];

y = [kni;kr;gt;hb];
x = categorical({'kni', 'kr', 'gt', 'hb'});
bar(x,y)
title('MI Using XY Alignment')
ylabel('MI (bits)')
legend('Fig. 8', 'Table 1', 'my calc.')
grid on

%% Plot comparisons (2D)
% Fig. 12 | my calc.
kni_kr = [3.26, 2.89];
kni_gt = [3.19, 3.22];
kni_hb = [3.15, 3.13];
kr_gt = [3.17, 3.10];
kr_hb = [3.45, 3.41];
gt_hb = [3.34, 3.40];

y = [kni_kr; kni_gt; kni_hb; kr_gt; kr_hb; gt_hb];
x = categorical({'kni/kr', 'kni/gt', 'kni/hb', 'kr/gt', 'kr/hb', 'gt/hb'});
bar(x,y)
title('MI Pairs Using XY Alignment')
ylabel('MI (bits')
legend('Fig. 12', 'my calc.')
grid on

