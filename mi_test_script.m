%% Setup

clear

load('gap_data_raw_dorsal_wt_XY-Aligned.mat')

gene1 = Y_align(:,:,1);
[nX,nEmbryos] = size(gene1);

x = 1/nX:1/nX:1;

clear d

nTrials = 1;
nBoots = 1e2;
binCounts = [10:2:50];
subSamps = [0.50, 0.75, 0.80, 0.85, 0.90, 0.95];
% subSamps = 0.95
% binCounts = 50

genes = {'Hb','Kr','Gt','Kni'};

[MI_est] = ...
    midirectestimate2(Y_align(:,:,[2,4]), nTrials, nBoots, binCounts, subSamps);

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

%% Plot comparison
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
