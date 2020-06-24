
clear

load('gap_data_raw_dorsal_wt.mat_XY-Aligned.mat')

% Y_aligned = Y_yAlign2;
gene1 = Y_align(:,:,1);
[nX,nEmbryos] = size(gene1);

% padSize = 0.5*sum(sum(isnan(gene1ppp)))/nEmbryos;
% idxLow = padSize + 1;
% idxHigh = nX - padSize;
% Y_aligned = Y_aligned(idxLow:idxHigh,:,:);
% gene1 = Y_aligned(:,:,1);
% [nX,nEmbryos] = size(gene1);
x = 1/nX:1/nX:1;

% d.Hb = Y_aligned(:,:,1);
clear d

% data = d.psmad;
nTrials = 1;
nBoots = 1e2;
binCounts = [10:2:50];
subSamps = [0.50, 0.75, 0.80, 0.85, 0.90, 0.95];


[MI_est] = ...
    midirectestimate2(Y_align(:,:,:), nTrials, nBoots, binCounts, subSamps);

extrapolatepi