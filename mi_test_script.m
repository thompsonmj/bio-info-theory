%% Setup
oworkdir = pwd;
load('gap_data_raw_dorsal_wt.mat')
nX = numel(data(1).Hb(1,:)); 
nEmb = numel(data); % Number of embryos (N in text).
d.Hb = zeros(nX,nEmb);
for iEmb = 1:nEmb
    d.Hb(:,iEmb) = transpose( data(iEmb).Hb(1:nX) );
end
d.Hb = d.Hb(nX*0.1:nX*0.9,:);
meanHb = mean(d.Hb,2);
d.Hb = (d.Hb - min(meanHb))/(max(meanHb)-min(meanHb));
meanHb = mean(d.Hb,2);
x = 0.1:0.001:0.9;
nX = numel(x);

% normout = normdat(d.Hb);
% d.Hb = normout.chiSq;

profileData = d.Hb;
xData = x;
nTrials = 3;
nBoots = 100;
% binCountArray = [10:2:50];
% subsampleArray = [0.50, 0.75, 0.80, 0.85, 0.90, 0.95];

binCountArray = [40:2:50];
subsampleArray = [0.85, 0.95];

I_est = ...
    midirectestimate(profileData, xData, nTrials, nBoots, binCountArray, subsampleArray);