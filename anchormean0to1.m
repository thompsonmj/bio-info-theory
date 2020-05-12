function [Y_anchor,Y_mean] = anchormean0to1(Y)
%ANCHORMEAN0TO1 anchors the mean of a profile set from 0 to 1.
%
% Input:
%   > Y: profiles to anchor (nPts x nSamples x nTargets).
%
% Output:
%   > Y_anchor: profiles scaled such that the mean ranges from 0 to 1.
%   > Y_mean: mean profile, min = 0 and max = 1.

Y_mean = nanmean(Y,2);
Y_anchor = (Y - min(Y_mean))./(max(Y_mean) - min(Y_mean));
Y_mean = nanmean(Y_anchor,2);

end