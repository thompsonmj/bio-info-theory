function [y, alpha, beta, gamma] = alignxy(Y_0, varargin)
%ALIGNXY jointly aligns data along X- and Y- axes of data in a 3D matrix.
% 
% Input:
%   > Y_0: initial profile data. 3D (no. data points) x (no. replicates) x
%   (no. targets).
%   > varargin: NaN pad size as a percentage of the total length. NaN pad
%   is used as a buffer to shift profiles left and right along the X-axis
%   without causing data points to wrap around to the opposite side. If
%   left empty, a default value of 0.05 is used.
% Output:
%   > y: aligned output data. 3D. Same size as input (pads removed).
%   > alpha: additive parameter (1 per profile per target)
%   > beta: scaling/multiplicative parameter (1 per profile per target)
%   > gamma: shifting parameter (1 per replicate, all targets within a
%   replicate are shifted by the same amount)
tic

[nPoints, nReplicates, nTargets] = size(Y_0);

%% Apply NaN padding to either side of each profile.
defaultNanPadSize = 0.05;
switch numel(varargin)
    case 0
        nanPadSize = defaultNanPadSize;
    case 1
        nanPadSize = varargin{1};
end

nNansInPad = round(nPoints * nanPadSize);
Y_0 = padarray(Y_0, nNansInPad, NaN, 'both');

% nNansInPad = sum(isnan(Y_0(:,1,1)))/2;

%% Initialize alignment parameters.
alpha_0 = zeros(nTargets, nReplicates);
beta_0 = ones(nTargets, nReplicates);
gamma_0 = zeros(1, nReplicates);

nAlpha = numel(alpha_0);
nBeta = numel(beta_0);
nGamma = numel(gamma_0);

options = optimoptions('ga', ...
    'UseVectorized', false);

nP = nAlpha + nBeta + nGamma;

%% Specify indices to rereference parameters in cost function.
for iT = 1:nTargets
    alphaStart = nAlpha*(iT/nTargets) - nReplicates + 1;
    alphaEnd = alphaStart + nReplicates - 1;
    idcs.alpha{iT} = alphaStart:alphaEnd;
    
    betaStart = nAlpha + nBeta*(iT/nTargets) - nReplicates + 1;
    betaEnd = betaStart + nReplicates - 1;
    idcs.beta{iT} = betaStart:betaEnd;
end
idcs.gamma = nP - nReplicates + 1: nP;


f = @(P)chisq(P, Y_0, idcs); 

%% Set alignment parameter constraints.
IntCon = (nP - nGamma + 1):nP;

lb_alpha = -ones(nAlpha,1);
ub_alpha = ones(nAlpha,1);

lb_beta = zeros(nBeta,1);
ub_beta = 2*ones(nBeta,1);

lb_gamma = repmat(-nNansInPad/2,nReplicates,1);
ub_gamma = repmat(nNansInPad/2,nReplicates,1);

lb = [lb_alpha; lb_beta; lb_gamma];
ub = [ub_alpha; ub_beta; ub_gamma];

%% Optimization.
% Need to replace GA with something more efficient.
[P_final,~,exitFlag,~,~] = ...
    ga(f,nP,[],[],[],[],lb,ub,[],IntCon,options);

exitFlag

%% Extract parameters from P_final.
alpha = zeros(nTargets,nReplicates);
beta = zeros(nTargets,nReplicates);
gamma = zeros(1,nReplicates);

for iR = 1:nReplicates
    for iT = 1:nTargets
        alpha(iT,iR) = P_final(idcs.alpha{iT}(iR));
        beta(iT,iR) = P_final(idcs.beta{iT}(iR));
    end
    gamma(iR) = P_final(idcs.gamma(iR));
end

%% Apply parameters alpha, beta, and gamma to profile data.
Y = zeros(size(Y_0));
y = zeros(size(Y));
for iR = 1:nReplicates
    Y(:,iR,:) = circshift( Y_0(:,iR,:),gamma(iR),1 );
end
for iT = 1:nTargets
    alpha_tmp = alpha(iT,:);
    beta_tmp = beta(iT,:);
    y(:,:,iT) = (Y(:,:,iT) - alpha_tmp)./beta_tmp;
end

shiftAmt = -round((mean(gamma)));
for iR = 1:nReplicates
    y(:,iR,:) = circshift(y(:,iR,:),shiftAmt,1 );
end

%% Remove NaN padding from either side of each profile.
[nPoints, nReplicates, nTargets] = size(y);

nFrontNans = round( sum(isnan(y(1:round(nPoints/2),:,:)),'all') / ...
    (nReplicates*nTargets) );
nEndNans = round( sum(isnan(y(round(nPoints/2):end,:,:)),'all') / ...
    (nReplicates*nTargets) );

y_temp = y( (nFrontNans : (nPoints - nEndNans + 1)), :, : );

% Ensure the size of the output matches the size of the input.
[n_points_temp,~,~] = size(y_temp);
if n_points_temp == nPoints
    y = y_temp;
elseif n_points_temp == nPoints - 2
    nFrontNans = nFrontNans - 1;
    nEndNans = nEndNans - 1;
    
    y = y( (nFrontNans : (nPoints - nEndNans + 1)), :, : );
elseif n_points_temp == nPoints - 1
    nFrontNans = nFrontNans - 1;
    y = y( (nFrontNans : (nPoints - nEndNans + 1)), :, : );
end

toc

end

%% Cost function
function chi2 = chisq(P, Y, idcs)
    [~,nReplicates,nTargets] = size(Y);

    alpha = zeros(nTargets, nReplicates);
    beta = zeros(nTargets, nReplicates);
    gamma = zeros(1, nReplicates);
        
    for iE = 1:nReplicates
        for iT = 1:nTargets
            alpha(iT,iE) = P(idcs.alpha{iT}(iE));
            beta(iT,iE) = P(idcs.beta{iT}(iE));
        end
        gamma(iE) = P(idcs.gamma(iE));
        Y(:,iE,:) = circshift(Y(:,iE,:), gamma(iE),1);
    end
    
    y_mean = squeeze(nanmean(Y, 2));
    
    chi2_Y = zeros(nTargets, 1);
    for iT = 1:nTargets
        Y_tmp = squeeze(Y(:,:,iT));
        alpha_tmp = alpha(iT,:);
        beta_tmp = beta(iT,:);
        chi2_Y(iT) = ...
            nanmean(nansum( (Y_tmp - (alpha_tmp + beta_tmp.*y_mean(:,iT))).^2 ));
    end
    
    chi2 = mean(chi2_Y);
end