function I_est = ...
    minaiveestimate(yData,xData,nTrials,nBoots,binCountArray,subsampleArray)
% MINAIVEESTIMATE computes a naive estimate of mutual information.
% 
% Inputs
%   > yData:
%   > xData:
%   > nTrials:
%   > nBoots:
%   > binCountArray: 
%   > subsampleArray:
% Outputs
%   > I_est: nTrials x nBinSizes x nSubsamplesets x nBoots array of naive
%   estimates. To be extrapolated to infinitte data and 0 bin size.
 
[nX,nY] = size(yData);

assert(nX == numel(xData));

m = round(subsampleArray*nY);
nSubsamplesets = numel(m);

replacement = false;

yIdx = [1:1:nY];

nBinSizes = numel(binCountArray);

pxBin = 1/nX; % Assume uniform distribution
px = repmat(pxBin,nX,1); % independent variable distribution

I_est = zeros(nTrials, nBinSizes, nSubsamplesets, nBoots);

myCluster = parcluster('local');
p = gcp('nocreate');
if isempty(p)
    parpool(myCluster.NumWorkers);
else
    % Parpool already running.
end

for iTrial = 1:nTrials
    disp(['Beginning trial ',num2str(iTrial),'/',num2str(nTrials)])
    tic
    iBinSizeIdx = 0;
    for iBinSize = 1:nBinSizes
        nBinsY = binCountArray(iBinSize);
        iBinSizeIdx = iBinSizeIdx + 1;
        for iSubsampleIdx = 1:nSubsamplesets
            % Initialize distributions for subsample size.
            mSubsamps = m(iSubsampleIdx);
            py = zeros(1,nBinsY); % signal distribution         
            pyxJ = zeros(nX, nBinsY); % joint distribution
            for iBoot = 1:nBoots
                k = m(iSubsampleIdx);
                subIdx = randsample(yIdx,k,replacement);
                sub = yData(:,subIdx);
                
                py = histcounts(sub, ... 
                    nBinsY, ...
                    'Normalization', 'probability');
                
                pyxJ = histcounts2(repmat(xData',1,mSubsamps),sub, ...
                    [nX, nBinsY], ...
                    'Normalization', 'probability');
                
                I_estTemp = sum(nansum( pyxJ.*log2(pyxJ./(px.*py)) ));
                
                I_est(iTrial, iBinSize, iSubsampleIdx, iBoot) = I_estTemp;
                
            end
        end
    end
    toc
    disp([num2str(round(100*iTrial/nTrials,2)),'%'])
end

save(['direct-mi-estimate_'],'I_est','nBoots','binCountArray','subsampleArray')