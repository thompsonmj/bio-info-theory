function MI_est = ...
    midirectestimate2(yData,nTrials,nBoots,binCounts,subSamps)
%MIDIRECTESTIMATE computes a direct estimate of mutual information through
%iterative naive estimates.
% 
% Inputs
%   > yData:
%   > xData:
%   > nTrials:
%   > nBoots:
%   > binCountArray: 
%   > subsampleArray:
% Outputs
%   > MI_est: nTrials x nBinSizes x nSubsamplesets x nBoots array of naive
%   estimates. To be extrapolated to infinitte data and 0 bin size.
%
%%% Working on implementation for >1 gene 
 
[nX,nY,nG] = size(yData);

%% Trim away NaNs used for padding and x alignment if present
padSize = unique(sum(isnan(yData(:,:,1))))/2;
assert(numel(padSize) == 1);

yData = yData(padSize:nX-padSize+1,:,[1,4]);
[nX,nY,nG] = size(yData);
x = 1/nX:1/nX:1;

assert(nX == numel(x));

% nX = 100; % For pSmad data with such small sample size.

m = round(subSamps*nY);
nSubsamplesets = numel(m);

replacement = false;

yIdx = [1:1:nY];

nBinSizes = numel(binCounts);

pxBin = 1/nX; % Assume uniform distribution
% px = repmat(pxBin,nX,1); % independent variable distribution


for iTrial = 1:nTrials
    % Initialize for trial.
    MI_est(iTrial).binSizes = binCounts;
    MI_est(iTrial).subSamps = subSamps;
    MI_est(iTrial).nBoots = nBoots;
    
    MI_est(iTrial).naiveEst_means = zeros(nBinSizes,nSubsamplesets);
    MI_est(iTrial).naiveEst_stds = zeros(nBinSizes,nSubsamplesets);
    if nTrials > 1
        disp(['Beginning trial ',num2str(iTrial),'/',num2str(nTrials)])
        tic
    end
    iBinSizeIdx = 0;
    for iBinSize = 1:nBinSizes
        if nTrials == 1
            disp(['Beginning bin ',num2str(iBinSize),'/',num2str(nBinSizes)])
            tic
        end
        nBinsY = binCounts(iBinSize);
        
        nX = nBinsY;
        
%         nBinsY = round(sqrt(nX));
        pxBin = 1/nX; % Assume uniform distribution
        px = repmat(pxBin,nX,1); % independent variable distribution
        
        iBinSizeIdx = iBinSizeIdx + 1;
        for iSubsampleIdx = 1:nSubsamplesets
            % Initialize distributions for subsample size.
            mSubsamps = m(iSubsampleIdx);
            py = zeros(nBinsY,nBinsY); % signal distribution         
            pyxJ = zeros(nX, nBinsY,nBinsY); % joint distribution
            MI_estTemp = zeros(nBoots,1);
            for iBoot = 1:nBoots
                k = m(iSubsampleIdx);
                subIdx = randsample(yIdx,k,replacement);
                sub1 = yData(:,subIdx,1);
                sub2 = yData(:,subIdx,2);

                sub1 = reshape(sub1,numel(sub1),1);
                sub2 = reshape(sub2,numel(sub2),1);

                xMat = repmat(x',1,mSubsamps);
                xColumn = reshape(xMat,numel(xMat),1);
                
                %% 1
                py = histcountsn([sub1,sub2], ... 
                    [nBinsY,nBinsY]);
                
%                 py = histcounts2(
                
%                 pyxJ1 = histcounts2(xColumn,sub,[nX,nBinsY],'Normalization','probability');
                
                %% 2
                pyxJ2 = histcountsn([xColumn,sub1,sub2],[nX,nBinsY,nBinsY]);

                pyxJ = pyxJ2;
                px = repmat(pxBin,nX,nBinsY); % independent variable distribution

                MI_estTemp(iBoot) = sum(sum(nansum( pyxJ.*log2(pyxJ./(px*py)) )));
                                  
            end
            
            MI_est(iTrial).naiveEst_means(iBinSize,iSubsampleIdx) = ...
                mean(MI_estTemp);
            MI_est(iTrial).naiveEst_stds(iBinSize,iSubsampleIdx) = ...
                std(MI_estTemp,0);
            
        end
        
        if nTrials == 1
            toc
            disp([num2str(round(100*iBinSize/nBinSizes,2)),'%'])
        end
    end
    
    if nTrials > 1
        toc
        disp([num2str(round(100*iTrial/nTrials,2)),'%'])
    end
end

save(['direct-mi-estimate'],'MI_est','nBoots','binCounts','subSamps','m','nY')
