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
% Based on Eq. 7 and Eq. 24 in Tkacik et al, 2015.
 
[nBinsX,nE,nG] = size(yData);
x = 1/nBinsX:1/nBinsX:1;

assert(nBinsX == numel(x));

% nX = 100; % For pSmad data with such small sample size.

m = round(subSamps*nE);
nSubsamplesets = numel(m);

replacement = false;

embryoIdcs = [1:1:nE];

nBinSizes = numel(binCounts);

for iTrial = 1:nTrials
    % Initialize this trial.
    MI_est(iTrial).binSizes = binCounts;
    MI_est(iTrial).subSamps = subSamps;
    MI_est(iTrial).nBoots = nBoots;
    
    MI_est(iTrial).naiveEst_means = zeros(nBinSizes,nSubsamplesets);
    MI_est(iTrial).naiveEst_stds = zeros(nBinSizes,nSubsamplesets);
    if nTrials > 1
%         disp(['Beginning trial ',num2str(iTrial),'/',num2str(nTrials)])
%         tic
    end
    iBinSizeIdx = 0;
    for iBinSize = 1:nBinSizes
        if nTrials == 1
%             disp(['Beginning bin ',num2str(iBinSize),'/',num2str(nBinSizes)])
%             tic
        end
        nBinsG = binCounts(iBinSize);
        
        nBinsX = nBinsG;
        pxBin = 1/nBinsX; % Assume uniform distribution
        
        iBinSizeIdx = iBinSizeIdx + 1;
        for iSubsampleIdx = 1:nSubsamplesets
            % Initialize distributions for subsample size.
            mSubsamps = m(iSubsampleIdx);
            pg = zeros(nBinsG,nBinsG); % signal distribution         
            pgxJ = zeros(nBinsX, nBinsG,nBinsG); % joint distribution
            MIEq24 = zeros(nBoots,1);
            for iBoot = 1:nBoots
                k = m(iSubsampleIdx);
                subIdx = randsample(embryoIdcs,k,replacement);
                sub1 = yData(:,subIdx,1);
                sub2 = yData(:,subIdx,2);

                sub1 = reshape(sub1,numel(sub1),1);
                sub2 = reshape(sub2,numel(sub2),1);

                xMat = repmat(x',1,mSubsamps);
                xColumn = reshape(xMat,numel(xMat),1);
                
                %% Distributions
                pgxJ = histcountsn([xColumn,sub1,sub2],[nBinsX,nBinsG,nBinsG], ...
                    'Normalization', 'probability');
                pgxJ = pgxJ + eps;
                pgxJ_perm = permute(pgxJ,[2,3,1]);

                % P(g1,g2): Joint PDF of {g1,g2}
                pg = squeeze(sum(pgxJ,1));

                % P(x): PDF of x (uniform)
                px = repmat(pxBin,nBinsX,1);                                  
                
                % P({g1,g2} | x): Joint PDF estimate of {g1,g2} conditional on x 
                pg_given_x = pgxJ ./ px;
                pg_given_x_perm = permute(pg_given_x,[2,3,1]);
                
                %% Eq. 7
                % Inner summation/intergraion term
                logValue = pg_given_x_perm .* log2(pg_given_x_perm ./ pg);
                
                MIEq7_temp = sum(px .* sum(logValue,3),'all');
                
                MIEq7(iBoot) = MIEq7_temp;
                 
                %% Eq.24
                MIEq24_temp = ...
                    sum( pgxJ_perm .* log2(pgxJ_perm ./ (px .* pg)), 'all' );  
                
                MIEq24(iBoot) = MIEq24_temp;
                
                %% Validation
                assert(abs(MIEq7_temp - MIEq24_temp) < 1e-10)
                             
            end
                        
            MI_est(iTrial).naiveEst_means(iBinSize,iSubsampleIdx) = ...
                mean(MIEq24);
            MI_est(iTrial).naiveEst_stds(iBinSize,iSubsampleIdx) = ...
                std(MIEq24,0);
            
        end
        
        if nTrials == 1
%             toc
%             disp([num2str(round(100*iBinSize/nBinSizes,2)),'%'])
        end
    end
    
    if nTrials > 1
%         toc
%         disp([num2str(round(100*iTrial/nTrials,2)),'%'])
    end
end

%save(['direct-mi-estimate'],'MI_est','nBoots','binCounts','subSamps','m','nE')
