function MI_est = ...
    midirectestimate(yData,nTrials,nBoots,binCounts,subSamps)
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

[nX,nE,nG] = size(yData);
x = 1/nX:1/nX:1;

assert(nX == numel(x));

% nX = 100; % For pSmad data with such small sample size.

m = round(subSamps*nE);
nSubsamplesets = numel(m);

replacement = false;

yIdx = [1:1:nE];

nBinSizes = numel(binCounts);

pxBin = 1/nX; % Assume uniform distribution
px = repmat(pxBin,nX,1); % independent variable distribution


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
        nBinsG = binCounts(iBinSize);
        
        nBinsX = nBinsG;
        pxBin = 1/nBinsX; % Assume uniform distribution
        px = repmat(pxBin,nBinsX,1); % independent variable distribution

        iBinSizeIdx = iBinSizeIdx + 1;
        for iSubsampleIdx = 1:nSubsamplesets
            % Initialize distributions for subsample size.
            mSubsamps = m(iSubsampleIdx);
            pg = zeros(1,nBinsG); % signal distribution         
            pgxJ = zeros(nBinsX, nBinsG); % joint distribution
            MI_estTemp = zeros(nBoots,1);
            for iBoot = 1:nBoots
                k = m(iSubsampleIdx);
                subIdx = randsample(yIdx,k,replacement);
                sub = yData(:,subIdx);

                sub = reshape(sub,numel(sub),1);

                xMat = repmat(x',1,mSubsamps);
                xColumn = reshape(xMat,numel(xMat),1);
                
                %% Joint distribution of signal(s) and position.
                pgxJ = histcountsn([xColumn, sub],[nBinsX,nBinsG], ...
                    'Normalization','probability');
                pgxJ = pgxJ + eps;  
                
                pg = squeeze(sum(pgxJ,1));
                
                %% Mutual information calculation using Eq.7 from Tkacik et al 2015
                %%% SATYA
                pgxJ_satya = histcountsn([xColumn, sub],[nBinsX, nBinsG], ...
                    'Normalization','probability'); % P(x, g1, g2): The joint PDF for (x, g1, g2)
                pgxJ_satya= pgxJ_satya+eps; % Add eps so that 0*log2(0) is replaced by eps*log2(eps). Also, matlab has its own eps.

                % P(x) = PDF of x
                px_satya= squeeze(sum(pgxJ_satya, [2 3]));
                % P(g1, g2) = (joint) PDF of {g1,g2}
                pg_satya= squeeze(sum(pgxJ_satya, 1)); 
                % P(g1,g2 | x)= Joint PDF estimate of {g1,g2} conditional on x 
                pg_given_x_satya= nan(size(pgxJ_satya)); 
                
                MIest_satya = 0; 
                for iX = 1:nBinsX
                    
                    % P(g1,g2 | x) for x= xData(xVar)
                    temp_conditional_pdf = ...
                        squeeze(pgxJ_satya(iX, :, :)) / px_satya(iX); 
                    
                    pg_given_x_satya(iX, :, :) = temp_conditional_pdf;

                    temp_logValue = ...
                        temp_conditional_pdf .* log2( temp_conditional_pdf ./ pg_satya); % The inner summation/intergraion term in Eq. 7 
                    MIest_satya = ...
                        MIest_satya + px_satya(iX) * sum(temp_logValue(:));
                end %%% WORKING
                
                %%% MINE
                px = squeeze(sum(pgxJ, [2 3]));
                assert(sum(px - px_satya) == 0);
                
                pg_given_x = pgxJ ./ px;
                assert(sum(pg_given_x - pg_given_x_satya,'all') == 0);
                assert(sum(pg - pg_satya) == 0);
                
                logValue = pg_given_x .* log2(pg_given_x ./ pg);
                assert(sum(logValue(end,:) - temp_logValue) < 1e-10);
                
                MIest = sum(px .* logValue,'all');
                assert(sum(MIest_satya - MIest) < 1e-10)
                
                %% Mutual information calculation using Eq.24 from Tkacik et al 2015

                MI_estTemp(iBoot) = ...
                    nansum( pgxJ.*log2(pgxJ./(px.*pg)), 'all');
                
                assert(sum(MI_estTemp(iBoot) - MIest) < 1e-10)
                                  
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

save(['direct-mi-estimate'],'MI_est','nBoots','binCounts','subSamps','m','nG')
