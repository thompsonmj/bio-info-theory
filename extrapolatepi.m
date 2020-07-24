%% Note
% Minimally working version for estimating results from midirectestimate.m
% and midirectestimate2.m

%%% Still needs error bars

% MI_est.trialData(nTrials).naiveEst_means annd .naiveEst_stds

nBinSizes = numel(MI_est(1).binSizes);
nSubsamps = numel(MI_est(1).subSamps);
nTrials = numel(MI_est);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expData = d.Hb;
% expData = yData;
% expData = Y_aligned;
% nY = numel(expData);
m = round(MI_est(1).subSamps*nEmbryos);

invSamps = 1./m;
invSamps = invSamps*3;
invBinsSq = 1./(MI_est(1).binSizes.^2);
invBinsSq = transpose(invBinsSq);
xRegress = [0:0.001:max(invSamps)];
intercepts_infiniteData = zeros(nBinSizes,nTrials);
intercepts_infiniteBins = zeros(1,nTrials);

themean = zeros(nBinSizes,nSubsamps);

% For each trial, extrapolate to infinite data.
for iT = 1:nTrials
    for iB = 1:nBinSizes
        coeffs_thisB = ...
            polyfit(invSamps, MI_est(iT).naiveEst_means(iB,:),1);
%         f = polyval(coeffs_thisB,xRegress);
%         thisIntercept = polyval(coeffs_thisB,0);
        thisIntercept = coeffs_thisB(2);
        intercepts_infiniteData(iB,iT) = thisIntercept;
    end
    %%% Now: use data from each trial in 'intercepts_infiniteData' to
    %%% extrapolate to infinite bins to get stats on this point. plot the curve
    %%% fit using the mean infinite sample intercepts
    coeffs_thisT = ...
        polyfit(invBinsSq, intercepts_infiniteData(:,iT),1);
    thisIntercept = coeffs_thisT(2);
    intercepts_infiniteBins(iT) = thisIntercept;
end

intercepts_infiniteData_means = mean(intercepts_infiniteData,2);
intercepts_infiniteData_stds = std(intercepts_infiniteData,0,2);

% Get stats across mean value of all boots for each trial.
%   Populate a new matrix to use the data in the struct array.
naiveEst_trialMeans = mean(cat(3,MI_est(:).naiveEst_means),3);
naiveEst_trialStds = std(cat(3,MI_est(:).naiveEst_means),0,3); %%% want error propagation

%% Plot results
% Plot 1: MI vs. 1/(no. samples)
% % % % % % figure
% % % % % % hold on
for iB = 1:nBinSizes
    % Get the data from each trial into one spot
    MI = zeros(nBinSizes,nSubsamps,nTrials);
    for iT = 1:nTrials
        MI(:,:,iT) = MI_est(iT).naiveEst_means;
    end
    % Take the stats across trials
    MI_means = mean(MI,3);
    MI_stds = std(MI,0,3);
    % Plot the matrix
% % % % % %     for iB = 1:nBinSizes
% % % % % %         scatter(1./round(MI_est(1).subSamps*nEmbryos),MI_means(iB,:))
% % % % % %     end
end
% % % % % % xlim([0,0.1])
% % % % % % ylim([1.5,2.5])
% % % % % % xticks([0,0.05,0.1])
% % % % % % yticks([1.5,2,2.5])

%% Working to here
%   
%%
% Plot 2: MI vs. 1/(no. bins)^2
% % % % % % figure
% % % % % % hold on

% plot(repmat(invSamps,nBinSizes,1),IMat,'o','k')
for iB = 1: nBinSizes
    coeffs_thisB = polyfit(invSamps,MI_means(iB,:),1);
    f = polyval(coeffs_thisB,xRegress);
    intercept = polyval(coeffs_thisB,0);
    intercepts_infiniteData(iB) = intercept;
%     scatter(invSamps,meanI(iB,:),'o','k')
% % % % % %     e = errorbar(invSamps,MI_means(iB,:),MI_stds(iB,:));
% % % % % %     e.Marker = 'o';
%     e.MarkerColor = 'k';
% % % % % %     e.Color = 'k';
% % % % % %     plot(xRegress,f,'-k')
% % % % % %     scatter(0,intercept,'o','b')
%     for iT = 1:nTrials
%         scatter(invSamps,MI_est(iB,:,iT),'k','.')
%     end
end


meanIntercepts = mean(intercepts_infiniteData,2);
stdIntercepts = std(intercepts_infiniteData,1,2);

% xlim([-0.001,0.1])
% ylim([1.5,2.5])
% % % % % % xlabel('1/(# samples)')
% % % % % % ylabel('I [bits]')

invBinsSq = 1./(binCounts.^2);

p2 = polyfit(invBinsSq',meanIntercepts,1);

f2 = polyval(p2,[0:0.001:max(invBinsSq)]);

intercept2 = polyval(p2,0);

MI_ext = intercept2;
disp(MI_ext)

% % % figure
% % % % % % hold on
% scatter(invBins,meanIntercepts,'b')
% % % % % % scatter(invBinsSq,meanIntercepts,'r')
% e = errorbar(invBins,meanIntercepts,stdIntercepts);

% plot([0:0.001:max(invBins)],f,'-b')
% % % % % % plot([0:0.001:max(invBinsSq)],f2,'-r')

% scatter(0,intercept,'sq','b')
% % % % % % scatter(0,intercept2,'sq','r')
% % % % % % xlim([-0.1*10^-3,3*10^-3])
% % % % % % ylim([1,3])
% % % % % % xticks([0,1,2,3]*1e-3)
% % % % % % yticks([1,1.5,2,2.5,3])
% % % % % % xlabel('1/(# bins)^2')
% % % % % % ylabel('I [bits]')

% cd(fullfile('..','norm_sandbox'))