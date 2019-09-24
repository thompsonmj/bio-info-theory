% clear

% load('at100percent.mat')
% load('parat100percent.mat')

% I_est: nTrials x nBinSizes x nSubsamplesets x nBoots


% I_est1: nBinSizes x nSubsamplesets x nBoots
I_estMEAN1 = squeeze(mean(I_est,1));
I_estSTD1 = squeeze(std(I_est,1,1));

% I_est2: nBinSizes x nSubsamplesets
I_estMEAN2 = squeeze(mean(I_estMEAN1,3));
I_estSTD2 = squeeze(std(I_estMEAN1,1,3));

pCell = cell(nBinSizes,1);

invSamps = 1./m;
invBins2 = 1./(binCounts.^2);
invBins2 = transpose(invBins2);
xRegress = [0:0.001:max(invSamps)];
intercepts = zeros(nBinSizes,1);


figure
hold on
% plot(repmat(invSamps,nBinSizes,1),IMat,'o','k')
for iB = 1: nBinSizes
    pTemp = polyfit(invSamps,I_estMEAN2(iB,:),1);
    f = polyval(pTemp,xRegress);
    intercept = polyval(pTemp,0);
    intercepts(iB) = intercept;
%     scatter(invSamps,meanI(iB,:),'o','k')
    e = errorbar(invSamps,I_estMEAN2(iB,:),I_estSTD2(iB,:));
    e.Marker = 'o';
%     e.MarkerColor = 'k';
    e.Color = 'k';
    plot(xRegress,f,'-k')
    scatter(0,intercept,'o','b')
%     for iT = 1:nTrials
%         scatter(invSamps,I(iB,:,iT),'k','.')
%     end
end


meanIntercepts = mean(intercepts,2);
stdIntercepts = std(intercepts,1,2);

% xlim([-0.001,0.1])
% ylim([1.5,2.5])
xlabel('1/(# samples)')
ylabel('I [bits]')

invBins = 1./(binCounts);
invBins2 = 1./(binCounts.^2);

p = polyfit(invBins',meanIntercepts,1);
p2 = polyfit(invBins2',meanIntercepts,1);

f = polyval(p,[0:0.001:max(invBins)]);
f2 = polyval(p2,[0:0.001:max(invBins2)]);

intercept = polyval(p,0);
intercept2 = polyval(p2,0);

figure
hold on
scatter(invBins,meanIntercepts,'b')
scatter(invBins2,meanIntercepts,'r')
% e = errorbar(invBins,meanIntercepts,stdIntercepts);

plot([0:0.001:max(invBins)],f,'-b')
plot([0:0.001:max(invBins2)],f2,'-r')

scatter(0,intercept,'sq','b')
scatter(0,intercept2,'sq','r')
xlim([-0.1*10^-3,3*10^-3])
ylim([1,3])
xlabel('1/(# bins)^2')
ylabel('I [bits]')