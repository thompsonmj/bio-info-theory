%% 0. 1D
%%% WORKING
clear
nBins = 20;
X = randn(1000,1);

figure
subplot(1,2,1)
histogram(X,nBins,'Normalization','probability')
title('histogram')

[N,XEdges] = histcountsn(X,nBins);
subplot(1,2,2)
histogram('BinEdges',XEdges{1},'BinCounts',N)
title('histcountsn')
%%% WORKING
%% 1. Plot results from histcounts2 to match histogram2.
%%% WORKING
clear
nBins = 20;
X = randn(1000,2);

figure
subplot(1,2,1)
histogram2(X(:,1),X(:,2),nBins)
title('histogram2')

[N,Xedges,Yedges] = histcounts2(X(:,1),X(:,2),nBins);
subplot(1,2,2)
histogram2('BinCounts',N,'XBinEdges',Xedges,'YBinEdges',Yedges)
title('histcounts2')
%%% WORKING
%% 2. Get histcounts2 to provide the same results as histcn with 2D input.

clear
nBins = 20;
X = randn(1000,2);

[N,Xedges,Yedges] = histcounts2(X(:,1),X(:,2),nBins);
figure
subplot(1,2,1)
histogram2('BinCounts',N,'XBinEdges',Xedges,'YBinEdges',Yedges)
title('histogram2')

[count edges mid loc] = histcn(X,nBins-1,nBins-1);
[~,nDims] = size(X);
for iD = 1:nDims
    d = mean(diff(edges{iD}));
    edges{iD} = edges{iD} - d/2;
    edges{iD} = [edges{iD},max(edges{iD}) + d];
end

subplot(1,2,2)
histogram2('BinCounts',count,'XBinEdges',edges{1},'YBinEdges',edges{2})
title('histcn')

%%% NOT QUITE WORKING

%% 3. Make histcountsn match histcounts2 by wrapping histcounts.
%%% WORKING
clear
nBins = 10;
X = randn(1000,2);

[N,Xedges,Yedges] = histcounts2(X(:,1),X(:,2),nBins);
figure
subplot(1,2,1)
histogram2('BinCounts',N,'XBinEdges',Xedges,'YBinEdges',Yedges)
title('histcounts2')

[p,edges] = histcountsn(X,[nBins,nBins]);
subplot(1,2,2)
histogram2('BinCounts',p,'XBinEdges',edges{1},'YBinEdges',edges{2});
title('histcountsn')
%%% WORKING
%% 4. Explore results for 3+ dimensions using histcountsn.
clear
nBins = 30;
nDims = 3;
X = randn(5000,nDims);

%%% HISTCOUNTSN
[p,edges] = histcountsn(X,[nBins,nBins,nBins]);

figure
for ii=1:5
    subplot(1,5,ii)
    histogram2('BinCounts',squeeze(p(6*ii,:,:)),'XBinEdges',edges{1},'YBinEdges',edges{2}, ...
        'DisplayStyle','tile', ...
        'ShowEmptyBins','on')
end
figure
for ii=1:5
    subplot(1,5,ii)
    histogram2('BinCounts',squeeze(p(:,6*ii,:)),'XBinEdges',edges{1},'YBinEdges',edges{2}, ...
        'DisplayStyle','tile', ...
        'ShowEmptyBins','on')
end
figure
for ii=1:5
    subplot(1,5,ii)
    histogram2('BinCounts',squeeze(p(:,:,6*ii)),'XBinEdges',edges{1},'YBinEdges',edges{2}, ...
        'DisplayStyle','tile', ...
        'ShowEmptyBins','on')
end
figure
scatter3(X(:,1),X(:,2),X(:,3))