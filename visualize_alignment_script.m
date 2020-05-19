load('gap_data_raw_dorsal_wt_XY-Aligned.mat')

%% Visualize results
[nX,~,~] = size(Y_raw);
xStart = 0.1;
xEnd = 0.9;
x = xStart + (xEnd - xStart)*(1/nX:1/nX:1);
figure
for iG = 1:nGenes
    subplot(2,2,iG)
    hold on
    plot(x,Y_raw(:,:,iG),'color',[0.5,0.5,0.5],'lines',':')
    plot(x,Y_align(:,:,iG),'color','k')
    
    plot(x,squeeze(nanmean(Y_raw(:,:,iG),2)),'color',[0.5,0.5,0.5],'linew',3)
    plot(x,squeeze(nanmean(Y_align(:,:,iG),2)),'color','k','linew',3)
    grid on
    xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    yticks([0,0.2,0.4,0.6,0.8,1])
    xlim([-0.1 1.1])
    ylim([-0.1 1.1])
    sgtitle('Raw vs. Fully Aligned')
end
