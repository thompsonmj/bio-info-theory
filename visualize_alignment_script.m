clear

load('gap_data_raw_dorsal_wt_XY-Aligned.mat')

%% Visualize results
[n_x_raw,~,~] = size(Y_raw);
[n_x_align,~,~] = size(Y_align);
xStart = 0.1;
xEnd = 0.9;
x_raw = xStart + (xEnd - xStart)*(1/n_x_raw : 1/n_x_raw : 1);
x_align = xStart + (xEnd - xStart)*(1/n_x_align : 1/n_x_align : 1);

figure
for iG = 1:nGenes
    subplot(2,2,iG)
    hold on
    plot(x_raw,Y_raw(:,:,iG),'color',[0.5,0.5,0.5],'lines',':')
    plot(x_align,Y_align(:,:,iG),'color','k')
    
    plot(x_raw,squeeze(nanmean(Y_raw(:,:,iG),2)),'color',[0.5,0.5,0.5],'linew',3)
    plot(x_align,squeeze(nanmean(Y_align(:,:,iG),2)),'color','k','linew',3)
    grid on
    xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    yticks([0,0.2,0.4,0.6,0.8,1])
    xlim([-0.1 1.1])
    ylim([-0.1 1.1])
    sgtitle('Raw vs. Fully Aligned')
end
