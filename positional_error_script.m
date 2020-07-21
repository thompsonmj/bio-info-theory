
load('gap_data_raw_dorsal_wt_XY-Aligned.mat')

%% Positional error
% [nSamplePts,~,~] = size(Y_align);
% sigma_x = zeros(nSamplePts-1,nGenes);
% figure
% hold on
% for iG = 1:nGenes
%     sigma_x(:,iG) = positionalerror(Y_align(:,:,iG));
%     
%     subplot(2,2,iG)
%     
%     yyaxis left
%     plot(Y_align(:,:,iG),'-')
%     ylabel('Relative Fluorescence')
%     
%     yyaxis right
%     plot(sigma_x(:,iG))
%     ylim([0,0.1])
%     ylabel('Positional Error (\sigma_x/L)')
%     title(genes{iG})
% end

%% PE for multiple genes simultaneously
sigma_x = positionalerrorn(Y_align);

