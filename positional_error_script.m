clear
load('gap_data_raw_dorsal_wt_XY-Aligned.mat','Y_align')
genes = {'Hb','Kr','Gt','Kni'};
[nX,nE,nG] = size(Y_align);
idx = nX*0.1;
Y_align = Y_align(idx:end-idx,:,:);
[nX,nE,nG] = size(Y_align);

x = 1/nX:1/nX:1;

%% Positional error for single genes.
figure
title('Positional Error for Single Genes')
for iG = 1:nG
    sigma_x.(genes{iG}) = positionalerror(Y_align(:,:,iG));
    subplot(3,2,iG)
    hold on
    yyaxis left
    plot(x,Y_align(:,:,iG),'-','color',[0.5,0.5,0.5])
    ylabel('Relative Expression Level')
    yyaxis right
    plot(x(1:end-1),sigma_x.(genes{iG}),'linew',2)
    ylabel('Positional Error (\sigma_x/L)')
    ylim([0, 0.1])
    title(genes{iG})
end

%% PE for multiple genes simultaneously
sx = positionalerrorn(Y_align)/nX;
subplot(3,2,[5,6])
hold on
for iG = 1:nG
    yyaxis left
    plot(x,Y_align(:,:,iG),'-')
end
yyaxis right
plot(x,sx,'linew',2)
ylabel('Positional Error (\sigma_x/L)')
ylim([0, 0.1])


