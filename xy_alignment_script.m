clear
owd = pwd;
%%
nEmbryosUsed = 20
nEmbryos = 20;

% inputData = 'gap_data_raw_dorsal_wt';
inputData = 'cochlea-data_2020-03-03_14-24-55';

inputFile = [inputData,'.mat']; % variable stored: 'data'


%%
load(inputFile);
totalEmbryos = numel(data);
count = 0;
smooth_TF = true;

for iW = 1:10
    clearvars -except iW smooth_TF nEmbryos inputFile inputData nEmbryosUsed owd data totalEmbryos
    clear d
    count = 0;
    for iE = 1:totalEmbryos
        if data(iE).age == 12.5 && ...
                data(iE).flag && ...
                ~isempty(data(iE).psmad) && ...
                ~isempty(data(iE).topro3)% && ...
    %             ~isempty(data(iE).sox2)

            count = count + 1;

            topro3 = uniformsample(data(iE).topro3,1000);
            psmad = uniformsample(data(iE).psmad,1000);
    %         sox2 = data(iE).sox2;

            if smooth_TF
                win = iW*80;
                method = 'lowess';
                d(count).topro3 = smoothrawdata(topro3,method,win);
                d(count).psmad = smoothrawdata(regressnoise(topro3,psmad),method,win);
%                 d(count).sox2 = smoothrawdata(regressnoise(topro3,sox2),method,win);  
            else
                d(count).topro3 = topro3;
                d(count).psmad = regressnoise(topro3,psmad);
            end

        end
    end
    %%
    assert(nEmbryos >= nEmbryosUsed);
    % nEmbryos = numel(data);
    nEmbryos = numel(d)

    % nSamplePts = numel(data(1).Hb);
    nSamplePts = numel(d(1).topro3);
    x = (1/nSamplePts : 1/nSamplePts : 1);
    idx_low = round(nSamplePts)*0.1;
    idx_high = round(nSamplePts)*0.9;
    x = x(idx_low:idx_high);

    % genes = {'Hb','Kr','Gt','Kni'};
    genes = {'psmad'};
    nGenes = numel(genes);

    %% Set up raw data
    Y_raw = zeros(nSamplePts,nEmbryosUsed,nGenes);
    idcs = sort(randperm(nEmbryos,nEmbryosUsed));

    for iG = 1:nGenes
        for iE = 1:nEmbryos
    %         Y_raw(:,iE,iG) = data(iE).(genes{iG});
            Y_raw(:,iE,iG) = d(iE).(genes{iG});
        end
    end
    Y_raw = Y_raw(idx_low:idx_high,idcs,:);
    [nPts,~,~] = size(Y_raw);
    padFraction = 0.05;
    padSize = round(nPts*padFraction);
    Y_raw_pad = padarray(Y_raw,padSize,NaN,'both');

    %% 1 Align Y alone (optimize alpha and beta).
    [Y_raw_pad,~] = anchormean0to1(Y_raw_pad);
    Y_yAlign1_pad = zeros(size(Y_raw_pad));
    for iG = 1:nGenes
        y = squeeze(Y_raw_pad(:,:,iG));
        [Y_yAlign1_pad(:,:,iG),~] = aligny(y);
    end

    %% 2 Numerically estimate optimal joint XY alignment (alpha, beta, gamma).
    UseParallel_TF = true;
    pool = gcp;
    if UseParallel_TF
        if isempty(gcp('nocreate'))
            parpool(pool.NumWorkers)
        end
    end
    [Y_xyAlign_pad,alpha,beta,gamma] = alignxy(Y_yAlign1_pad,UseParallel_TF);

    %% 3 Align Y alone once more for closed-form Y alignment.
    Y_yAlign2_pad = zeros(size(Y_yAlign1_pad));
    for iG = 1:nGenes
        y = squeeze(Y_xyAlign_pad(:,:,iG));
        [Y_yAlign2_pad(:,:,iG),~] = aligny(y);
    end

    Y_aligned_pad = Y_yAlign2_pad;

    [nX0,~] = size(Y_aligned_pad);
    idxLow = padSize + 1;
    idxHigh = nX0 - padSize;
    Y_yAlign1 = Y_yAlign1_pad(idxLow:idxHigh,:,:);
    Y_aligned = anchormean0to1(Y_aligned_pad(idxLow:idxHigh,:,:));
    gene1 = Y_aligned(:,:,1);
    [nX,~] = size(Y_aligned);
    x = 0.1+((0.9-0.1)/nX : (0.9-0.1)/nX : (0.9-0.1));
    X = repmat(x,nEmbryosUsed,1)';

    %% Visualize results
    fh = figure;
    for iG = 1:nGenes
        hold on
        plot(X,Y_aligned(:,:,iG),'color',[0.25,0.25,0.25],'linew',2)
        plot(x,nanmean(Y_aligned(:,:,iG),2),'color','k','linew',4)

        xlim([0,1])
        ylim([-0.1,1.2])
        xticks([0,0.2,0.4,0.6,0.8,1.0])
        yticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
        grid on
        title(genes{iG})
    end

    if smooth_TF
        fname = [inputData,'_',method,'-sm_',num2str(iW/10),'-win_',num2str(padFraction),'-pad_','XY-Aligned'];
    else
        fname = [inputData,'_raw_',num2str(padFraction),'-pad_','XY-Aligned'];
    end
    

%     cd(fullfile('..','mi_sandbox'))
    
%     mi_test_script

%     extrapolatepi
    
%     PI_bits = intercept2;
    
%     title([genes{iG},' ',num2str(PI_bits),' bits'])
    
    cd(fullfile('..','norm_sandbox',method))
    savefig(fh,[fname,'.fig'])
    save([fname,'.mat'])
    cd('..')
    
    close(fh)
    
end