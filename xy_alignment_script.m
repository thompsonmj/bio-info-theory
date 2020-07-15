% % % clear
% % % 
% % % %%
% % % nEmbryosUsed = 102 % 102 fly embryos total
% % % inputData = 'gap_data_raw_dorsal_wt'; % variable stored: 'data'
% % % inputFile = [inputData,'.mat'];
% % % 
% % % %%
% % % load(inputFile);
% % % nEmbryos = numel(data);
% % % assert(nEmbryos >= nEmbryosUsed);
% % % 
% % % % Authors' description of the data: "Each profile vector is of length 1000
% % % % pixels, where 0 corresponds to the anterior (A) and 1000 to the posterior
% % % % (P) of the embryo."
% % % % To mitigate 'edge effects,' we extract only the middle 80% of the data to 
% % % % use throughout the analysis.
% % % 
% % % nSamplePts = numel(data(1).Hb);
% % % x = (1/nSamplePts : 1/nSamplePts : 1);
% % % idx_low = round(nSamplePts)*0.1;
% % % idx_high = round(nSamplePts)*0.9;
% % % x = x(idx_low:idx_high);
% % % 
% % % genes = {'Hb','Kr','Gt','Kni'};
% % % nGenes = numel(genes);
% % % 
% % % %% Set up raw data
% % % Y_raw = zeros(nSamplePts,nEmbryosUsed,nGenes);
% % % idcs = sort(randperm(nEmbryos,nEmbryosUsed));
% % % 
% % % for iG = 1:nGenes
% % %     for iE = 1:nEmbryos
% % %         Y_raw(:,iE,iG) = data(iE).(genes{iG});
% % %     end
% % % end
% % % Y_raw = Y_raw(idx_low:idx_high,idcs,:);
% % % [nPts,~,~] = size(Y_raw);
% % % padSize = round(nSamplePts*0.05);
% % % Y_raw = padarray(Y_raw,padSize,NaN,'both');

%% Cochlea data
load('F:\projects\cochlea\data\cochlea-data_2020-07-15_11-36-59.mat')
[fsa,ft] = filtercochleadata(data,12.5,true,{'bmp4_mrna','topro3'});
clear Y_raw
clear Y_yAlign1
clear Y_xAlign
clear Y_align
genes = {'bmp4_mrna'};
nGenes = numel(genes);

gene1 = ft.bmp4_mrna;
gene2 = ft.topro3;

[nSamplePts,nC,~] = size(gene1);

for iC = 1:nC
    gene1(:,iC) = regressnoise(gene2(:,iC), gene1(:,iC));
%     gene2(:,iC) = regressnoise(gene3(:,iC), gene2(:,iC));
    
    Y_raw(:,iC,1) = smoothrawdata(gene1(:,iC),'movmean',50);
%     Y_raw(:,iC,2) = smoothrawdata(gene2(:,iC),'movmean',50);
end

% Y_raw(:,2,:) = [];
% Y_raw(:,5,:) = [];
% Y_raw(:,6,:) = [];
[nSamplePts,nC,~] = size(Y_raw(:,:,1));

% Remove first and last 10%.
nSamplePts = numel(Y_raw(:,1,1));
x = (1/nSamplePts : 1/nSamplePts : 1);
idx_low = round(nSamplePts)*0.1;
idx_high = round(nSamplePts)*0.9;
x = x(idx_low:idx_high);

% Pad ends with NaNs.
Y_raw = Y_raw(idx_low:idx_high,:,:);
[nPts,~,~] = size(Y_raw);
padSize = round(nSamplePts*0.1);
Y_raw = padarray(Y_raw,padSize,NaN,'both');

%% Comments
% Up to this point, we have simply deleted data outside of the middle 80%
% of the domain (replacing it with NaNs to act as a buffer to shift the
% profiles to the left/right during X-alignment).
%% 1 Align Y alone (optimize alpha and beta).
[Y_raw,~] = anchormean0to1(Y_raw);
Y_yAlign1 = zeros(size(Y_raw));
for iG = 1:nGenes
    y = squeeze(Y_raw(:,:,iG));
    [Y_yAlign1(:,:,iG),~] = aligny(y);
end

figure
hold on
plot(Y_raw(:,:,1))
% plot(Y_raw(:,:,2))

%% 2 Numerically estimate optimal joint XY alignment (alpha, beta, gamma).
% UseParallel_TF = false;
% pool = gcp;
% if UseParallel_TF
%     if isempty(gcp('nocreate'))
%         parpool(pool.NumWorkers)
%     end
% end
[Y_xAlign,alpha,beta,gamma] = alignxy(Y_yAlign1,false);
% Should explore a 'better' solution than genetic optimization. This was
% chosen because I shift profiles in x using their integer indices, and
% MATLAB's genetic optimizer allowed for easily implemented integer
% constraints on one parameter (gamma) and continuous variation on others
% (alpha and beta). This approach does not yield the same result twice, but
% reliably yields a better x-alignment than the data came with originally.
% However, it needs tweaking for the data at hand (e.g. the NaN buffer size
% on either side of each profile should be large enough, but should not be
% too large). Additionally, y-alignment is poor when done jointly using
% this method, so an additional chi^2 minimization must be run afterward.

%% Trim away NaNs used for padding and x alignment.
[nX,~,~] = size(Y_xAlign);
Y_xAlign = Y_xAlign(padSize:nX-padSize+1,:,:);

%% 3 Align Y alone once more for closed-form Y alignment.
Y_align = zeros(size(Y_xAlign));
for iG = 1:nGenes
    y = squeeze(Y_xAlign(:,:,iG));
    [Y_align(:,:,iG),~] = aligny(y);
end

figure
hold on
plot(Y_align(:,:,1),'k')
% plot(Y_align(:,:,2),'r')

figure(2)
hold on
plot(mean(Y_align(:,:,2),2),'r','linew',2)

% figure
% hold on
% for iC = 1:nC
%     plot(Y_align(:,iC,1),'b')
%     plot(Y_align(:,iC,2),'r')
%     title(iC)
%     pause(4)
% end

% % % %% Save final result
% % % save([inputData,'_XY-Aligned.mat'])