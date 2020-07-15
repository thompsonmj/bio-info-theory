function [g,alpha,beta,gamma] = alignxy(G_0,UseParallel_TF)
tic

disp('Beginning joint XY alignment with parallelization set to:')
UseParallel_TF

[~,nE,nG] = size(G_0);

padSize = sum(isnan(G_0(:,1,1)))/2;

alpha_0 = zeros(nG,nE);
beta_0 = ones(nG,nE);
gamma_0 = zeros(1,nE);

nAlpha = numel(alpha_0);
nBeta = numel(beta_0);
nGamma = numel(gamma_0);

options = optimoptions('ga', ...
    'UseParallel', UseParallel_TF, ...
    'UseVectorized', false);

nP = nAlpha + nBeta + nGamma;

%% Specify indices to rereference parameters in cost function.
for iG = 1:nG
    alphaStart = nAlpha*(iG/nG) - nE + 1;
    alphaEnd = alphaStart + nE - 1;
    idcs.alpha{iG} = alphaStart:alphaEnd;
    
    betaStart = nAlpha + nBeta*(iG/nG) - nE + 1;
    betaEnd = betaStart + nE - 1;
    idcs.beta{iG} = betaStart:betaEnd;
end
idcs.gamma = nP - nE + 1: nP;

f = @(P)chisq(P, G_0, idcs); 

IntCon = (nP - nGamma + 1):nP;

lb_alpha = -ones(nAlpha,1);
ub_alpha = ones(nAlpha,1);

lb_beta = zeros(nBeta,1);
ub_beta = 2*ones(nBeta,1);

lb_gamma = repmat(-padSize/2,nE,1);
ub_gamma = repmat(padSize/2,nE,1);

lb = [lb_alpha; lb_beta; lb_gamma];
ub = [ub_alpha; ub_beta; ub_gamma];

[P_final,~,exitFlag,~,~] = ...
    ga(f,nP,[],[],[],[],lb,ub,[],IntCon,options);

exitFlag

%% Extract parameters from P_final.
alpha = zeros(nG,nE);
beta = zeros(nG,nE);
gamma = zeros(1,nE);

for iE = 1:nE
    for iG = 1:nG
        alpha(iG,iE) = P_final(idcs.alpha{iG}(iE));
        beta(iG,iE) = P_final(idcs.beta{iG}(iE));
    end
    gamma(iE) = P_final(idcs.gamma(iE));
end

%% Apply parameters alpha, beta, and gamma to G_0.
G = zeros(size(G_0));
g = zeros(size(G));
for iE = 1:nE
    G(:,iE,:) = circshift( G_0(:,iE,:),gamma(iE),1 );
end
for iG = 1:nG
    alpha_tmp = alpha(iG,:);
    beta_tmp = beta(iG,:);
    g(:,:,iG) = (G(:,:,iG) - alpha_tmp)./beta_tmp;
end

[g,~] = anchormean0to1(g);

shiftAmt = -round((mean(gamma)));
for iE = 1:nE
    g(:,iE,:) = circshift(g(:,iE,:),shiftAmt,1 );
end
toc

end

%% Cost function
function chi2 = chisq(P, G, idcs)
    [~,nE,nG] = size(G);

    alpha = zeros(nG,nE);
    beta = zeros(nG,nE);
    gamma = zeros(1,nE);
        
    for iE = 1:nE
        for iG = 1:nG
            alpha(iG,iE) = P(idcs.alpha{iG}(iE));
            beta(iG,iE) = P(idcs.beta{iG}(iE));
        end
        gamma(iE) = P(idcs.gamma(iE));
        G(:,iE,:) = circshift(G(:,iE,:),gamma(iE),1);
    end
        
    g_mean = squeeze(nanmean(G,2));
    
    chi2_G = zeros(nG,1);
    for iG = 1:nG
        G_tmp = squeeze(G(:,:,iG));
        alpha_tmp = alpha(iG,:);
        beta_tmp = beta(iG,:);
        chi2_G(iG) = ...
            nanmean(nansum( (G_tmp - (alpha_tmp + beta_tmp.*g_mean(:,iG))).^2 ));
    end
    
    chi2 = mean(chi2_G);
end