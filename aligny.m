function [g,chi2] = aligny(G_0)

disp('Beginning Y alignment')

[~, N] = size(G_0);
 
alpha_0 = zeros(1,N);
beta_0 = ones(1,N);
G_0 = (G_0 - alpha_0)./beta_0;
G_mean_0 = nanmean(G_0,2);

p_0 = [alpha_0;beta_0];

options = optimset( ...
    'MaxIter', 1e6, ...
    'MaxFunEvals', 1e6, ...
    'Display','notify-detailed');

lb = [-ones(1,N);zeros(1,N)];
ub = [ones(1,N);inf(1,N)];

f = @(p_0)chisq(p_0, G_0, G_mean_0);
[p,chi2,exitFlag,~] = fmincon(f,p_0,[],[],[],[],lb,ub,[],options);

exitFlag

alpha = p(1,:);
beta = p(2,:);

g = (G_0 - alpha)./beta;

end

%% Cost function
function chi2 = chisq(p,G,g_mean)
   alpha = p(1,:);
   beta = p(2,:);
   chi2 = nanmean(nansum( (G - (alpha + beta.*g_mean)).^2 ));
end