function [y, alpha, beta, chi2] = aligny(Y_0)
%ALIGNY aligns data along the Y-axis of data in a 2D matrix.
% 
% Input:
%   > Y_0: initial profile data. 2D (no. data points) x (no. replicates).
%
% Output:
%   > y: aligned output data. 2D. Same size as input.
%   > alpha: additive parameter (1 per profile)
%   > beta: scaling/multiplicative parameter (1 per profile)
%   > chi2: final squared error

disp('Beginning Y alignment')

[nPoints, nReplicates] = size(Y_0);

%% Initialize alignment parameters.
alpha_0 = zeros(1,nReplicates);
beta_0 = ones(1,nReplicates);

Y_0 = (Y_0 - alpha_0)./beta_0;
Y_mean_0 = nanmean(Y_0,2);

p_0 = [alpha_0;beta_0];

%% Set optimization options.
options = optimset( ...
    'MaxIter', 1e6, ...
    'MaxFunEvals', 1e6, ...
    'Display','notify-detailed');

%% Set optimization parameter constraints.
lb = [-ones(1,nReplicates);zeros(1,nReplicates)];
ub = [ones(1,nReplicates);inf(1,nReplicates)];

%% Optimization.
f = @(p_0)chisq(p_0, Y_0, Y_mean_0);
[p,chi2,exitFlag,~] = fmincon(f,p_0,[],[],[],[],lb,ub,[],options);

exitFlag

%% Apply final optimization parameters to profile data.
alpha = p(1,:);
beta = p(2,:);

y = (Y_0 - alpha)./beta;

end

%% Cost function
function chi2 = chisq(p,Y,y_mean)
   alpha = p(1,:);
   beta = p(2,:);
   chi2 = nanmean(nansum( (Y - (alpha + beta.*y_mean)).^2 ));
end