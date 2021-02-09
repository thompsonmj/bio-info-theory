function sigma_x = positionalerrorn(Y)
%POSITIONALERRORN calculates the local positional error for a set of N
%profiles.

% Input:
%   > Y_in: size: (number of x-points) x (number of replicates) x (number of
%   independent signals)
% Output:
%   > sigma_x: 1-D array of the positional error (not normalized)

% Based on Eq. 11 in Dubuis et al. 2013.
% ref. https://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm

[n_points, n_replicates, n_targets] = size(Y);

Y_mean = nanmean(Y,2);

%% Calculate the local slope of the mean of each profile
dgdx = zeros(n_points-1,n_targets);
for i_t = 1:n_targets
    dgdx(:,i_t) = diff(Y_mean(:,:,i_t))./(1/(n_points));
end

%% Estimate the nG x nG covariance matrix at each position
Y_mean = squeeze(mean(Y,2)); % n_points x n_targets

C = zeros(n_points, n_targets, n_targets, n_replicates);
for i_x = 1:n_points
    for i_t = 1:n_targets
        Y_i_mean = mean(Y(i_x,:,i_t),2);
        for j_t = 1:n_targets
            Y_j_mean = mean(Y(i_x,:,j_t),2);
            for i_r = 1:n_replicates
                y_i_alpha = Y(i_x,i_r,i_t);
                y_j_alpha = Y(i_x,i_r,j_t);
                C(i_x, i_t, j_t, i_r) = ...
                    (y_i_alpha - Y_i_mean)* ...
                    (y_j_alpha - Y_j_mean);   
            end
        end
    end
end
C_mean = squeeze(nanmean(C,4)); % n_bins_x X n_targets X n_targets
C_inv = zeros(n_points, n_targets, n_targets);
for i_x = 1:n_points
    C_inv(i_x,:,:) = inv(squeeze(C_mean(i_x,:,:)));
end

%% Calculate the positional error
var_inv = zeros(n_points,1);
sigma_x = zeros(n_points,1);
for i_x = 1:n_points-1
    for i_t = 1:n_targets
        for j_t = 1:n_targets
            var_inv(i_x) = var_inv(i_x) + ...
                dgdx(i_x,i_t) * C_inv(i_x,i_t,j_t) * dgdx(i_x,j_t);
        end
    end
    sigma_x(i_x) = ((var_inv(i_x)))^(-1/2)/numel(dgdx);
    
end

end