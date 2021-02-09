function MI_est = ...
    midirectestimaten(Y, n_trials, n_boots, bin_counts, subsamps)
%MIDIRECTESTIMATEN computes a direct estimate of mutual information through
%iterative naive estimates.
% 
% Inputs
%   > Y: 3D matrix containing profile data. 
%       Dims: (no. sampled points) x (no. replicates) x (no. targets)
%   > n_trials: No. trials for bootstrapping
%   > n_boots: number of 
%   > binCountArray: 
%   > subsampleArray:
% Outputs
%   > MI_est: nTrials x nBinSizes x nSubsamplesets x nBoots array of naive
%   estimates. To be extrapolated to infinitte data and 0 bin size.
%
% Based on Eq. 7 and Eq. 24 in Tkacik et al, 2015.

[n_points, n_replicates, n_targets] = size(Y);
x = 1/n_points:1/n_points:1;
% nX = 100; % For pSmad data with such small sample size.

m = round(subsamps*n_replicates);
n_subsample_sets = numel(m);

replacement = false;

replicate_idcs = [1:1:n_replicates];

n_bin_sizes = numel(bin_counts);

perm_order = [linspace(2,n_targets+1,n_targets),1];

for i_trial = 1:n_trials
    % Initialize struct array entry for this trial
    MI_est(i_trial).binSizes = bin_counts;
    MI_est(i_trial).subSamps = subsamps;
    MI_est(i_trial).nBoots = n_boots;
    
    MI_est(i_trial).naiveEst_means = zeros(n_bin_sizes,n_subsample_sets);
    MI_est(i_trial).naiveEst_stds = zeros(n_bin_sizes,n_subsample_sets);
    
    % Begin bootstrapping across bin sizes
    for i_bin_idx = 1:n_bin_sizes
        n_bins_y = bin_counts(i_bin_idx);
        n_bins_x = n_bins_y;
        
        py_input = repmat(n_bins_y, 1, n_targets);
        py = zeros(py_input);
        px = 1/n_bins_x; % Assume uniform distribution.
        px = repmat(px,[n_bins_x,py_input]);

        for i_subsamp_idx = 1:n_subsample_sets
            m_subsamps = m(i_subsamp_idx);
            pyx_joint = zeros([n_bins_x,py_input]);
            MIEq24 = zeros(n_boots,1);
            for i_b = 1:n_boots
                k = m(i_subsamp_idx);
                %%%Need to generalize sub to match input dims
                subsample_idcs = randsample(replicate_idcs, k, replacement);
                subsampled_replicates_temp = Y(:,subsample_idcs,:);
                % ???
                subsampled_replicates = ...
                    reshape(subsampled_replicates_temp, ...
                    numel(subsampled_replicates_temp(:,:,1)), ...
                    n_targets);    
                
                x_mat = repmat(x',1,m_subsamps);
                x_column = reshape(x_mat,numel(x_mat),1);
                % Distributions
                pyx_joint = histcountsn([x_column, subsampled_replicates], ...
                    [n_bins_x, py_input], ...
                    'Normalization','probability');
                pyx_joint = pyx_joint + eps;

                pyx_joint_perm = permute(pyx_joint,perm_order);


                py = squeeze(sum(pyx_joint,1));
                py_given_x = pyx_joint ./ px;
                py_given_x_perm = permute(py_given_x,perm_order);

                log_value = py_given_x_perm .* log2(py_given_x_perm ./ py);
                MIEq7_temp = sum(px .* log_value, 'all');
                MIEq7(i_b) = MIEq7_temp;
                
                MIEq24_temp = ...
                    sum( pyx_joint_perm .* log2(pyx_joint_perm ./ (px .* py)),'all' );
                MIEq24(i_b) = MIEq24_temp;
                
                assert(abs(MIEq7_temp - MIEq24_temp) < 1e-10)
            end
            MI_est(i_trial).naiveEst_means(i_bin_idx, i_subsamp_idx) = ...
                mean(MIEq24);
            MI_est(i_trial).naiveEst_stds(i_bin_idx, i_subsamp_idx) = ...
                std(MIEq24,0);
        end
    end
end

%% More than 2 genes
% To extend this method tractably to more than two genes, one needs to
% resort to approximations for P({gi}|x), the simplest one of which is the
% so-called Gaussian approximation, shown in Equation 4.
  
% Estimate the nTargets x nTargets covariance matrix at each position.


% C = zeros(n_bins_x, n_targets, n_targets, n_replicates);
% 
% for i_x = 1:n_bins_x
%     for i_t = 1:n_targets
%         Y_i_mean = mean(Y(i_x,:,i_t),2);
%         for j_t = 1:n_targets
%             Y_j_mean = mean(Y(i_x,:,j_t),2);
%             for i_r = 1:n_replicates
%                 y_i_alpha = Y(i_x,i_r,i_t);
%                 y_j_alpha = Y(i_x,i_r,j_t);
% 
%                 C(i_x, i_t, j_t, i_r) = ...
%                     (y_i_alpha - Y_i_mean)* ...
%                     (y_j_alpha - Y_j_mean);
%             end
%         end
%     end
% end
%     % C_mean: n_bins_x X n_targets X n_targets
% C_mean = squeeze(nanmean(C,4));
% C_inv = zeros(n_bins_x, n_targets, n_targets);
% C_det = zeros(n_bins_x,1);
% for i_x = 1:n_bins_x
%     C_inv(i_x,:,:) = inv(squeeze(C_mean(i_x,:,:)));
%     C_det(i_x) = det(squeeze(C_mean(i_x,:,:)));
% end
% 
% for i_x = 1:n_bins_x
%     py_given_x = (2*pi)^(n_targets/2)*(C_mean(i_x)).^(-1/2)* ...
%         exp(-1/2*());
% end