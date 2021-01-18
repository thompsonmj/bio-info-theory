function P_map_alpha = decodingmapn(Y, n_bins_x, n_bins_y)

[n_points, n_replicates, n_targets] = size(Y);
y_min = min(Y,[],'all');
y_max = max(Y,[],'all');
Y_mean = mean(Y,2);

Z = 1; %n_points;

%% Px(x)
px_bin = (1/n_bins_x);
px = repmat(px_bin, n_bins_x, 1);

x = 1:n_points;
X = repmat(x', 1, n_replicates*n_targets);
x_column = reshape(X,numel(X),1);

Y_column = reshape(Y,numel(Y),1);

%% Calculate joint distribution P(g,x)
[pyx_joint, pyx_joint_edges] = ...
    histcountsn([x_column, Y_column],[n_bins_x, n_bins_y]);
pyx_joint = pyx_joint + eps;

bin_width_x = mean(diff(pyx_joint_edges{1}));
bin_width_y = mean(diff(pyx_joint_edges{2}));

%% dyP(y|x) = 1
for i_x = 1:n_bins_x
    pyx_joint(i_x,:) = ...
        pyx_joint(i_x,:)/sum(pyx_joint(i_x,:));
end
%% Calculate P(y)
% Sum at each level across all x.
py = sum(pyx_joint,1);
% ?dyP(y) = 1
py = py/sum(py);

if n_targets == 1
    %% Single Gene Case
    disp('Mapping single gene') 
    P_map_alpha = zeros(n_bins_x, n_bins_x, n_replicates);
    Y_std = std(Y,0,2); 
        % Arg. 2: weight. See "Bessel's correction"
        %   0 -- Normalize by N-1, N is no. observations.
        %   1 -- Normalize by N
        %   Vector made up of nonnegative scalar weights corresponding to
        %   the dimension of Y along which the standard deviation is
        %   calculated
    for i_r = 1:n_replicates
        for i_x_imp = 1:n_bins_x % Implied position (imp)
            
            y_mean_imp = Y_mean(i_x_imp);
            y_std_imp = Y_std(i_x_imp);
                        
            for i_x_act = 1:n_bins_x % Actual position (act)
                
                y_alpha_act = Y(i_x_act,i_r);
                
%                 idx_y = max([round( (y_alpha_act - y_min)/bin_width_y),1]);
%                 Z = py(idx_y);
                
                chi2_imp = (y_alpha_act - y_mean_imp)^2 / (y_std_imp^2);
                
                p_y_given_imp = (2 * pi * y_std_imp^2)^(-1/2)* ...
                    exp(-chi2_imp/2);
                
                p_imp_given_y_alpha = 1/Z * p_y_given_imp * px(i_x_imp); 
                
                p_imp_given_act_alpha = p_imp_given_y_alpha;
                
                P_map_alpha(i_x_imp,i_x_act,i_r) = p_imp_given_act_alpha;
                
            end
            P_temp = P_map_alpha(i_x_imp,:,i_r);
            P_temp = P_temp/sum(P_temp)/(bin_width_x*bin_width_y);

            P_map_alpha(i_x_imp,:,i_r) = P_temp;
        end
    end
    
elseif n_targets > 1
    %% Multiple Gene Case
    disp(['Mapping ',num2str(n_targets),' genes'])   
    
    Y_mean = squeeze(mean(Y,2)); % n_points x n_targets
    
    % Estimate the nTargets x nTargets covariance matrix at each position.
    P_map_alpha = zeros(n_bins_x, n_bins_x, n_replicates);
    
    C = zeros(n_bins_x, n_targets, n_targets, n_replicates);
    
    for i_x = 1:n_bins_x
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
    
    % n_bins_x X n_targets X n_targets
    C_mean = squeeze(nanmean(C,4));
    C_inv = zeros(n_bins_x, n_targets, n_targets);
    C_det = zeros(n_bins_x,1);
    for i_x = 1:n_bins_x
        C_inv(i_x,:,:) = inv(squeeze(C_mean(i_x,:,:)));
        C_det(i_x) = det(squeeze(C_mean(i_x,:,:)));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_map_alpha = zeros(n_bins_x, n_bins_x, n_replicates);
    Y_std = std(Y,0,2); 
        % Arg. 2: weight. See "Bessel's correction"
        %   0 -- Normalize by N-1, N is no. observations.
        %   1 -- Normalize by N
        %   Vector made up of nonnegative scalar weights corresponding to
        %   the dimension of Y along which the standard deviation is
        %   calculated
    for i_r = 1:n_replicates
        for i_x_imp = 1:n_bins_x % Implied position (imp)
            
            y_mean_imp = Y_mean(i_x_imp);
            y_std_imp = Y_std(i_x_imp);
                        
            for i_x_act = 1:n_bins_x % Actual position (act)
                
                chi2_imp = zeros(n_targets, n_targets);
                
               
                % ?^2 for multiple genes
                for i_t = 1:n_targets
                    y_i_alpha_act = Y(i_x_act,i_r,i_t);
                    y_i_mean_imp = Y_mean(i_x_imp,i_t);
                    for j_t = 1:n_targets
                        y_j_alpha_act = Y(i_x_act,i_r,j_t);
                        y_j_mean_imp = Y_mean(i_x_imp,j_t);
                        chi2_imp(i_t,j_t) = ...
                            (y_i_alpha_act - y_i_mean_imp)* ...
                            C_inv(i_x_imp,i_t,j_t)* ...
                            (y_j_alpha_act - y_j_mean_imp);
                    end
                end
                chi2_imp = sum(chi2_imp,'all');
                
                p_y_given_imp = ((2*pi)^n_targets * C_det(i_x_imp))^(-1/2)* ...
                    exp(-chi2_imp/2);
                
                p_imp_given_y_alpha = 1/Z * p_y_given_imp * px(i_x_imp);
                
                p_imp_given_act_alpha = p_imp_given_y_alpha;
                
                P_map_alpha(i_x_imp,i_x_act,i_r) = p_imp_given_act_alpha;
                
            end
            P_temp = P_map_alpha(i_x_imp,:,i_r);
            P_temp = P_temp/sum(P_temp)/(bin_width_x*bin_width_y);

            P_map_alpha(i_x_imp,:,i_r) = P_temp;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end