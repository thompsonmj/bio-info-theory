function sigma_x = positionalerror(Y_in)
%POSITIONALERROR calculates the local positional error for a single
%profile.

% Input:
%   > Y_in: size: (number of x-points) x (number of replicates)
% Output:
%   > sigma_x: 1-D array of the positional error (not normalized)

% Based on Eq. 11 in Dubuis et al. 2013.

[n_points, n_replicates] = size(Y_in);

g = mean(Y_in,2);

dgdx = diff(g)./(1/(n_points));

sigma_g = std(Y_in,0,2);
sigma_g = sigma_g(1:end-1);

sigma_x = sigma_g.*(dgdx.^-1);

end