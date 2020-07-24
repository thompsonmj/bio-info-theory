function sigma_x = positionalerror(Y_in)
%POSITIONALERROR calculates the local positional error for a single
%profile.

% Input:
%   > Y_in: size: (number of x-points) x (number of replicates)
% Output:
%   > sigma_x: 1-D array of the positional error (not normalized)

% Based on Eq. 11 in Dubuis et al. 2013.
% ref. https://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm

g = mean(Y_in,2);

dgdx = diff(g);

sigG = std(Y_in,1,2);
sigG = sigG(1:end-1);

sigma_x = sigG.*abs(dgdx).^-1/numel(dgdx);

end