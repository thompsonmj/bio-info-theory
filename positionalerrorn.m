function sigma_x = positionalerrorn(Y_in)

% Based on Eq. 11 in Dubuis et al. 2013.
% ref. https://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm

[nX,nE,nG] = size(Y_in);

Y_mean = nanmean(Y_in,2);

%% Calculate the local slope of the mean of each profile
dgdx = zeros(nX-1,nG);
for iG = 1:nG
    dgdx(:,iG) = diff(Y_mean(:,:,iG));
end

%% Calculate the nG x nG covariance matrix at each position
C = zeros(nG,nG,nX);
for iX = 100:nX
    for iG = 1:nG
        for jG = 1:nG
            if iG == jG
                C(iG,jG,iX) = nanvar(Y_in(iX,:,iG),0,'all');
            else
                C(iG,jG,iX) = mean(Y_in(iX,:,iG) .* Y_in(iX,:,jG)) - ...
                    mean(Y_in(iX,:,iG)) * mean(Y_in(iX,:,jG));
            end
            % Alternatively
            % ref. https://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
            C(iG,jG,iX) = ...
                (1/(nE-1))*sum((Y_in(iX,:,iG)-mean(Y_in(iX,:,iG))).*(Y_in(iX,:,jG)-mean(Y_in(iX,:,jG))));
        end
    end
end

%% Calculate the positional error
C_inv = zeros(nG,nG,nX);
var_inv = zeros(nX,1);
sigma_x = zeros(nX,1);
for iX = 1:nX-1
    C_inv(:,:,iX) = inv(C(:,:,iX));
    for iG = 1:nG
        for jG = 1:nG
            var_inv(iX) = var_inv(iX) + ...
                dgdx(iX,iG) * C_inv(iG,jG,iX) * dgdx(iX,jG);
        end
    end
    sigma_x(iX) = sqrt(1/var_inv(iX))/nX;
end

end