function sigma_x = positionalerror(profileSet)

% sigX/L = sigG*abs(dg/dx)^-1
% profileSet: nX x nSamples (vertical)

g = mean(profileSet,2);

dgdx = diff(g);

sigG = std(profileSet,1,2);
sigG = sigG(1:end-1);

sigma_x = sigG.*abs(dgdx).^-1/numel(dgdx);

end