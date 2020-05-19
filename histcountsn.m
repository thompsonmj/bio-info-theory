function [p,edges] = histcountsn(x,nBins,varargin)
%HISTCOUNTSN N-variate histogram bin counts.
%   [N,X1EDGES,X2EDGES,...,XNEDGES] = HISTCOUNTSN(X) partitions the values
%   in columns of X into bins, and returns the count in each bin, as well
%   as the bin edges. HISTCOUNTSN determines the bin edges using an
%   automatic binning algorithm that returns uniform bins chosen to cover
%   the range of values in each column of X and reveal the shape of the
%   underlying distribution.
%%
%   N an I1-by-I2-by-...-by-IN matrix where I1 through IN are the number of
%   bins along the X1 through XN dimensions respectively. ...

%
%   N is an I-by-J matrix where I and J are the number of bins along the
%   X and Y dimensions respectively. N(i,j) will count the value [X(k),Y(k)]
%   if XEDGES(i) <= X(k) < XEDGES(i+1) and YEDGES(j) <= Y(k) < YEDGES(j+1).
%   The last bins in the X and Y dimensions will also include the upper
%   edge. For example, [X(k),Y(k)] will fall into the i-th bin in the last
%   row if XEDGES(end-1) <= X(k) <= XEDGES(end) &&
%   YEDGES(i) <= Y(k) < YEDGES(i+1).
%%
%   [N,X1EDGES,X2EDGES,...,XNEDGES] = HISTCOUNTSN(X,NBINS) where NBINS is a
%   scalar or N-element vector, specifies the number of bins to use. A
%   scalar specifies the same number of bins in each dimension, whereas the
%   N-element vector [nbinsx1 nbinsx2 ... nbinsxn] specifies a different
%   number of bins for the X1 through XN dimensions.

%%
% Filter NaN values from input.
% pyxJ2 = histcounts2(x(:,1),x(:,2),nBins);
% pyxJ2_norm = histcounts2(x(:,1),x(:,2),nBins,'Normalization','probability');


x(any(isnan(x),2),:) = [];
% nnn = numel(x(:,1));

[~,nDims] = size(x);
subs = [];
binIdcs = cell(nDims,1);
edges = cell(nDims,1);

if numel(nBins) == 0
    % Allow histcounts to autobin.
    nBins = [];
elseif nBins == 1
    % One bin size for all dimensions.
    nBins = repmat(nBins,nDims,1);
elseif nBins > 1
    % Read bin size for each dimension.
end

sz = zeros(1,nDims);

for iD = 1:nDims
    [~,edges{iD},binIdcs{iD}] = histcounts(x(:,iD),nBins(iD));
    subs_tmp = binIdcs{iD};
    % Filter out-of-range data (bin index = 0).
    subs(any(subs_tmp==0,2),:) = [];
    subs = [subs, subs_tmp];
    sz(iD) = repmat(nBins(iD),1,1);
end

if nDims == 1
    sz = [sz 1];
end

n = accumarray(subs,ones(size(subs,1),1),sz);
p = n/numel(subs(:,1));

end