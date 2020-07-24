function [n,edges,binIdcs] = histcountsn(x,nBins,varargin)
%HISTCOUNTSN N-variate histogram bin counts.
%   [N,EDGES] = HISTCOUNTSN(X) partitions the values
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

%% Note
% Minimally working version for use with midirectestimate.m and
% midirectestimate2.m.

%%
nBinsIN = nBins;

opts = parseinput(varargin);
% Filter NaN values from input.
x(any(isnan(x),2),:) = [];
% nnn = numel(x(:,1));

[~,nDims] = size(x);
subs = [];
binIdcs = cell(nDims,1);
edges = cell(nDims,1);

if numel(nBinsIN) == 1
    % One bin size for all dimensions.
    nBins = repmat(nBins,nDims,1);
end

sz = zeros(1,nDims);

for iD = 1:nDims
    
    if numel(nBinsIN) == 0
        % Allow histcounts to autobin.
        [~,edges{iD},binIdcs{iD}] = histcounts(x(:,iD));
        nBins = [nBins, numel(edges{iD}) - 1];
    elseif numel(nBinsIN) == 1
        % One bin size for all dimensions.
        [~,edges{iD},binIdcs{iD}] = histcounts(x(:,iD),nBins(iD));
    elseif numel(nBinsIN) > 1
        % Read bin size for each dimension.
       [~,edges{iD},binIdcs{iD}] = histcounts(x(:,iD),nBins(iD));
    end

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

%% Normalization options
switch opts.Normalization
    
    case 'countdensity'
        edgeL = zeros(nDims,1);
        for iD = 1:nDims
            edgeL(iD) = mean(double(diff(edges{iD})));
        end
        binVolumeN = prod(edgeL);
        n = n / binVolumeN;
        
    case 'cumcount'
        for iD = 1:nDims
            n = cumsum(n,iD);
        end
        
    case 'probability'
        n = n/numel(subs(:,1));
        
    case 'pdf'
        edgeL = zeros(nDims,1);
        for iD = 1:nDims
            edgeL(iD) = mean(double(diff(edges{iD})));
        end
        binVolumeN = prod(edgeL);
        n = n/numel(subs(:,1)) / binVolumeN;
        
    case 'cdf'
        n = n/numel(subs(:,1));
        for iD = 1:nDims
            n = cumsum( n, iD );
        end
end

end

%% LOCAL FUNCTIONS
%%% TAKEN VERBATIM FROM HISTCOUNTS2
function opts = parseinput(input) % Input is varargin (inputs 2+)

% opts = struct('NumBins',[],'BinEdges',{},'BinLimits',{},'BinWidth', ...
%     'Normalization','count','BinMethod','auto');

opts = struct('NumBins',[],'XBinEdges',[],'YBinEdges',[],'XBinLimits',[],...
    'YBinLimits',[],'BinWidth',[],'Normalization','count','BinMethod','auto');
funcname = mfilename;

% Parse third and fourth input in the function call
inputlen = length(input);
if inputlen > 0
    in = input{1};
    inputoffset = 0;
    if isnumeric(in) || islogical(in)
        if inputlen == 1 || ~(isnumeric(input{2}) || islogical(input{2}))
            % Numbins
            if isscalar(in)
                in = [in in];
            end
            validateattributes(in,{'numeric','logical'},{'integer', 'positive', ...
                'numel', 2, 'vector'}, funcname, 'm', inputoffset+3)
            opts.NumBins = in;
            input(1) = [];
            inputoffset = inputoffset + 1;
        else
            % XBinEdges and YBinEdges
            in2 = input{2};
            validateattributes(in,{'numeric','logical'},{'vector', ...
                'real', 'nondecreasing'}, funcname, 'xedges', inputoffset+3)
            if length(in) < 2
                error(message('MATLAB:histcounts2:EmptyOrScalarXBinEdges'));
            end
            validateattributes(in2,{'numeric','logical'},{'vector', ...
                'real', 'nondecreasing'}, funcname, 'yedges', inputoffset+4)
            if length(in2) < 2
                error(message('MATLAB:histcounts2:EmptyOrScalarYBinEdges'));
            end
            opts.XBinEdges = in;
            opts.YBinEdges = in2;
            input(1:2) = [];
            inputoffset = inputoffset + 2;
        end
        opts.BinMethod = [];
    end
    
    % All the rest are name-value pairs
    inputlen = length(input);
    if rem(inputlen,2) ~= 0
        error(message('MATLAB:histcounts2:ArgNameValueMismatch'))
    end
    
    for i = 1:2:inputlen
        name = validatestring(input{i}, {'NumBins', 'XBinEdges', ...
            'YBinEdges','BinWidth', 'BinMethod', 'XBinLimits', ...
            'YBinLimits','Normalization'}, i+2+inputoffset);
        
        value = input{i+1};
        switch name
            case 'NumBins'
                if isscalar(value)
                    value = [value value]; %#ok
                end
                validateattributes(value,{'numeric','logical'},{'integer', ...
                    'positive', 'numel', 2, 'vector'}, funcname, 'NumBins', i+3+inputoffset)
                opts.NumBins = value;
                if ~isempty(opts.XBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedXBinInputs'))
                elseif ~isempty(opts.YBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedYBinInputs'))
                end
                opts.BinMethod = [];
                opts.BinWidth = [];
            case 'XBinEdges'
                validateattributes(value,{'numeric','logical'},{'vector', ...
                    'real', 'nondecreasing'}, funcname, 'XBinEdges', i+3+inputoffset);
                if length(value) < 2
                    error(message('MATLAB:histcounts2:EmptyOrScalarXBinEdges'));
                end
                opts.XBinEdges = value;
                % Only set NumBins field to empty if both XBinEdges and
                % YBinEdges are set, to enable BinEdges override of one
                % dimension
                if ~isempty(opts.YBinEdges)
                    opts.NumBins = [];
                    opts.BinMethod = [];
                    opts.BinWidth = [];
                end
                opts.XBinLimits = [];
            case 'YBinEdges'
                validateattributes(value,{'numeric','logical'},{'vector', ...
                    'real', 'nondecreasing'}, funcname, 'YBinEdges', i+3+inputoffset);
                if length(value) < 2
                    error(message('MATLAB:histcounts2:EmptyOrScalarYBinEdges'));
                end                
                opts.YBinEdges = value;
                % Only set NumBins field to empty if both XBinEdges and
                % YBinEdges are set, to enable BinEdges override of one
                % dimension
                if ~isempty(opts.XBinEdges)
                    opts.BinMethod = [];
                    opts.NumBins = [];
                    %opts.BinLimits = [];
                    opts.BinWidth = [];
                end
                opts.YBinLimits = [];
            case 'BinWidth'
                if isscalar(value)
                    value = [value value]; %#ok
                end
                validateattributes(value, {'numeric','logical'}, {'real', 'positive',...
                    'finite','numel',2,'vector'}, funcname, ...
                    'BinWidth', i+3+inputoffset);
                opts.BinWidth = value;
                if ~isempty(opts.XBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedXBinInputs'))
                elseif ~isempty(opts.YBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedYBinInputs'))
                end
                opts.BinMethod = [];
                opts.NumBins = [];
            case 'BinMethod'
                opts.BinMethod = validatestring(value, {'auto','scott',...
                    'fd','integers'}, funcname, 'BinMethod', i+3+inputoffset);
                if ~isempty(opts.XBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedXBinInputs'))
                elseif ~isempty(opts.YBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedYBinInputs'))
                end
                opts.BinWidth = [];
                opts.NumBins = [];
            case 'XBinLimits'
                validateattributes(value, {'numeric','logical'}, {'numel', 2, ...
                    'vector', 'real', 'finite','nondecreasing'}, funcname, ...
                    'XBinLimits', i+3+inputoffset)
                opts.XBinLimits = value;
                if ~isempty(opts.XBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedXBinInputs'))
                end
            case 'YBinLimits'
                validateattributes(value, {'numeric','logical'}, {'numel', 2, ...
                    'vector', 'real', 'finite','nondecreasing'}, funcname, ...
                    'YBinLimits', i+3+inputoffset)
                opts.YBinLimits = value;
                if ~isempty(opts.YBinEdges)
                    error(message('MATLAB:histcounts2:InvalidMixedYBinInputs'))
                end
            otherwise % 'Normalization'
                opts.Normalization = validatestring(value, {'count', 'countdensity', 'cumcount',...
                    'probability', 'pdf', 'cdf'}, funcname, 'Normalization', i+3+inputoffset);
        end
    end
end
end