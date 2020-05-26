# Generated with SMOP  0.41
from smop.libsmop import *
# anchormean0to1.m

    
@function
def anchormean0to1(Y=None,*args,**kwargs):
    varargin = anchormean0to1.varargin
    nargin = anchormean0to1.nargin

    #ANCHORMEAN0TO1 anchors the mean of a profile set from 0 to 1.
    
    # Input:
#   > Y: profiles to anchor (nPts x nSamples x nTargets).
    
    # Output:
#   > Y_anchor: profiles scaled such that the mean ranges from 0 to 1.
    
    Y_mean=nanmean(Y,2)
# anchormean0to1.m:10
    Y_anchor=(Y - min(Y_mean)) / (max(Y_mean) - min(Y_mean))
# anchormean0to1.m:11
    Y_mean=nanmean(Y_anchor,2)
# anchormean0to1.m:12
    return Y_anchor,Y_mean
    
if __name__ == '__main__':
    pass
    