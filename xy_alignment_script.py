# Generated with SMOP  0.41
from smop.libsmop import *

# xy_alignment_script.m

#clear

    
nEmbryosUsed=102
# xy_alignment_script.m:4

inputData='gap_data_raw_dorsal_wt'
# xy_alignment_script.m:5

inputFile=concat([inputData,'.mat'])

# xy_alignment_script.m:6

#load(inputFile)

nEmbryos=numel(data) 


# xy_alignment_script.m:10

assert_(nEmbryos >= nEmbryosUsed)

    # Authors' description of the data: "Each profile vector is of length 1000
# pixels, where 0 corresponds to the anterior (A) and 1000 to the posterior
# (P) of the embryo."
# To mitigate 'edge effects,' we extract only the middle 80# of the data to 
# use throughout the analysis.
    
nSamplePts=numel(data(1).Hb)
# xy_alignment_script.m:19

x=(arange(1 / nSamplePts,1,1 / nSamplePts))
# xy_alignment_script.m:20

idx_low=dot(round(nSamplePts),0.1)
# xy_alignment_script.m:21

idx_high=dot(round(nSamplePts),0.9)
# xy_alignment_script.m:22

x=x(arange(idx_low,idx_high))
# xy_alignment_script.m:23

genes=cellarray(['Hb','Kr','Gt','Kni'])
# xy_alignment_script.m:25

nGenes=numel(genes)
# xy_alignment_script.m:26
    ## Set up raw data
    
Y_raw=zeros(nSamplePts,nEmbryosUsed,nGenes)
# xy_alignment_script.m:29

idcs=sort(randperm(nEmbryos,nEmbryosUsed))
# xy_alignment_script.m:30

for iG in arange(1,nGenes).reshape(-1):
    for iE in arange(1,nEmbryos).reshape(-1):
        Y_raw[arange(),iE,iG]=getattr(data(iE),(genes[iG]))
# xy_alignment_script.m:34
    
    Y_raw=Y_raw(arange(idx_low,idx_high),idcs,arange())
# xy_alignment_script.m:37
    
    nPts,__,__=size(Y_raw,nargout=3)
# xy_alignment_script.m:38
    
    padSize=round(dot(nSamplePts,0.1))
# xy_alignment_script.m:39
    
    Y_raw=padarray(Y_raw,padSize,NaN,'both')
# xy_alignment_script.m:40
    
    ## Comments
# Up to this point, we have simply deleted data outside of the middle 80#
# of the domain (replacing it with NaNs to act as a buffer to shift the
# profiles to the left/right during X-alignment).
## 1 Align Y alone (optimize alpha and beta).
    
    Y_raw,__=anchormean0to1(Y_raw,nargout=2)
# xy_alignment_script.m:47
    
    Y_yAlign1=zeros(size(Y_raw))
# xy_alignment_script.m:48
    
    for iG in arange(1,nGenes).reshape(-1):
        y=squeeze(Y_raw(arange(),arange(),iG))
# xy_alignment_script.m:50
        
        
       # Y_yAlign1(arange(),arange(),iG),__ = aligny(y,nargout=2)
       
        
# xy_alignment_script.m:51
    
    ## 2 Numerically estimate optimal joint XY alignment (alpha, beta, gamma).
    UseParallel_TF=copy(true)
    
# xy_alignment_script.m:55
    pool=copy(gcp)
    
# xy_alignment_script.m:56
    if UseParallel_TF:
        if isempty(gcp('nocreate')):
            parpool(pool.NumWorkers)
    
    Y_xyAlign,alpha,beta,gamma=alignxy(Y_yAlign1,UseParallel_TF,nargout=4)
    
# xy_alignment_script.m:62
    # Should explore a 'better' solution than genetic optimization. This was
# chosen because I shift profiles in x using their integer indices, and
# MATLAB's genetic optimizer allowed for easily implemented integer
# constraints on one parameter (gamma) and continuous variation on others
# (alpha and beta). This approach does not yield the same result twice, but
# reliably yields a better x-alignment than the data came with originally.
# However, it needs tweaking for the data at hand (e.g. the NaN buffer size
# on either side of each profile should be large enough, but should not be
# too large). Additionally, y-alignment is poor when done jointly using
# this method, so an additional chi^2 minimization must be run afterward.
    
    ## 3 Align Y alone once more for closed-form Y alignment.
    
    Y_align=zeros(size(Y_yAlign1))
# xy_alignment_script.m:75
    
    for iG in arange(1,nGenes).reshape(-1):
        y=squeeze(Y_xyAlign(arange(),arange(),iG))
# xy_alignment_script.m:77
        
        
       # Y_align(arange(),arange(),iG),__=aligny(y,nargout=2)
        
        
# xy_alignment_script.m:78
    
    save(concat([inputData,'_XY-Aligned.mat']))