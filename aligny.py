# Generated with SMOP  0.41
from smop.libsmop import *
# aligny.m

    
@function
def aligny(G_0=None,*args,**kwargs):
    varargin = aligny.varargin
    nargin = aligny.nargin

    tic
    disp('Beginning Y alignment')
    __,N=size(G_0,nargout=2)
# aligny.m:6
    alpha_0=zeros(1,N)
# aligny.m:8
    beta_0=ones(1,N)
# aligny.m:9
    G_0=(G_0 - alpha_0) / beta_0
# aligny.m:10
    G_mean_0=nanmean(G_0,2)
# aligny.m:11
    p_0=concat([[alpha_0],[beta_0]])
# aligny.m:13
    options=optimset('MaxIter',1000000.0,'MaxFunEvals',1000000.0,'Display','notify-detailed')
# aligny.m:15
    lb=concat([[- ones(1,N)],[zeros(1,N)]])
# aligny.m:20
    ub=concat([[ones(1,N)],[inf(1,N)]])
# aligny.m:21
    f=lambda p_0=None: chisq(p_0,G_0,G_mean_0)
# aligny.m:23
    p,chi2,exitFlag,__=fmincon(f,p_0,[],[],[],[],lb,ub,[],options,nargout=4)
# aligny.m:24
    exitFlag
    alpha=p(1,arange())
# aligny.m:28
    beta=p(2,arange())
# aligny.m:29
    g=(G_0 - alpha) / beta
# aligny.m:31
    g,__=anchormean0to1(g,nargout=2)
# aligny.m:33
    toc
    return g,chi2
    
if __name__ == '__main__':
    pass
    
    ## Cost function
    
@function
def chisq(p=None,G=None,g_mean=None,*args,**kwargs):
    varargin = chisq.varargin
    nargin = chisq.nargin

    alpha=p(1,arange())
# aligny.m:40
    beta=p(2,arange())
# aligny.m:41
    chi2=nanmean(nansum((G - (alpha + multiply(beta,g_mean))) ** 2))
# aligny.m:42
    return chi2
    
if __name__ == '__main__':
    pass
    