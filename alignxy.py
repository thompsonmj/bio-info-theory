import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, NonlinearConstraint
from anchormean0to1 import anchormean0to1

    
def alignxy(G_0):

    print("Beginning joint XY alignment")

    padSize=sum(np.isnan(G_0[:,0])) / 2
    
    
    nE = G_0.shape[1]
    nG = G_0.shape[2]
    
    alpha_0 = np.zeros([nG*nE,1]) # One alpha for each profile
    beta_0 = np.ones([nG*nE,1]) # One beta for each profile
    gamma_0 = np.zeros([nE,1]) # One gamma for each embryo
    
    #G_mean_0 = np.nanmean(G_0,1)
    
    nAlpha = np.size(alpha_0)
    nBeta = np.size(beta_0)
    nGamma = np.size(gamma_0)
    
    # Combine parameters into a column array
    p_0 = np.vstack([alpha_0,beta_0,gamma_0])
    
    fun = lambda params: chisq(params, G_0)
    solution = minimize(fun, p_0, method='Nelder-Mead', tol=1e-6)
    
    
def chisq(params, G):
    nP = params.shape[0]
    
    
    # G: 1000 x 102 x 4
    nE = G.shape[1]
    nG = G.shape[2]
    
    alpha = params[0 : nG*nE]
    beta = params[nG*nE : 2*nG*nE]
    gamma = params[2*nG*nE : nP]
    
    iP = 0
    for iE in range(0, nE):
        G[:,iE,:] = np.roll(G[:,iE,:], gamma[iE])
        
    g_mean = np.squeeze(np.nanmean(G,axis=1))
    
    chi2_G = np.zeros([nG,1])
    for iE in range(0, nE):
        for iG in range(0, nG):
            G_tmp = np.squeeze(G[:,:,iG])
            
            # Working on getting param indexing right
            alpha_tmp = alpha[iE*nG - nG : iE*nG]
            beta_tmp = beta[iE*nG - nG : iE*nG]
            
            chi2_G[iG] = np.nanmean(np.nansum( G_tmp - (alpha_tmp + beta_tmp*g_mean[:,iG])**2 ))
            

            
            
    
    #np.roll == circshift
    
    
    
    
    

# =============================================================================
# # alignxy.m:13
#     nAlpha=numel(alpha_0)
# # alignxy.m:15
#     nBeta=numel(beta_0)
# # alignxy.m:16
#     nGamma=numel(gamma_0)
# # alignxy.m:17
#     options=optimoptions('ga','UseParallel',UseParallel_TF,'UseVectorized',false)
# # alignxy.m:19
#     nP=nAlpha + nBeta + nGamma
# # alignxy.m:23
#     ## Specify indices to rereference parameters in cost function.
#     for iG in arange(1,nG).reshape(-1):
#         alphaStart=dot(nAlpha,(iG / nG)) - nE + 1
# # alignxy.m:27
#         alphaEnd=alphaStart + nE - 1
# # alignxy.m:28
#         idcs.alpha[iG]=arange(alphaStart,alphaEnd)
# # alignxy.m:29
#         betaStart=nAlpha + dot(nBeta,(iG / nG)) - nE + 1
# # alignxy.m:31
#         betaEnd=betaStart + nE - 1
# # alignxy.m:32
#         idcs.beta[iG]=arange(betaStart,betaEnd)
# # alignxy.m:33
#     
#     idcs.gamma = copy(arange(nP - nE + 1,nP))
# # alignxy.m:35
#     f=lambda P=None: chisq(P,G_0,idcs)
# # alignxy.m:37
#     IntCon=arange((nP - nGamma + 1),nP)
# # alignxy.m:39
#     lb_alpha=- ones(nAlpha,1)
# # alignxy.m:41
#     ub_alpha=ones(nAlpha,1)
# # alignxy.m:42
#     lb_beta=zeros(nBeta,1)
# # alignxy.m:44
#     ub_beta=dot(2,ones(nBeta,1))
# # alignxy.m:45
#     lb_gamma=repmat(- padSize,nE,1)
# # alignxy.m:47
#     ub_gamma=repmat(padSize,nE,1)
# # alignxy.m:48
#     lb=concat([[lb_alpha],[lb_beta],[lb_gamma]])
# # alignxy.m:50
#     ub=concat([[ub_alpha],[ub_beta],[ub_gamma]])
# # alignxy.m:51
#     P_final,__,exitFlag,__,__=ga(f,nP,[],[],[],[],lb,ub,[],IntCon,options,nargout=5)
# # alignxy.m:53
#     exitFlag
#     ## Extract parameters from P_final.
#     alpha=zeros(nG,nE)
# # alignxy.m:59
#     beta=zeros(nG,nE)
# # alignxy.m:60
#     gamma=zeros(1,nE)
# # alignxy.m:61
#     for iE in arange(1,nE).reshape(-1):
#         for iG in arange(1,nG).reshape(-1):
#             alpha[iG,iE]=P_final(idcs.alpha[iG](iE))
# # alignxy.m:65
#             beta[iG,iE]=P_final(idcs.beta[iG](iE))
# # alignxy.m:66
#         gamma[iE]=P_final(idcs.gamma(iE))
# # alignxy.m:68
#     
#     ## Apply parameters alpha, beta, and gamma to G_0.
#     G=zeros(size(G_0))
# # alignxy.m:72
#     g=zeros(size(G))
# # alignxy.m:73
#     for iE in arange(1,nE).reshape(-1):
#         G[arange(),iE,arange()]=circshift(G_0(arange(),iE,arange()),gamma(iE),1)
# # alignxy.m:75
#     
#     for iG in arange(1,nG).reshape(-1):
#         alpha_tmp=alpha(iG,arange())
# # alignxy.m:78
#         beta_tmp=beta(iG,arange())
# # alignxy.m:79
#         g[arange(),arange(),iG]=(G(arange(),arange(),iG) - alpha_tmp) / beta_tmp
# # alignxy.m:80
#     
#     g,__=anchormean0to1(g,nargout=2)
# # alignxy.m:83
#     shiftAmt=- round((mean(gamma)))
# # alignxy.m:85
#     for iE in arange(1,nE).reshape(-1):
#         g[arange(),iE,arange()]=circshift(g(arange(),iE,arange()),shiftAmt,1)
# # alignxy.m:87
#     
#     toc
#     return g,alpha,beta,gamma
#     
# if __name__ == '__main__':
#     pass
#     
#     ## Cost function
#     
# def chisq(P=None,G=None,idcs=None,*args,**kwargs):
#     varargin = chisq.varargin
#     nargin = chisq.nargin
# 
#     __,nE,nG=size(G,nargout=3)
# # alignxy.m:95
#     alpha=zeros(nG,nE)
# # alignxy.m:97
#     beta=zeros(nG,nE)
# # alignxy.m:98
#     gamma=zeros(1,nE)
# # alignxy.m:99
#     for iE in arange(1,nE).reshape(-1):
#         for iG in arange(1,nG).reshape(-1):
#             alpha[iG,iE]=P(idcs.alpha[iG](iE))
# # alignxy.m:103
#             beta[iG,iE]=P(idcs.beta[iG](iE))
# # alignxy.m:104
#         gamma[iE]=P(idcs.gamma(iE))
# # alignxy.m:106
#         G[arange(),iE,arange()]=circshift(G(arange(),iE,arange()),gamma(iE),1)
# # alignxy.m:107
#     
#     
#     g_mean=squeeze(nanmean(G,2))
# # alignxy.m:110
#     chi2_G=zeros(nG,1)
# # alignxy.m:112
#     for iG in arange(1,nG).reshape(-1):
#         G_tmp=squeeze(G(arange(),arange(),iG))
# # alignxy.m:114
#         alpha_tmp=alpha(iG,arange())
# # alignxy.m:115
#         beta_tmp=beta(iG,arange())
# # alignxy.m:116
#         chi2_G[iG]=nanmean(nansum((G_tmp - (alpha_tmp + multiply(beta_tmp,g_mean(arange(),iG)))) ** 2))
# # alignxy.m:117
#     
#     
#     chi2=mean(chi2_G)
# # alignxy.m:121
#     return chi2
#     
# if __name__ == '__main__':
#     pass
# =============================================================================
    