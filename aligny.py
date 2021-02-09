import numpy as np
import matplotlib.pyplot as plt
#import scipy
from scipy.optimize import minimize
#from scipy.optimize import NonlinearConstraint
from anchormean0to1 import anchormean0to1

def aligny(G_0):
    
    print("Beginning Y alignment")
    
    N = G_0.shape[1]
    alpha_0 = np.zeros([1,N])
    beta_0 = np.ones([1,N])
    G_0 = (G_0 - alpha_0)/beta_0
    G_mean_0 = np.nanmean(G_0,1)
    

    p_0 = np.squeeze(np.array([np.transpose(alpha_0),np.transpose(beta_0)]))
    #p_0 = ([[alpha_0],[beta_0]])
    
    # alpha lower bound is -1, upper bound is 1
    # beta lower bound is 0,   upper bound is inf

    lb_alpha = np.ones([N,1])*-1
    lb_beta = np.zeros([N,1])
    ub_alpha = np.ones([N,1])
    ub_beta = np.ones([N,1])*np.inf
    
    lb = np.vstack([lb_alpha, lb_beta])
    ub = np.vstack([ub_alpha, ub_beta])
    
    bnds = np.hstack([lb,ub])
    
    ### Working to here -MJT
    fun = lambda params: chisq(params, G_0, G_mean_0, N)
    solution = minimize(fun, p_0, method='SLSQP', bounds=bnds, constraints=(), 
                        tol=1e-6, callback=None, options=None)
    #                    options={'gtol': 1e-6, 'disp': True})
    
    #solution = minimize(chisq(p_0,G_0,G_mean_0), p_0, method='SLSQP', bounds=None, constraints=(), tol=None, callback=None, options=None)
    
    # -> line 618 _minimize.py -> line 399 slsqp.py -> line 327 optimize.py -> line 26 aligny.py -> line 52 aligny.py
    # *** IndexError: too many indices for array
    
    #f = chisq(p_0, G_0, G_mean_0)
    #[p,chi2,exitFlag,~] = fmincon(f,p_0,[],[],[],[],lb,ub,[],options);
    
    #scipy.optimize.NonLinearConstraint(f,lb,ub)


    alpha = solution.x[0:N]
    beta = solution.x[N:2*N]

    g = (G_0 - alpha)/beta
    
    return g



#Cost Function

def chisq(params,G,g_mean,N):
    alpha = params[0:N]
    beta = params[N:2*N]
    chi2 = np.nanmean(np.nansum ((G - (alpha + beta))**2))
    
    return chi2