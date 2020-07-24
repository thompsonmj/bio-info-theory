import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy.optimize import minimize


def anchormean0to1(Y):
    
#ANCHORMEAN0TO1 anchors the mean of a profile set from 0 to 1.
#
# Input:
#   > Y: profiles to anchor (nPts x nSamples x nTargets).
# Output:
#   > Y_anchor: profiles scaled such that the mean ranges from 0 to 1.
#   > Y_mean: mean profile, min = 0 and max = 1.
    
    Y_mean = np.nanmean(Y,axis=1)
    Y_anchor = ((Y - np.nanmin(Y_mean))/((np.nanmax(Y_mean)) - (np.nanmin(Y_mean))))    
    Y_mean = np.nanmean(Y_anchor,axis=1)


    return Y_anchor,Y_mean


def aligny(G_0):
    #tic
    
    print("Beginning Y alignment")
    
    G_0 = np.swapaxes(G_0,0,1)
    
    N = G_0.shape[1] #N=embryos
    
    alpha_0 = np.zeros((N), dtype=int)
    beta_0 = np.ones((N), dtype=int)
    G_0 = (G_0 - alpha_0)/beta_0 # G_0 should be 801xN
    G_mean_0 = np.nanmean(G_0,axis=1)
  #  G_mean_0 = G_mean_0[:,None] #G_mean_0 should be 801x1
    
    p_0 = np.array([[alpha_0],[beta_0]])
    p_0 = np.squeeze(p_0,axis=1) #p_0 should be 2x3 array
   
    print(p_0)
  
    
    lb = ([[np.negative(np.ones((1,N)))],[np.zeros((1,N))]]) 
    #ub = ([[np.ones((1,N))],[(np.ones((1,N)))*np.inf]])
    ub = ([[np.ones((1,N))],[None]])
    
#    print('g_mean_0')
#    print(G_mean_0.shape)

    f = lambda p_0: chisq(p_0, G_0, G_mean_0)
    

    #cons = ({}) constraints dict
    #opt = ({'maxiter' : 1e6}) options dict
    
    #result = scipy.optimize.minimize(f, p_0, args=(), jac=None, hess=None, hessp=None, bounds=((lb,ub)), constraints=(), tol=None, callback=None, options=None)
    #result = scipy.optimize.minimize(f, p_0, args=(), method=None , jac=None, hess=None, hessp=None, constraints=(), tol=None, callback=None, options=None)
    result = scipy.optimize.minimize(f, p_0, args=(), method=None, jac=None, bounds=None, constraints=(), tol=None, callback=None, options={'disp': True})
  #  result = _solve__(opt_problem={}, sens_type='FD', store_sol=True, disp_opts=False, store_hst=False, hot_start=False, sens_mode='', sens_step={}, *args, **kwargs)

    print('x')
    print(result.x)

    alpha = result.x[0] #1xN
    beta = result.x[1] ##1xN
    chi2 = result.fun
 
    g = (G_0 - alpha)/beta

    g,g_ignore = anchormean0to1(g)
    
    
    g = np.swapaxes(g,0,1)
    #toc
    return g, chi2

#Cost Function

def chisq(p,G,g_mean):
    
    print(p.shape)
    
    alpha = p[0]
    beta = p[1]
    chi2 = np.nanmean(np.nansum ((G - (alpha + beta * g_mean[:,None] ))**2))
    
    
    
    return chi2
    

inputData = 'gap_data_raw_dorsal_wt' 
inputFile = inputData + '.mat'
mat = scipy.io.loadmat(inputFile)

Hblist = [val for sublist in mat['data']['Hb'][0][0:102] for val in sublist]
Krlist = [val for sublist in mat['data']['Kr'][0][0:102] for val in sublist]
Gtlist = [val for sublist in mat['data']['Gt'][0][0:102] for val in sublist]
Knilist = [val for sublist in mat['data']['Kni'][0][0:102] for val in sublist]

gapDataRawDorsalArray = np.array([Hblist, Krlist, Gtlist, Knilist])

nGenes=gapDataRawDorsalArray.shape[0] # xy_alignment_script.m:26
nEmbryos= gapDataRawDorsalArray.shape[1] #  = 102
nEmbryosUsed = 3 #102 fly embryos total

assert(nEmbryos >= nEmbryosUsed) # xy_alignment_script.m:11

"""
Authors' description of the data: "Each profile vector is of length 1000
pixels, where 0 corresponds to the anterior (A) and 1000 to the posterior
(P) of the embryo."
To mitigate 'edge effects,' we extract only the middle 80% of the data to 
use throughout the analysis.
"""

nSamplePts = gapDataRawDorsalArray.shape[2] # = 1000
x=(np.arange(1 / nSamplePts,1,1 / nSamplePts)) 
idx_low=np.dot(round(nSamplePts),0.1)
idx_high=np.dot(round(nSamplePts),0.9)
x=np.arange(int(idx_low),int(idx_high))

## Set up raw data xy_alignment_script.m:28

idcs=np.sort(np.random.permutation(np.random.randint(1, nEmbryos,size= nEmbryosUsed)))# xy_alignment_script.m:30
Y_raw = gapDataRawDorsalArray.astype(float) # xy_alignment_script.m:32-36

Y_raw = Y_raw[:,idcs,:]
Y_raw = Y_raw[:,:,x]
padSize=round((nSamplePts*0.1))
Y_raw=np.pad(Y_raw,[(0,0),(0,0),(padSize,padSize)],mode='constant', constant_values=(np.nan))

"""
Up to this point, we have simply deleted data outside of the middle 80%
of the domain (replacing it with NaNs to act as a buffer to shift the
profiles to the left/right during X-alignment).
"""

#1 Align Y alone (optimize alpha and beta).

#Y_raw,Y_rawMean = anchormean0to1(Y_raw)
Y_yAlign1 = np.zeros(Y_raw.shape)

for iG in np.arange(0,nGenes):
    y=np.squeeze(Y_raw[iG,:,:])
    Y_yAlign1[iG,:,:],chi2Ignored = aligny(y)
    


#Plot the charts

x = np.arange(1000)

Y_anchor,Y_mean = anchormean0to1(Y_raw)

#for i in np.arange(1,102): 
for i in np.arange(1,3): 
   plt.plot(Y_anchor[0,i,:],color='grey', linewidth=0.5)
plt.plot(Y_mean[0,:],color='black')
plt.title('yalign before')
plt.show()


Y_mean = np.nanmean(Y_yAlign1,axis=1)
#for i in np.arange(1,102):    
for i in np.arange(1,3): 
   plt.plot(Y_yAlign1[0,i,:],color='grey', linewidth=0.5)
plt.plot(Y_mean[0,:],color='black')
plt.title('yalign after')
plt.show()