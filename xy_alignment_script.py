import scipy.io
import matplotlib.pyplot as plt
import numpy as np

# Is there a cleaner way to use these functions?
from anchormean0to1 import anchormean0to1
from aligny import aligny,chisq
from alignxy import alignxy#,chisq

inputData = 'gap_data_raw_dorsal_wt' 
inputFile = inputData + '.mat'
mat = scipy.io.loadmat(inputFile)

Hblist = [val for sublist in mat['data']['Hb'][0][0:102] for val in sublist]
Krlist = [val for sublist in mat['data']['Kr'][0][0:102] for val in sublist]
Gtlist = [val for sublist in mat['data']['Gt'][0][0:102] for val in sublist]
Knilist = [val for sublist in mat['data']['Kni'][0][0:102] for val in sublist]

# 3D array with dimensions nSamplePoints x nEmbryos x nGenes
gapDataRawDorsalArray = np.transpose(np.array([Hblist, Krlist, Gtlist, Knilist]))

nGenes=gapDataRawDorsalArray.shape[2] # xy_alignment_script.m:26
nEmbryos= gapDataRawDorsalArray.shape[1] #  = 102
nEmbryosUsed = 102 #102 fly embryos total

assert(nEmbryos >= nEmbryosUsed) # xy_alignment_script.m:11

"""
Authors' description of the data: "Each profile vector is of length 1000
pixels, where 0 corresponds to the anterior (A) and 1000 to the posterior
(P) of the embryo."
To mitigate 'edge effects,' we extract only the middle 80% of the data to 
use throughout the analysis.
"""

nSamplePts = gapDataRawDorsalArray.shape[0] # = 1000
x=(np.arange(1 / nSamplePts,1,1 / nSamplePts)) 
idx_low=np.dot(round(nSamplePts),0.1)
idx_high=np.dot(round(nSamplePts),0.9)
x=np.arange(int(idx_low),int(idx_high))

## Set up raw data xy_alignment_script.m:28

idcs=np.sort(np.random.permutation(np.random.randint(1, nEmbryos,size= nEmbryosUsed)))# xy_alignment_script.m:30
Y_raw = gapDataRawDorsalArray.astype(float) # xy_alignment_script.m:32-36

Y_raw = Y_raw[:,idcs,:]
Y_raw = Y_raw[x,:,:]
padSize=round((nSamplePts*0.1))
Y_raw=np.pad(Y_raw,[(padSize,padSize),(0,0),(0,0)],mode='constant', constant_values=(np.nan))

"""
display anchormean0to1 for Hb gene
"""

x = np.arange(1000)

Y_mean = np.nanmean(Y_raw,axis=1) #"in" plot for Hb
#for i in np.arange(1,102): 
#    plt.plot(gapDataRawDorsalArray[:,i,0],color='grey', linewidth=0.5)
#plt.plot(Y_mean[:,0],color='black')
#plt.title('Python in')
#plt.show()

Y_anchor,Y_mean = anchormean0to1(Y_raw)

#for i in np.arange(1,102): 
#   plt.plot(Y_anchor[:,i,0],color='grey', linewidth=0.5)
   
#plt.plot(Y_mean[:,0],color='black')
#plt.title('Python out')
#plt.show()


# =============================================================================
# ## Initial Y alignment
# Y_yAlign1 = np.zeros(Y_anchor.shape)
# for iG in np.arange(0,nGenes):
#     y = np.squeeze(Y_anchor[:,:,iG])
#     #plt.plot(y,color='grey',linewidth=0.5)
#     #plt.show()
#     
#     y_align_temp = aligny(y)
#     Y_yAlign1[:,:,iG] = y_align_temp
#     
#     #plt.plot(y_align_temp,color='black',linewidth=0.5)
#     #plt.show()
# =============================================================================
    
## XY joint alignment
Y_xyAlign = np.zeros(Y_anchor.shape)

Y_xyAlign = alignxy(Y_anchor) # Pass all genes simultaneously
