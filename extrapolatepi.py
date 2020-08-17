import numpy as np

nBinSizes = np.size(MI_est[0].binSizes)
nSubsamps = np.size(MI_est[0].subSamps)
nTrials = np.size(MI_est)

m= round(MI_est[0].subSamps*nEmbryos)

invSamps = 1/m
invSamps = invSamps*3
invBinsSq = 1/(MI_est[0].binSizes**2)
invBinsSq = invBinsSq.np.transpose()
xRegress = np.linspace(0, max(invSamps), 0.001)
intercepts_infiniteData = np.zeros(nBinSizes,nTrials)
intercepts_infiniteBins = np.zeros(1,nTrials)

themean = np.zeros(nBinSizes,nSubsamps)

for iT in range(0, nTrials):
    for iB in range(0, iBinSizes):
        coeffs_thisB = np.polyfit(invSamps, MI_est[iT].naiveEst_means[iB,:],1)
        thisIntercept = coeffs_thisB[2]
        intercepts_infiniteData[iB,iT] = thisIntercept
        
    coeffs_thisT = np.polyfit(invBinsSq, intercepts_infiniteData[:,iT],1)
    thisIntercept = coeffs_thisT[2]
    intercepts_infiniteBins[iT] = thisIntercept
    
intercepts_infiniteData_means = np.mean(intercepts_infiniteData,2)
intercepts_infiniteData_stds = np.std(intercepts_infiniteData,0,2)

naiveEst_trialMeans = np.mean(np.concatenate(3,MI_est[:].naiveEst_means),3)
naiveEst_trialStds = np.std(np.concatenate(3,MI_est[:].naiveEst_means),0,3)

#naiveEst_trialMeans = np.mean(cat(3,MI_est(:).naiveEst_means),3)
#naiveEst_trialStds = np.std(cat(3,MI_est(:).naiveEst_means),0,3)


#plot results[]
for iB in range(0, iBinSizes):
    MI = np.zeros(nBinSizes,nSubsamps,nTrials)
    
    for iT in range(0, iTrials):
        MI[:,:,iT] = MI_est[iT].naiveEst_means
        
    MI_means = np.mean(MI,3)
    MI_stds = np.std(MI,0,3)   
    
for iB in range(0, nBinSizes):
    coeffs_thisB = np.polyfit(invSamps,MI_means[iB,:],1)
    f = np.polyval(coeffs_thisB, 0)
    intercepts_infiniteData[iB] = intercept
    
    meanIntercepts = np.mean(intercepts_infiniteData, 2)
    stdIntercepts = np.std(intercepts_infiniteData, 1, 2)
    
    invBinSq = 1/(binCounts**2)
    
    p2 = np.polyfit(invBinSq.conj().transpose(), meanIntercepts, 1)
    f2 = np.polyval(p2, np.linspace(0, max(invBinSq), 0.001))
 
    intercept2 = np.polyval(p2,0)

    MI_ext = intercept2
    print(MI_ext)

    
