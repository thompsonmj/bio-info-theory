import numpy as np
import random
import math

def midirectestimate2(yData,nTrials,nBoots,binCounts,subSamps):
    
    nBinsX = yData.shape[0]
    nE = yData.shape[1]
    nG = yData.shape[2]

    x=np.arange(1/nBinsX, 1, 1/nBinsX) # start: 1/nBinsX, stop: 1, step: 1/nBinsX

    assert (nBinsX == np.size(x))

    # nX = 100; % For pSmad data with such small sample size.

    m = round(subSamps*nE)
    nSubsamplesets = np.size(m)

    replacement = False

    embryoIdcs = [np.arange(1, nE, 1)]

    nBinSizes = np.size(binCounts)
    
    MI_est = {}

    for iTrial in range(0,nTrials):
        
        MI_est[iTrial] = {'binSizes':binCounts,'subSamps':subSamps,'nBoots':nBoots,'naiveEst_means': (np.zeros(nBinSizes,nSubsamplesets)),'naiveEst_stds':(np.zeros(nBinSizes,nSubsamplesets))}
        
    
        if nTrials > 1:
            print("Beginning trial " + str(iTrial) + "/" + str(nBinSizes))
            #tic
            
            iBinSizeIdx =0
            
            for iBinSize in range(0,nBinSizes):
                if nTrials == 1:
                    print("Beginning bin " + str(iBinSize) + "/" + str(nBinSizes))
                    #tic
            
            nBinsG = binCounts[iBinSize]
        
            nBinsX = nBinsG
            pxBin = 1/nBinsX
        
            iBinSizeIdx = iBinSizeIdx + 1
            
            for iSubsampleIdx in range(0, nSubsamplesets):
                mSubsamps = m[iSubsampleIdx]
                pg = np.zeros(nBinsG, nBinsG)
                pgxJ = np.zeros(nBinsX, nBinsG,nBinsG)
                MIEq24 = np.zeros(nBoots, 1)
                for iBoot in range(0, nBoots):
                    subIdx = random.choices(embryoIdcs, weights=None, k=mSubsamps)
                    sub1 = yData[:,subIdx, 1]
                    sub2 = yData[:,subIdx, 2]
                    
                    sub1 = np.reshape(sub1, np.size(sub1))
                    sub2 = np.reshape(sub2, np.size(sub2))
                    
                    xMat = np.tile(x.conj().transpose(), (1, mSubsamps))
                    xColumn = np.reshape(xMat, np.size(xMat))
                    
                    #Distributions
                    #pgxJ = histcountsn([xColumn, sub],[nBinsX,nBinsG], ...
                    #'Normalization','probability'); %joint probability distribution of G and X
                    #Use new histcounts function
                    
                    
                    pgxJ =pgxJ + eps
                    pgxJ_perm = np.transpose(pgxJ(2,3,1))
                    
                    pg = np.sum(pgxJ,1).squeeze()
                    px = np.tile(pxBin, (nBinsX, 1))
                    
                    pg_given_x = pgxJ / px
                    pg_given_x_perm = np.transpose(pg_given_x,(2,3,1))
                    
                    #Eq. 7
                    logValue = pg_given_x_perm  * math.log2((pg_given_x_perm /pg))
                    
                    MIEq7_temp = np.sum(px * np.sum(logValue,3))
                
                    MIEq7[iBoot] = MIEq7_temp
                    
                    #Eq. 24
                    MIEq24_temp = np.sum(pgxJ_perm * math.log2(pgxJ_perm / (px * pg)))
                
                    MIEq24[iBoot] = MIEq24_temp
                    
                    #Validation
                    assert(abs(MIEq7_temp - MIEq24_temp) < 1e-10)
                
            
            if nTrials ==1:
                print(str(round(100*iBinSize/nBinSizes,2)))
                
        if nTrials>1:
            print(str(100*iTrial/nTrials,2))
            
    #save file
                    
                    
        
    
    return MI_est