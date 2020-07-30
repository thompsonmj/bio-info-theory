import numpy as np
import matplotlib.pyplot as plt

import csv

##################### Import Data, Clean Nan's ################################
def sigma_x(filename):
    Y_in = []
    position = []
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
        for row in reader: # each row is a list
            Y_in.append(row)
    Y_in = np.array(Y_in)
    
    Y_in= Y_in[~np.isnan(Y_in).any(axis=1)]
    
    for i in range(0,len(Y_in)):
        placeHold = i / len(Y_in)
        position.append(placeHold)
    
    ##################### Row mean of the Array ###############################
    mean = np.mean(Y_in,axis=1);
    normalizedMean = [];
    minMean = min(mean);
    maxMean = max(mean);
    #print("Length of Mean Array: "+ str(len(mean)))
    for i in range(0,len(mean)):
        placeHold = (mean[i] - minMean) / maxMean
        normalizedMean.append(placeHold)
        
    #################### Calculate the Difference, Create Array ##################
    dgdx = []
    posErrorX = []
    count = 0
    for i in range(0,len(Y_in)-1):
        count = count
        d = (mean[i+1] - mean[i])
        d = abs(d**-1)
        dgdx.append(d)
        posErrorX.append(count/(len(mean)-1))
        count = count + 1
        
    #print("Length of dgdx: " + str(len(posErrorX)))
    # should be N-1, where N is the len of mean 
    #print(posErrorX)
        
    ##################### Calculate Standard Deviation ############################
    sigG = []
    SD = np.std(Y_in, axis= 1)
    for i in range(0,len(Y_in)-1):
        sigG.append(SD[i])
    
    #print("Length of sigG: " + str(len(sigG)))
    ##################### Calculate SigmaX ########################################
    sigmaX = np.multiply(sigG,dgdx) / len(dgdx);
    
    plt.figure(1)
    plt.plot(position,normalizedMean)
    
    plt.figure(2)
    plt.plot(posErrorX,sigmaX)
    plt.ylim(0,0.1)
    
sigma_x("Hb_raw_data.csv")