import numpy as np
import math as ma

import csv
import pandas as pd 

Y_in = []

with open("BIOTHEORY.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        Y_in.append(row)

##################### Row mean of the dataframe ###############################
mean = np.mean(Y_in,axis=1)
print(len(mean))

# ################### Calculate the Difference, Create Array ##################
dgdx = []
for i in range(0,len(Y_in)-1):
    d = (mean[i+1] - mean[i])
    d = abs(d**-1)
    dgdx.append(d)
#print((dgdx))
# should be N-1, where N is the len of mean 
    
##################### Calculate Standard Deviation ############################
sigG = []
SD = np.std(Y_in, axis= 1)
for i in range(0,len(Y_in)-1):
    sigG.append(SD[i])

#print(stDev)
##################### Calculate SigmaX ########################################
sigmaX = np.multiply(sigG,dgdx) / len(dgdx);
