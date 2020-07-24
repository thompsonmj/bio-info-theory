import scipy.io
import numpy as np

inputData = 'gap_data_raw_dorsal_wt' 
inputFile = inputData + '.mat'
data = scipy.io.loadmat(inputFile)

Hb_list = [val for sublist in data['data']['Hb'][0][0:102] for val in sublist]
Kr_list = [val for sublist in data['data']['Kr'][0][0:102] for val in sublist]
Gt_list = [val for sublist in data['data']['Gt'][0][0:102] for val in sublist]
Kni_list = [val for sublist in data['data']['Kni'][0][0:102] for val in sublist]

# 3D array with dimensions nSamplePoints x nEmbryos x nGenes
Y_raw = np.transpose(np.array([Hb_list, Kr_list, Gt_list, Kni_list]))

np.savetxt('Hb_raw_data.csv',Y_raw[:,:,0],delimiter=',',fmt='%f')

