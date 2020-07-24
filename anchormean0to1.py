import numpy as np
import matplotlib.pyplot as plt


def anchormean0to1(Y):
    
    """Anchor the mean of a profile set from 0 to 1.
    
    Input:
      > Y: profiles to anchor (nPts x nSamples x nTargets).
    
    Output:
      > Y_anchor: profiles scaled such that the mean ranges from 0 to 1.
      > Y_mean: mean profile, min = 0 and max = 1.
    """
    
    Y_mean = np.nanmean(Y,axis=1)
    Y_anchor = ((Y - np.nanmin(Y_mean))/((np.nanmax(Y_mean)) - (np.nanmin(Y_mean))))    
    Y_mean = np.nanmean(Y_anchor,axis=1)



    #for i in np.arange(1,102): 
    #   plt.plot(Y_anchor[0,i,:],color='grey', linewidth=0.5)
       
    #plt.plot(Y_mean[0,:],color='black')
    #plt.title('Python out')
    #plt.show()

    return Y_anchor,Y_mean