import numpy as np
import matplotlib.pyplot as plt

from anchormean0to1 import anchormean0to1
from aligny import aligny,chisq

## Load data
Hb = np.genfromtxt('Hb_raw_data.csv', delimiter=",")

# Normalize
Hb_anchor, Hb_mean = anchormean0to1(Hb)
nSamplePoints, nEmbryos = Hb.shape

# Perform Y-alignment
Hb_yAlign = aligny(Hb_anchor)

# Plot results
for iE in range(0,nEmbryos):
    plt.plot(Hb_anchor[:,iE],color='grey',linewidth=0.5)
    
for iE in range(0,nEmbryos):
    plt.plot(Hb_yAlign[:,iE],color='black',linewidth=0.5)


