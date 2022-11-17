import numpy as np
import matplotlib.pyplot as plt


import math


import matplotlib as mpl

#matplotlib.verbose.level = 'debug-annoying'
mpl.rcParams.update({'font.size': 16})
import sys

meansFile = sys.argv[1]
saveFolder = sys.argv[2]

# print(meansFile)
# print(saveFolder)

file = open(meansFile)
integralMeans = np.loadtxt(meansFile, skiprows=1)[:, 1]
integralMeanErr = np.loadtxt(meansFile, skiprows=1)[:, 2]
angles = np.loadtxt(meansFile, skiprows=1)[:, 0]

integralMeansSum = np.sum(integralMeans)

plt.errorbar(angles, integralMeans, yerr=integralMeanErr)
plt.xlabel("Angle assigned to Channel (°)")
plt.ylabel("Integral Mean normalised to Means' Sum")
plt.xticks(angles)

plt.savefig(saveFolder + "/propagatedIntegralMeans.png", bbox_inches="tight", dpi=300)





