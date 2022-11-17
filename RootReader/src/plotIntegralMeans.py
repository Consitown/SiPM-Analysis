import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
from scipy.odr import Model, ODR, Data, RealData
import copy
import matplotlib as mpl
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{nicefrac} \usepackage[separate-uncertainty]{siunitx} \usepackage{physics}')
#matplotlib.verbose.level = 'debug-annoying'
mpl.rcParams.update({'font.size': 16})
import sys

meansFile = "/mnt/d/Programme/RootAnalysis/RootAnalysis/integralAnalysis/0_pos_no_ang_no/integralMeans.txt" #sys.argv[1]

file = open(meansFile)
integralMeans = np.loadtxt(meansFile, skiprows=1)[:, 0]
integralMeanErr = np.loadtxt(meansFile, skiprows=1)[:, 1]
StdDev = np.loadtxt(meansFile, skiprows=1)[:, 2]
StdDevErr = np.loadtxt(meansFile, skiprows=1)[:, 3]

print(integralMeans)





