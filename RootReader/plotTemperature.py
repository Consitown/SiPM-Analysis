import numpy as np
import matplotlib.pyplot as plt


import math
import csv
import pandas as pd
import plotly.express as px

import matplotlib as mpl

#matplotlib.verbose.level = 'debug-annoying'
mpl.rcParams.update({'font.size': 16})
import sys

temperatureFile = "/mnt/d/RUNDATA_UNI/temperature/Temperature Detector Lab-data-as-seriestocolumns-2021-04-23 10 24 45.csv" #sys.argv[1]
#saveFolder = sys.argv[2]

df = pd.read_csv(temperatureFile)
df.head()


