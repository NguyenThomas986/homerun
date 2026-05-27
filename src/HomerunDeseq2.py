import pandas as pd # gives excel like functions
import seaborn as sns # pretty plot package
from matplotlib import pyplot as plt # basic plotting package

import numpy as np # numerical functions
import os  # operation system functions
import sys # system functions
import scipy.stats as stat
import matplotlib.lines as mlines # for midlines
sys.setrecursionlimit(3000)

def quickDeseq2(workingPath):
    os.chdir(workingPath + '/analysis/')
    !mkdir DESeq2
    os.chdir(workingPath + '/analysis/DESeq2')

    raw_TSRs = pd.read_csv('../peakCalling/keyFiles/allTSSmerged_anoRaw.txt', sep='\t') # filename might need changed
    raw_TSRs.columns = [x.split('_')[1] if 'Total' in x  else x for x in raw_TSRs.columns.values]
    raw_TSRs = raw_TSRs.rename(columns={raw_TSRs.columns[0]: 'TSR_ID'})
    raw_TSRs = raw_TSRs.iloc[:,np.r_[0,19:29]]

    conditions = raw_TSRs.columns # probably doesn't work like this
    conditions_string = ' '.join(conditions)

    !getDiffExpression.pl ../peakCalling/keyFiles/allTSSmerged_anoRaw.txt $conditions_string -export ALL_conditions -AvsA > ALL_Conditions_DeSeq2.tsv