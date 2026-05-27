#!/usr/bin/env python3

# native imports that don't need installed
import os # operation system functions
import sys # system functions
import subprocess # subproccess system package
import time # native import capable of keeping track of time for me
import random as rand # random number generator and related packages
import re # regex functions

#start tracking the time the program takes to run
timerStart = time.perf_counter()

#try to import all packages, install if not
isPass = False
while isPass == False:
    try:
        # import bunch of programs needed
        import re
        import multiprocessing # multiprocess package
        from multiprocessing.connection import wait
        import glob
        import csv
        import gc
        import math
        from collections import Counter
        import argparse
        # modules for the umapping
        import numpy as np
        from sklearn.datasets import load_digits
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import StandardScaler
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
        # configure where plots are displayed
        #%matplotlib inline 
        import numpy as np # numerical functions
        import scipy.stats as stats # statistical engine
        import matplotlib.lines as mlines # for midlines
        sys.setrecursionlimit(3000)
        # exit the loop
        isPass = True
    # split up the error data and grab just the module name, then install that module with mamba
    except ModuleNotFoundError:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        splitObjects = str(exception_object).split("'")
        # now try to install it
        os.system("mamba install "+splitObjects[1])
        
# Parser ------------------------------------------------------------------------------------
# parse variables from user inputs into the items above
parser=argparse.ArgumentParser(prog = 'Homerun',
                               description = 'takes fastq files and creates TagDirs from them. Can then run them through QC and/or peakcalling, TSR grouping, and stability calling',
                               epilog = 'In order, this code runs through trimming --> mapping --> Tagdirs --> QC --> Peak Calling --> TSR Stability --> iTSS calling --> iTSS Stability --> Tagdir Stats --> Annotation Stats --> Stability Annotation --> TSR Stats --> Mapping Stats'
                               )
parser.add_argument("-p",
                    "--species",
                    help="<REQUIRED> The species to analyze. EXPERIMENTAL: List multiple species with spaces inbetween",
                    nargs='+',
                    required = True
                   )
parser.add_argument("-f",
                    "--fastqPath",
                    help="<REQUIRED> for tagdir! The file path to the fastq files to be used",
                    #required = True
                   )
parser.add_argument("-g",
                    "--genomePath",
                    help="<REQUIRED> for tagdir! The file path to the directory containing the .fa genome file to use",
                    #required = True
                   )
parser.add_argument('-w',
                    '--workingPath',
                    help='The parent file path to be used for the experiment. This should be the main folder for the experiment you are running where all files will be within it. If not specified, will assume directory where program was called',
                    default = os.getcwd()
                   )
parser.add_argument('-n',
                    '--ntag',
                    help = "<REQUIRED> ntag value to use for findcsRNApeaks",
                    type = int,
                    default = 7
                   )
parser.add_argument('-m',
                    '--mode',
                    help = "<REQUIRED> star or hisat",
                    required = True
                   )
parser.add_argument('-t',
                    '--step',
                    help = "tagdir, qc, or stat",
                    required = True
                   )

args = parser.parse_args()
species_lst = args.species     # this treats elements that are seperated by a space in the "species" section as an element of a list
fastqPath = args.fastqPath
genomePath = args.genomePath   # the directory containing the genome file that the species will be compared against
workingPath = args.workingPath # working path calls present working dir when not specified
ntag = args.ntag
mode = args.mode               # star or hisat
step = args.step               # tagdir, qc, or stat (secret: prep, all)

# add '/' at the end of the working path if needed
if not workingPath.endswith('/'): workingPath = workingPath + '/'

# create an error file
with open('errFile.out','w') as errFile:
    if step == 'tagdir' or step == 'prep':
        if mode is not None and fastqPath is not None and genomePath is not None:
            # build directories
            import HomerunDirBuilder
            HomerunDirBuilder.dirBuild(workingPath, fastqPath, genomePath, species_lst, mode)

            # trim and align
            import HomerunTrim
            HomerunTrim.align(workingPath, mode, species_lst)
    
    if step == 'tagdir':
        if mode is not None and fastqPath is not None and genomePath is not None:
            # TagDir creation
            import HomerunTagdir
            HomerunTagdir.createTagDirs(workingPath)
        else: print('Please make sure that you have specified the genome path, the fastq path, and whether the mode is STAR or HISAT.')
    
    # QC
    elif step == 'qc':
        import HomerunQC
        HomerunQC.qc(workingPath, species_lst, mode, errFile)

    # Statistical analysis
    elif step == 'stat':
        import HomerunTSS
        # Find peaks
        HomerunTSS.peakfinding(workingPath, species_lst)

        # TSR Stability Counts
        HomerunTSS.TSRstabilityCount(workingPath, errFile)

        # TSS from reads
        HomerunTSS.TSSfromReads(workingPath, errFile)

        # +TSS stats
        HomerunTSS.plusTSSstats(workingPath, errFile)

        # Tag Dir Stats
        HomerunTSS.tagDirStats(workingPath, species_lst, errFile)

        # TSR Annotation
        HomerunTSS.TSRannotation(workingPath, errFile)

        # TSR Stability Annotation
        HomerunTSS.TSRstabAnno(workingPath, errFile)

        # TSR detail summary
        HomerunTSS.TSRsummary(workingPath, errFile)

        # Summary Mapping Stats
        HomerunTSS.mappingStats(workingPath, errFile)

        # Summary Tag Dir Stats
        HomerunTSS.tagDirStats(workingPath, errFile)

errFile.close()
timerEnd = time.perf_counter()
print(f'Program elapsed time of {timerEnd - timerStart}')
