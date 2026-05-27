import os # operation system functions
import subprocess
import re
import multiprocessing # multiprocess package
from multiprocessing.connection import wait

# creates tagdirs for data available in workingPath
def createTagDirs(workingPath):
    print('Beginning tagdir creation...')
    dataPath = f'{workingPath}data/'
    os.chdir(dataPath) # workingPath/data/
    for species in os.listdir():
        speciesPath = f'{dataPath}{species}/'
        os.chdir(f'{speciesPath}fastq') # workingPath/data/{species}/fastq/
        tagPath = f'{speciesPath}tagDirs/'
        samPath = os.getcwd()
        # make the two neccessary files: samNames.txt, sampleInfo.txt
        with open('samNames.txt','w') as samNames:
            for eachFile in os.listdir():
                if eachFile.endswith('.sam'):
                    samNames.write(f'{eachFile}\n')
                else:
                    continue
        samNames.close()
        # make sampleInfo.txt - this file consists of <full tagdir name and path>\t<full sam name and path>
        with open('sampleInfo.txt','w') as sampleInfo:
            for eachFile in os.listdir():
                if eachFile.endswith('.sam'):
                    samLine = f'{samPath}/{eachFile}'
                    tagName = eachFile.split('-r')[0]+'-r' + eachFile.split('-r')[-1][0]
                    tagLine = f'{tagPath}{tagName}'
                    writeLine = f'{tagLine}\t{samLine}\n'
                    sampleInfo.write(writeLine)
        sampleInfo.close()
        # run the batch build tagdir command using Popen
        genomefa =  f'{workingPath}genomes/' + species + '/*.fa'
        tagDirCmd = f'batchMakeTagDirectory.pl sampleInfo.txt -cpu 2 -genome {genomefa} -omitSN -checkGC -fragLength 150 -r'
        proc = subprocess.Popen([tagDirCmd],shell=True).wait()
    print('Tag dir process complete')