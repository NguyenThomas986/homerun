import os # operation system functions
import sys
import time
import re
import subprocess
import multiprocessing # multiprocess package
from multiprocessing.connection import wait

# trim / align script

def align(workingPath, mode, species_lst):
    iTSSPath = f'{workingPath}files/iTSS/'
    samLst = []
    dataPath = f'{workingPath}data/'
    filesPath = f'{workingPath}files/'
    TSRPath = f'{workingPath}files/TSR/'
    os.chdir(dataPath)
    print('Beginning trimming and alignment...')
    for species in species_lst:
        speciesPath = f'{dataPath}{species}/'
        os.chdir(f'{speciesPath}fastq')
        tagPath = f'{speciesPath}tagDirs/'
        samPath = os.getcwd()
        # loop through the sequencing type folders and run accordingly
        for seqType in os.listdir():
            if os.path.isdir(f'{samPath}/{seqType}') != True:
                continue
            elif os.path.isdir(f'{samPath}/{seqType}') == True:
                if (seqType == 'csRNA') or (seqType == 'sRNA'):
                    os.chdir(f'{samPath}/{seqType}')
                    
                    # check if the fastqs are unzipped, zip them for homer trim if not
                    zipLst = []
                    for eachFastq in os.listdir():
                        # check if the fastqs are zipped
                        if eachFastq.endswith('.fastq'):
                            zipLst.append(eachFastq)

                    # multiprocess the zipping
                    if len(zipLst) > 0:
                        procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                        procPool.map(gzipFunc, zipLst)
                        procPool.close() # close the pool
                        procPool.join()

                    # loop a second time time to find the correct file names and append them to a list
                    fastqLst = []
                    for eachFastq in os.listdir():
                        if eachFastq.endswith('.fastq.gz'):
                            fastqLst.append(eachFastq)
                        
                    # trim ... with a multiprocessing pool (this is the last block where hisat / star is the same)
                    if len(fastqLst) > 0:
                        procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                        procPool.map(homerTrimFunc, fastqLst)
                        procPool.close() # close the pool
                        procPool.join()
                        
                    # this is where HISAT and STAR processes diverge
                    if mode == 'star':
                        # don't multiprocess the star function as it eats too much memory and will almost always fail
                        for eachTrim in os.listdir():
                            if eachTrim.endswith('.fastq.gz.trimmed'):
                                if '_R1' in eachTrim:
                                    starFunc(eachTrim,species,workingPath)
                    
                    elif mode == 'hisat':
                        trimLst = []
                        for eachTrim in os.listdir():
                            if eachTrim.endswith('.fastq.gz.trimmed'):
                                trimLst.append([eachTrim,species, workingPath])
                        # run the hisat program through multiprocessing
                        if len(trimLst) > 0:
                            procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                            #procPool.starmap(hisatAlign, trimLst, workingPath)
                            procPool.starmap(hisatAlign, trimLst)
                            procPool.close() # close the pool
                            procPool.join()
                            
                # STAR totalRNA
                elif seqType == 'totalRNA' and mode == 'star':
                    os.chdir(f'{samPath}/{seqType}')
                    # start by unzipping files
                    unzipLst = []
                    for eachFile in os.listdir():
                        if eachFile.endswith('fastq.gz'):
                            unzipLst.append(eachFile) # append to a list to be unzipped  
                    if len(unzipLst) > 0:
                        procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                        procPool.map(gunzipFunc, unzipLst)
                        procPool.close() # close pool
                        procPool.join()
                    # run Star aligner on the totalRNA fastq files... one... at... a... time...
                    for eachTrim in os.listdir():
                        if eachTrim.endswith('.fastq'):
                            if '_R1' in eachTrim:
                                starFunc(eachTrim,species,workingPath)
                                
                # HISAT totalRNA
                elif seqType == 'totalRNA' and mode == 'hisat':
                    os.chdir(f'{samPath}/{seqType}')
                    unzipLst = [] # unzip the files that need it
                    for eachFastq in os.listdir():
                        # check if fastqs are zipped
                        if eachFastq.endswith('.fastq.gz'):
                            unzipLst.append(eachFastq)
                    # multiprocess zipping
                    if len(unzipLst) > 0:
                        procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                        procPool.map(gunzipFunc, unzipLst)
                        procPool.close() # close pool
                        procPool.join()
                    skewerLst = []
                    for eachFile in os.listdir():
                        pairMatch = False
                        # select only samples containing the R1
                        if ('_R1' in eachFile) and ('trimmed' not in eachFile):
                            # split out the SRR number and replace R1 with R2
                            SRRsplit = eachFile.split('_SRR')[-1]
                            SRR = '_SRR' + SRRsplit.split('_')[0]
                            sampleName = eachFile.replace(SRR, '')
                            sample2Name = sampleName.replace('_R1','_R2')
                            # loop through the file a second time stripping out the SRR number like before and check if the two match
                            for otherFile in os.listdir():
                                if '_R2' in otherFile:
                                    SRRsplit = otherFile.split('_SRR')[-1]
                                    SRR = '_SRR' + SRRsplit.split('_')[0]
                                    otherName = otherFile.replace(SRR, '')
                                    other1Name = otherName.replace('_R2','_R1')
                                    # check to see if they match
                                    if other1Name == sampleName:
                                        print(f'pairs found:\n{eachFile}\t{otherFile}\n')
                                        fastqPath = samPath
                                        trimmed_fastq_output = fastqPath + f'/{seqType}/'+ eachFile.split('_SRR')[0]
                                        trim = f'time -p skewer -m mp {eachFile} {otherFile} -t 40 -o {trimmed_fastq_output}'
                                        skewerLst.append(trim)
                                        pairMatch = True
                                        break
                            if pairMatch == False:
                                print(f'found no pairs for:\n{eachFile}\n')
                                fastqPath = samPath
                                trimmed_fastq_output = fastqPath + f'/{seqType}/' + eachFile.split('_SRR')[0]
                                trim = f'time -p skewer -m mp {eachFile} -t 40 -o {trimmed_fastq_output}'
                                skewerLst.append(trim)
                    # now run the multi
                    if len(skewerLst) > 0:
                        procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                        print('\nRunnning Trim Commands:')
                        procPool.map(genericFunc, skewerLst)
                        procPool.close() # close the pool
                        procPool.join()
                    # after trimming is complete, loop through files looking for the trimmed files and do the same
                    # except with the hisat aligner
                    trimLst1 = []
                    trimLst2 = []
                    for eachTrim in os.listdir():
                        if eachTrim.endswith('-trimmed-pair1.fastq'):
                            trimLst1.append([eachTrim,species,workingPath])
                        elif eachTrim.endswith('-trimmed.fastq'):
                            trimLst2.append([eachTrim,species,workingPath])
                    # run the hisat program through multiprocessing
                    if len(trimLst1) > 0:
                        print('\nRunning Align Commands:')
                        procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                        procPool.starmap(rnaHisatAlign, trimLst1)
                        procPool.close() # close the pool
                        procPool.join()
                    if len(trimLst2) > 0:
                        print('\nRunning Align Commands:')
                        procPool = multiprocessing.Pool(int(os.getenv('SLURM_CPUS_PER_TASK')))
                        procPool.starmap(rnaSingleHisat, trimLst2)
                        procPool.close() # close the pool
                        procPool.join()
                # run the Chipseq files too
                elif seqType == 'ChIPseq' and mode == 'hisat':
                    os.chdir(f'{samPath}/{seqType}')
                    skewerLst = []
                    # include a way of dropping out the ENCF identifier to find the other file
                    for eachFile in os.listdir():
                        if eachFile.endswith('R1.fastq'):
                            encodeID = eachFile.split('_')[4]
                            subread1 = eachFile.split(encodeID)[0]+eachFile.split(encodeID)[-1]
                            subread2 = eachFile.split('R1.fastq')[0] + 'R2.fastq'
                            # loop through again finding the one that matches
                            for otherFile in os.listdir():
                                encodeID2 = otherFile.split('_')[4]
                                posRead2 = otherFile.split(encodeID2)[0]+otherFile.split(encodeID2)[-1]
                                if subread1.split('R1.fastq')[0] == posRead2.split('R2.fastq')[0]:
                                    read1 = eachFile
                                    read2 = otherFile
                                    trimmed_fastq_output = os.getcwd() +'/'+ eachFile.split('_R1.fastq')[0]
                                    trim = f'time -p skewer -m mp {read1} {read2} -t 40 -o {trimmed_fastq_output}'
                                    skewerLst.append(trim)
                                else:
                                    continue
    # loop through the files again to grab the sam files and output them into the fastq parent folder
    if mode == 'hisat':
        os.chdir(dataPath)
        for species in species_lst:
            speciesPath = f'{dataPath}{species}/'
            os.chdir(f'{speciesPath}fastq')
            tagPath = f'{speciesPath}tagDirs/'
            samPath = os.getcwd()
            for seqType in os.listdir():
                os.chdir(f'{speciesPath}fastq')
                if os.path.isdir(f'{samPath}/{seqType}') != True:
                    continue
                elif os.path.isdir(f'{samPath}/{seqType}') == True:
                    os.chdir(seqType)
                    for eachSam in os.listdir():
                        if eachSam.endswith('.sam'):
                            mvCmd = f'mv {eachSam} ../'
                            os.system(mvCmd)
                        elif eachSam.endswith('mappingstats.txt'):
                            mvCmd = f'mv {workingPath}files/mappingStats/'
                            os.system(mvCmd)
    elif mode == 'star':
        os.chdir(dataPath)
        for species in species_lst:
            speciesPath = f'{dataPath}{species}/'
            os.chdir(f'{speciesPath}fastq')
            tagPath = f'{speciesPath}tagDirs/'
            samPath = os.getcwd()
            for seqType in os.listdir():
                os.chdir(samPath)
                if not os.path.isdir(f'{samPath}/{seqType}'):
                    if seqType.endswith('Log.final.out'):
                        mvCmd = f'mv {seqType} {workingPath}files/mappingStats/{seqType.split(".Log.final")[0]}.mappingstats.txt'
                        os.system(mvCmd)
                elif os.path.isdir(f'{samPath}/{seqType}') == True:
                    os.chdir(seqType)
                    for eachSam in os.listdir():
                        if eachSam.endswith('.sam'):
                            mvCmd = f'mv {eachSam} ../'
                            os.system(mvCmd)
                        elif eachSam.endswith('Log.final.out'):
                            mvCmd = f'mv {eachSam} {workingPath}files/mappingStats/{eachSam.split(".Log.final")[0]}.mappingstats.txt'
                            os.system(mvCmd)
    print("Alignment complete!")

# functions used in the big one above

def homerTrimFunc(fastqFile):
    # TO DO: re-implement trimming
    # cmd = f'homerTools trim -3 TATAAA -mis 2 -minMatchLength 4 -min 20 {fastqFile}'
    cmd = f'homerTools trim -mis 2 -minMatchLength 4 -min 20 {fastqFile}'
    print(cmd)
    process = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)
    process.communicate()

def starFunc(trimFile, species, workingPath):
    outName = '../' + trimFile.split('.fastq')[0]
    trimFile2 = trimFile.replace('_R1','_R2')
    starIndex = f'{workingPath}genomes/{species}/STARIndex'
    if os.path.exists(trimFile2):
        starCmd = f'STAR --genomeDir {starIndex} --runThreadN 20 --readFilesIn {trimFile} {trimFile2} --outFileNamePrefix {outName}. --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterMultimapNmax 10000 --limitOutSAMoneReadBytes 10000000'
    elif not os.path.exists(trimFile2):
        starCmd = f'STAR --genomeDir {starIndex} --runThreadN 20 --readFilesIn {trimFile} --outFileNamePrefix {outName}. --genomeLoad NoSharedMemory --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterMultimapNmax 10000 --limitOutSAMoneReadBytes 10000000'
    print(starCmd)
    process = subprocess.Popen([starCmd], shell=True, stdout=subprocess.PIPE)
    process.communicate()
    
def hisatAlign(csRNAtrim, species, workingPath):
    csRNApath = os.getcwd()
    output_sam = csRNAtrim.split('.fastq.gz.trimmed')[0] + '.sam'
    mapping_stats_path = f'{workingPath}files/mappingStats/'
    mapping_file = mapping_stats_path + csRNAtrim.split('.fastq.gz.trimmed')[0] + '_mappingstats.txt'
    hisat2_index =  f'{workingPath}genomes/' + species + '/Hisat2/index'
    map_SEreads = f'hisat2 -p 4 --rna-strandness RF --dta -x {hisat2_index} -U {csRNApath}/{csRNAtrim} -S ../{output_sam} 2> {mapping_file}'
    print(map_SEreads)
    process = subprocess.Popen([map_SEreads], shell=True, stdout=subprocess.PIPE)
    process.communicate()

def rnaHisatAlign(RNAtrim1, species, workingPath):
    output_sam = '../' + RNAtrim1.split('-trimmed-pair1.fastq')[0] + '.sam'
    RNAtrim2 = RNAtrim1.split('-trimmed-pair1.fastq')[0] + '-trimmed-pair2.fastq'
    mapping_stats_path = f'{workingPath}files/mappingStats/'
    mapping_file = mapping_stats_path + RNAtrim1.split('-trimmed-pair1.fastq')[0] + '_mappingstats.txt'
    hisat2_index = f'{workingPath}genomes/' + species + '/Hisat2/index'
    map_SEreads = map_RNA = f'hisat2 -p 4 --rna-strandness RF --dta -x {hisat2_index} -1 {RNAtrim1} -2 {RNAtrim2} -S {output_sam} 2> {mapping_file}'
    print(map_SEreads)
    process = subprocess.Popen([map_SEreads], shell=True, stdout=subprocess.PIPE)
    process.communicate()
    
def rnaSingleHisat(RNAtrim1,species, workingPath):
    output_sam = '../' + RNAtrim1.split('-trimmed.fastq')[0] + '.sam'
    mapping_stats_path = f'{workingPath}files/mappingStats/'
    mapping_file = mapping_stats_path + RNAtrim1.split('-trimmed.fastq')[0] + '_mappingstats.txt'
    hisat2_index =  f'{workingPath}genomes/' + species + '/Hisat2/index'
    map_SEreads = map_RNA = f'hisat2 -p 4 --rna-strandness RF --dta -x {hisat2_index} -U {RNAtrim1} -S {output_sam} 2> {mapping_file}'
    print(map_SEreads+'\n')
    process = subprocess.Popen([map_SEreads], shell=True, stdout=subprocess.PIPE)
    process.communicate()
    
def gzipFunc(file):
    gzipCmd = f'gzip {file}'
    process = subprocess.Popen([gzipCmd],shell=True, stdout=subprocess.PIPE)
    process.communicate()

def gunzipFunc(file):
    gunzipCmd = f'gunzip {file}'
    process = subprocess.Popen([gunzipCmd],shell=True, stdout=subprocess.PIPE)
    process.communicate()
    
def genericFunc(cmd):
    print(cmd)
    proc = subprocess.Popen([cmd],shell=True, stdout=subprocess.PIPE)
    proc.communicate()