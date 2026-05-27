#!/usr/bin/env python3
import os # operation system functions
import sys # system functions
import csv
import shutil # for file copying

def dirBuild(workingPath, fastqPath, genomePath, speciesList, mode):
    os.chdir(workingPath)
    print('Building directories...')
    
    # create highest level folders
    dirNames = ["analysis", "badButKeep", "data", "files", "genomes", "GEOsubmit", "scripts"]
    for subDir in dirNames:
        os.makedirs(subDir, exist_ok=True)
    
    # create species folders
    os.chdir("data/")
    for species in speciesList:
        if speciesList.index(species) != len(speciesList):
            os.makedirs(species, exist_ok=True)
            print("data/" + species + "/ created")
            
    # create species subfolders
    speciesSubdirNames = ["bam", "bedgraphs", "fastq", "igv", "tagDirs"]
    fastqSubdirNames = ["csRNA", "totalRNA", "sRNA"]
    for species in os.listdir():
        for subDir in speciesSubdirNames:
            os.makedirs(os.path.join(workingPath, "data/", species, subDir), exist_ok=True)
            print(os.path.join("data/", species, subDir) + " created")
        # populate the fastq folder
        for subDir in fastqSubdirNames:
            os.makedirs(os.path.join(workingPath, "data/", species, "fastq/", subDir), exist_ok=True)
            print(os.path.join("data/", species, "fastq/", subDir) + " created")
    
    # create files subfolders
    os.chdir(workingPath)
    filesSubdirNames = ["mappingStats", "QC", "TSR", "iTSS"]
    for subDir in filesSubdirNames:
        os.makedirs(os.path.join(workingPath, "files/", subDir), exist_ok=True)
        
    # genome folder
    newGenomePath = workingPath+'genomes/'
    genomeFa = ""
    os.chdir(newGenomePath)
    # make genome/species/ folder
    for species in speciesList:
        os.makedirs(species, exist_ok=True)
    # copy files from genome/species/
    os.chdir(workingPath)
    os.chdir(genomePath)
    for file in os.listdir():
        os.chdir(workingPath)
        os.chdir(genomePath)
        if os.path.isfile(file):
            shutil.copy(file, newGenomePath+species)
            print("Copied from genome folder: " + file)
    if mode == "star":
        # copy files from genome/species/STARIndex/
        for species in speciesList:
            os.chdir(os.path.join(newGenomePath,species))
            os.makedirs("STARIndex", exist_ok=True)
            os.chdir(workingPath) # because the genome path entered as an argument is relative
            os.chdir(genomePath+'/STARIndex')
            # copy original STARIndex to the new STARindex
            for file in os.listdir():
                shutil.copy(file, newGenomePath+species+'/STARIndex/')
                print("Copied from STARIndex folder: " + file)
        print("Copied " + genomePath + " to " + newGenomePath+species+'/')
    # TO DO: handle multiple species?
    
    # copy raw data to workingPath/data/{species}/fastq/
    os.chdir(workingPath) # because fastqPath is relative
    for species in speciesList:
        for fastqFile in os.listdir(fastqPath):
            if "csRNA" in fastqFile and fastqFile.endswith("fastq.gz"): shutil.copy(os.path.join(fastqPath, fastqFile), os.path.join(workingPath, "data/", species, "fastq/csRNA/"))
            elif "sRNA" in fastqFile and fastqFile.endswith("fastq.gz"): shutil.copy(os.path.join(fastqPath, fastqFile), os.path.join(workingPath, "data/", species, "fastq/sRNA/"))
            elif "RNA" in fastqFile and fastqFile.endswith("fastq.gz"): shutil.copy(os.path.join(fastqPath, fastqFile), os.path.join(workingPath, "data/", species, "fastq/totalRNA/"))
        # TO DO: handle multiple species?
        
    print('Directories built!')