import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# variety of quality control (QC) functions that can be performed once tagdirs are made

# primary function to be called
def qc(workingPath, species_lst, mode, errFile):
    print('Beginning QC')
    totalRnaTag = '_RNA'
    output = f'{workingPath}files/QC'
    if mode == 'star': totalRnaTag = '_RNA'
    elif mode == 'hisat': totalRnaTag = '_totalRNA'
    for species in species_lst:
        csRNAmedianTags(workingPath, species, errFile)                 # csRNA median tags
        csRNAtagsVsFrac(workingPath, species, errFile)                 # csRNA Tags vs Frac of pos.
        sRNAmedianTags(workingPath, species, errFile)                  # sRNA Median Tags
        sRNAtagsVsFrac(workingPath, species, errFile)                  # sRNA Tags vs Frac of pos.
        totalRNAmedianTags(workingPath, species, totalRnaTag, errFile) # totalRNA Median Tags
        totalRNAtagsVsFrac(workingPath, species, totalRnaTag, errFile) # totalRNA Tags vs Frac of pos.
        csRNAlengthPlot(workingPath, species, totalRnaTag, errFile)    # csRNA Length Plots
        csRNAcomboLengthPlot(workingPath, species, errFile)            # csRNA Combo Length Plots
        sRNAlengthPlot(workingPath, species, errFile)                  # sRNA Length Plots
        sRNAcomboLengthPlot(workingPath, species, errFile)             # sRNA Combo Length Plots
        ntPrefs(workingPath, species, errFile, output)                 # nt Preferences
        csRNAaPlot(workingPath, species, errFile)                      # csRNA A-plots
        csRNAcomboAPlot(workingPath, species, errFile)                 # csRNA Combo A-plots
        sRNAaPlot(workingPath, species, errFile)                       # sRNA A-plots
        sRNAcomboAPlot(workingPath, species, errFile)                  # sRNA Combo A-plots

def csRNAmedianTags(workingPath, species, errFile):
    parentDir = f'{workingPath}data/{species}/tagDirs'
    os.chdir(parentDir)
    pwd = os.getcwd()
    # loop through tag dirs
    for eachTagDir in os.listdir():
        if 'sRNA' in eachTagDir:
            # make a ucsc file of the tagdir
            cmd = f'makeUCSCfile {eachTagDir} -strand separate -style tss > ../bedgraphs/{eachTagDir}.ucsc.bedGraph'
            print(cmd)
            os.system(cmd)
    output = f'{workingPath}files/QC'
    os.chdir(f'{parentDir}')
    pwd = os.getcwd()
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sorted(sampFile):
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            if 'csRNA' not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of EmptyDataError
    try:
        patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
        patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Val'})
        patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]        
        my_dict = {"Library":[],"Median tags per tag position (should be =1)":[]}        
        for f in tagdirs:
            file = f +'/tagCountDistribution.txt'
            # read in file
            read_tagCountDist_file = pd.read_csv(file, sep='\t')        
            name = list(read_tagCountDist_file.columns)
            median_val = str(name).split('=')[1].split(',')[0]        
            plotName = f.split('/')[-1]
            my_dict["Library"].append(plotName)
            my_dict["Median tags per tag position (should be =1)"].append(median_val)
        Median_frame = pd.DataFrame(my_dict)
        Median_frame["Median tags per tag position (should be =1)"] = pd.to_numeric(Median_frame["Median tags per tag position (should be =1)"]) # convert column to numeric values for plotting
        graph = sns.barplot(data=Median_frame,  y = 'Library', x="Median tags per tag position (should be =1)", capsize=.4, errcolor=".5", linewidth=2, edgecolor=".5", facecolor=(0, 0, 0, 0),)
        graph.axvline(1, color='g', ls='--')
        graph.axvline(1.2, color='r')        
        sns.set(rc={'figure.figsize':(6,6)})
        plt.tight_layout()
        plt.savefig(f'{output}/csRNA_{species}_medianTagsPerPosition.png')
        plt.close()
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with csRNA present for medianTagsPerPosition\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(OSError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'File not Found {tagdirs[1]}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(pd.errors.EmptyDataError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'Empty File Found {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def csRNAtagsVsFrac(workingPath, species, errFile):
    parentDir = f'{workingPath}data/{species}/tagDirs'
    path = f'{parentDir}'
    output_dir = f'{workingPath}files/QC'
    os.chdir(f'{parentDir}/../fastq')
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sorted(sampFile):
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            if 'csRNA' not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of EmptyDataError or FileNotFoundError
    try:
        patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
        patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Tags per tag position'})
        patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]        
        for f in tagdirs:
            file = f +'/tagCountDistribution.txt'
            # read in file
            read_tagCountDist_file = pd.read_csv(file, sep='\t')
            # define variables for names
            median_value = str(list(read_tagCountDist_file)).split('=')[1].split(',')[0].split(' ')[1] # gets the median value out of the column name
            modName = f.split('/')[-1]
            column_name = modName + ' (' + median_value + ')'        
            # rename columns for merge
            read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[0]: 'Tags per tag position'})
            read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[1]: column_name})
            # cat all tagdir values together
            merged_frame = pd.merge(patens_tagCountDist_start, read_tagCountDist_file, left_on='Tags per tag position', right_on='Tags per tag position', how='left')
            patens_tagCountDist_start = merged_frame.copy()        
        patens_tagCountDist_start = patens_tagCountDist_start.set_index('Tags per tag position')
        patens_tagCountDist_start.to_csv(output_dir + f'/{species}_csRNA_sum_TagCountDist.txt', sep = '\t')
        # convert into a df for scatter plots
        combined_frames = patens_tagCountDist_start.replace(['0', 0], np.nan) # make 0 to NaN
        combined_frames_stacked = combined_frames.stack().reset_index() # stack and reset index
        combined_frames_stacked.columns = ['Tags per tag position','Library','Fraction of Positions'] # rename columns
        combined_frames_stacked_log = combined_frames_stacked
        combined_frames_stacked_log['Tags per tag position'] = np.log(combined_frames_stacked_log['Tags per tag position']) # log
        combined_frames_stacked_log['Fraction of Positions'] = np.log(combined_frames_stacked_log['Fraction of Positions']) # log
        g = sns.FacetGrid(combined_frames_stacked_log, col='Library',height=8, aspect=1, col_wrap=4)
        plots = g.map(sns.scatterplot, "Tags per tag position", "Fraction of Positions", alpha=.5)
        plt.savefig(f'{output}/csRNA_{species}_tagsPer_Vs_FracofPos.png') # save plot as a png
        plt.savefig(f'{output}/csRNA_{species}_tagsPer_Vs_FracofPos.svg') # and as an svg
        plt.close()
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with csRNA present for tagsPer_Vs_FracofPos\nErr. at line no. {exception_traceback.tb_lineno}\n')
    except(pd.errors.EmptyDataError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'Empty File Found {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def sRNAmedianTags(workingPath, species, errFile):
    parentDir = f'{workingPath}data/{species}/tagDirs'
    output = f'{workingPath}files/QC'
    os.chdir(f'{parentDir}')
    pwd = os.getcwd()
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sorted(sampFile):
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            if '_sRNA' not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
        patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Val'})
        patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]        
        my_dict = {"Library":[],"Median tags per tag position (should be =1)":[]}        
        for f in tagdirs:
            file = f +'/tagCountDistribution.txt'
            # read in file
            read_tagCountDist_file = pd.read_csv(file, sep='\t')        
            name = list(read_tagCountDist_file.columns)
            median_val = str(name).split('=')[1].split(',')[0]        
            plotName = f.split('/')[-1]
            my_dict["Library"].append(plotName)
            my_dict["Median tags per tag position (should be =1)"].append(median_val)        
        Median_frame = pd.DataFrame(my_dict)
        Median_frame["Median tags per tag position (should be =1)"] = pd.to_numeric(Median_frame["Median tags per tag position (should be =1)"]) # convert column to numeric values for plotting
        graph = sns.barplot(data=Median_frame,  y = 'Library', x="Median tags per tag position (should be =1)", capsize=.4, errcolor=".5", linewidth=2, edgecolor=".5", facecolor=(0, 0, 0, 0),)
        graph.axvline(1, color='g', ls='--')
        graph.axvline(1.2, color='r')        
        sns.set(rc={'figure.figsize':(6,6)})
        plt.tight_layout()
        plt.savefig(f'{output}/sRNA_{species}_medianTagsPerPosition.png') # save plot as a png
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with _sRNA present for medianTagsPerPosition\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def sRNAtagsVsFrac(workingPath, species, errFile):
    parentDir = f'{workingPath}data/{species}/tagDirs'
    path = f'{parentDir}'
    output_dir = f'{workingPath}files/QC'
    os.chdir(f'{parentDir}/../fastq')    
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sorted(sampFile):
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            if '_sRNA' not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
        patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Tags per tag position'})
        patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]        
        for f in tagdirs:
            file = f +'/tagCountDistribution.txt'
            # read in file
            read_tagCountDist_file = pd.read_csv(file, sep='\t')
            # define variables for names
            median_value = str(list(read_tagCountDist_file)).split('=')[1].split(',')[0].split(' ')[1] # gets the median value out of the column name
            modName = f.split('/')[-1]
            column_name = modName + ' (' + median_value + ')'        
            # rename columns for merge
            read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[0]: 'Tags per tag position'})
            read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[1]: column_name})
            # cat all tagDir values together
            merged_frame = pd.merge(patens_tagCountDist_start, read_tagCountDist_file, left_on='Tags per tag position', right_on='Tags per tag position', how='left')
            patens_tagCountDist_start = merged_frame.copy()        
        patens_tagCountDist_start = patens_tagCountDist_start.set_index('Tags per tag position')
        patens_tagCountDist_start.to_csv(output_dir + f'/{species}_sRNA_sum_TagCountDist.txt', sep = '\t')        
        # convert into a df for scatter plots
        combined_frames = patens_tagCountDist_start.replace(['0', 0], np.nan)  # make 0 to NaN 
        combined_frames_stacked = combined_frames.stack().reset_index() # stack and reset index
        combined_frames_stacked.columns = ['Tags per tag position','Library','Fraction of Positions'] # rename columns
        combined_frames_stacked_log = combined_frames_stacked
        combined_frames_stacked_log['Tags per tag position'] = np.log(combined_frames_stacked_log['Tags per tag position']) # log
        combined_frames_stacked_log['Fraction of Positions'] = np.log(combined_frames_stacked_log['Fraction of Positions']) # log
        g = sns.FacetGrid(combined_frames_stacked_log, col='Library',height=8, aspect=1, col_wrap=4)
        plots = g.map(sns.scatterplot, "Tags per tag position", "Fraction of Positions", alpha=.5)
        plt.savefig(f'{output}/sRNA_{species}_tagsPer_Vs_FracofPos.png') # save plot as a png
        plt.savefig(f'{output}/sRNA_{species}_tagsPer_Vs_FracofPos.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with _sRNA present for tagsPer_Vs_FracofPos\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def totalRNAmedianTags(workingPath, species, totalRnaTag, errFile):
    parentDir = f'{workingPath}data/{species}/tagDirs'
    output = f'{workingPath}files/QC'
    os.chdir(f'{parentDir}')
    pwd = os.getcwd()
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sorted(sampFile):
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            if totalRnaTag not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:    
        patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
        patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Val'})
        patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]        
        my_dict = {"Library":[],"Median tags per tag position (should be =1)":[]}
        for f in tagdirs:
            file = f +'/tagCountDistribution.txt'
            # read in file
            read_tagCountDist_file = pd.read_csv(file, sep='\t')        
            name = list(read_tagCountDist_file.columns)
            median_val = str(name).split('=')[1].split(',')[0]        
            plotName = f.split('/')[-1]
            my_dict["Library"].append(plotName)
            my_dict["Median tags per tag position (should be =1)"].append(median_val)        
        Median_frame = pd.DataFrame(my_dict)
        Median_frame["Median tags per tag position (should be =1)"] = pd.to_numeric(Median_frame["Median tags per tag position (should be =1)"]) # convert column to numeric values for plotting        
        graph = sns.barplot(data=Median_frame,  y = 'Library', x="Median tags per tag position (should be =1)", capsize=.4, errcolor=".5", linewidth=2, edgecolor=".5", facecolor=(0, 0, 0, 0),)
        graph.axvline(1, color='g', ls='--')
        graph.axvline(1.2, color='r')        
        sns.set(rc={'figure.figsize':(6,6)})
        plt.savefig(f'{output}/totalRNA_medianTagsPerPosition.png') # save plot as a png
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with '+totalRnaTag+' present for medianTagsPerPosition\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def totalRNAtagsVsFrac(workingPath, species, totalRnaTag, errFile):
    parentDir = f'{workingPath}data/{species}/tagDirs'
    path = f'{parentDir}'
    output_dir = f'{workingPath}files/QC'
    os.chdir(f'{parentDir}/../fastq')
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sorted(sampFile):
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            if totalRnaTag not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        patens_tagCountDist_start = pd.read_csv(tagdirs[1] + '/tagCountDistribution.txt', sep='\t')
        patens_tagCountDist_start = patens_tagCountDist_start.rename(columns={patens_tagCountDist_start.columns[0]: 'Tags per tag position'})
        patens_tagCountDist_start = patens_tagCountDist_start.iloc[:,[0]]        
        for f in tagdirs:
            file = f +'/tagCountDistribution.txt'
            # read in file
            read_tagCountDist_file = pd.read_csv(file, sep='\t')
            # define variables for names
            median_value = str(list(read_tagCountDist_file)).split('=')[1].split(',')[0].split(' ')[1] # gets the median value out of the column name
            modName = f.split('/')[-1]
            column_name = modName + ' (' + median_value + ')'        
            # rename columns for merge
            read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[0]: 'Tags per tag position'})
            read_tagCountDist_file = read_tagCountDist_file.rename(columns={read_tagCountDist_file.columns[1]: column_name})
            # cat all tagDir values together
            merged_frame = pd.merge(patens_tagCountDist_start, read_tagCountDist_file, left_on='Tags per tag position', right_on='Tags per tag position', how='left')
            patens_tagCountDist_start = merged_frame.copy()        
        patens_tagCountDist_start = patens_tagCountDist_start.set_index('Tags per tag position')
        patens_tagCountDist_start.to_csv(output_dir + f'/{species}_totalRNA_sum_TagCountDist.txt', sep = '\t')        
        # convert into a df for scatter plots
        combined_frames = patens_tagCountDist_start.replace(['0', 0], np.nan) # make 0 to NaN
        combined_frames_stacked = combined_frames.stack().reset_index() # stack and reset index
        combined_frames_stacked.columns = ['Tags per tag position','Library','Fraction of Positions'] # rename columns
        combined_frames_stacked_log = combined_frames_stacked
        combined_frames_stacked_log['Tags per tag position'] = np.log(combined_frames_stacked_log['Tags per tag position']) # log
        combined_frames_stacked_log['Fraction of Positions'] = np.log(combined_frames_stacked_log['Fraction of Positions']) # log
        g = sns.FacetGrid(combined_frames_stacked_log, col='Library',height=8, aspect=1, col_wrap=4)
        plots = g.map(sns.scatterplot, "Tags per tag position", "Fraction of Positions", alpha=.5)
        plt.savefig(f'{output}/totalRNA_{species}_tagsPer_Vs_FracofPos.png') # save plot as a png
        plt.savefig(f'{output}/totalRNA_{species}_tagsPer_Vs_FracofPos.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with'+totalRnaTag+'present for frac of pos.\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def csRNAlengthPlot(workingPath, species, totalRnaTag, errFile):
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            if 'ChIPseq' in eachLine:
                continue
            if totalRnaTag in eachLine:
                continue
            if '_sRNA' in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        tagLength_start = pd.read_csv(tagdirs[1] + '/tagLengthDistribution.txt', sep='\t')
        tagLength_start = tagLength_start.rename(columns={tagLength_start.columns[0]: 'Length (nt)'})
        tagLength_start = tagLength_start.iloc[:,[0]]
        for f in sorted(tagdirs):
            file = f +'/tagLengthDistribution.txt'
            # read in file
            read_tagLength_file = pd.read_csv(file, sep='\t')
            # define variables for names
            column_name = f.split('/')[-1]        
            # rename columns for merge
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[0]: 'Length (nt)'})
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[1]: column_name})
            # cat all tagDir values together
            merged_frame = pd.merge(tagLength_start, read_tagLength_file, left_on='Length (nt)', right_on='Length (nt)', how='left')
            tagLength_start = merged_frame.copy()        
        tagLength_start = tagLength_start.iloc[1: , :] # delete 0 because not a read after trimming
        tagLength_start = tagLength_start.set_index('Length (nt)')
        tagLength_start.to_csv(output_dir + f'/{species}_csRNA_sum_TagLengthDist.txt', sep = '\t')        
        combined_frames = tagLength_start.replace(['0', 0], np.nan) # make 0 to NaN 
        combined_frames_stacked = combined_frames.stack().reset_index() # stack and reset index
        combined_frames_stacked.columns = ['Length (nt)','Library','Fraction of Reads'] # rename columns
        ax = sns.lineplot(data=combined_frames_stacked, x="Length (nt)", y='Fraction of Reads', hue="Library", linewidth=2, alpha = 0.6)
        sns.set(rc={'figure.figsize':(20,15)}, font_scale=1)
        sns.move_legend(ax, "upper right")
        plt.savefig(f'{output}/csRNA_{species}_Length_plot.png') # save plot as a png
        plt.savefig(f'{output}/csRNA_{species}_Length_plot.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with not ChIPseq _sRNA or '+totalRnaTag+' present for length plots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def csRNAcomboLengthPlot(workingPath, species, errFile):
    tagDirectory = f'{workingPath}data/{species}/tagDirs/'
    tagdirs = []
    for eachFile in os.listdir(tagDirectory):
        # TO DO (?): code this for the csRNA
        if 'csRNA' not in eachFile:
            continue
        elif '-r' in eachFile:
            continue
        else:
            tagdirs.append(tagDirectory+'/'+ eachFile)        
    file = '/tagLengthDistribution.txt' # defining early in case of reference in FileNotFoundError
    try:
        tagLength_start = pd.read_csv(tagdirs[1] + '/tagLengthDistribution.txt', sep='\t')
        tagLength_start = tagLength_start.rename(columns={tagLength_start.columns[0]: 'Length (nt)'})
        tagLength_start = tagLength_start.iloc[:,[0]]
        for f in sorted(tagdirs):
            file = f +'/tagLengthDistribution.txt'
            # read in file
            read_tagLength_file = pd.read_csv(file, sep='\t')
            # define variables for names
            column_name = f.split('/')[-1]
            # rename columns for merge
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[0]: 'Length (nt)'})
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[1]: column_name})
            # cat all tagDir values together
            merged_frame = pd.merge(tagLength_start, read_tagLength_file, left_on='Length (nt)', right_on='Length (nt)', how='left')
            tagLength_start = merged_frame.copy()
        tagLength_start = tagLength_start.iloc[1: , :] # delete 0 cause not a read after trimming
        tagLength_start = tagLength_start.set_index('Length (nt)')
        tagLength_start.to_csv(output_dir + f'/{species}_csRNA_sum_TagLengthDist.txt', sep = '\t')
        combined_frames = tagLength_start.replace(['0', 0], np.nan) # make 0 to NaN 
        combined_frames_stacked = combined_frames.stack().reset_index() # stack and reset index
        combined_frames_stacked.columns = ['Length (nt)','Library','Fraction of Reads'] # rename columns
        ax = sns.lineplot(data=combined_frames_stacked, x="Length (nt)", y='Fraction of Reads', hue="Library", linewidth=2, alpha = 0.6)
        sns.set(rc={'figure.figsize':(20,15)}, font_scale=1)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f'{output}/csRNA_{species}_Combo_Length_plot.png') # save plot as a png
        plt.savefig(f'{output}/csRNA_{species}_Combo_Length_plot.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with not csRNA present for length plots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def sRNAlengthPlot(workingPath, species, errFile):
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            if 'ChIPseq' in eachLine or '_RNA' in eachLine or '_totalRNA' in eachLine or 'csRNA' in eachLine:
                continue
            elif '_sRNA' in eachLine:
                tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        tagLength_start = pd.read_csv(tagdirs[1] + '/tagLengthDistribution.txt', sep='\t')
        tagLength_start = tagLength_start.rename(columns={tagLength_start.columns[0]: 'Length (nt)'})
        tagLength_start = tagLength_start.iloc[:,[0]]        
        for f in sorted(tagdirs):
            file = f +'/tagLengthDistribution.txt'
            # read in file
            read_tagLength_file = pd.read_csv(file, sep='\t')
            # define variables for names
            column_name = f.split('/')[-1]        
            # rename columns for merge
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[0]: 'Length (nt)'})
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[1]: column_name})
            # cat all tagDir values together
            merged_frame = pd.merge(tagLength_start, read_tagLength_file, left_on='Length (nt)', right_on='Length (nt)', how='left')
            tagLength_start = merged_frame.copy()        
        tagLength_start = tagLength_start.iloc[1: , :] # delete 0 cause not a read after trimming
        tagLength_start = tagLength_start.set_index('Length (nt)')
        tagLength_start.to_csv(output_dir + f'/{species}_sRNA_sum_TagLengthDist.txt', sep = '\t')        
        combined_frames = tagLength_start.replace(['0', 0], np.nan) # make 0 to NaN 
        combined_frames_stacked = combined_frames.stack().reset_index() # stack and reset index
        combined_frames_stacked.columns = ['Length (nt)','Library','Fraction of Reads'] # rename columns
        ax = sns.lineplot(data=combined_frames_stacked, x="Length (nt)", y='Fraction of Reads', hue="Library", linewidth=2, alpha = 0.6)
        sns.set(rc={'figure.figsize':(20,15)}, font_scale=1)
        sns.move_legend(ax, "upper right")
        plt.savefig(f'{output}/sRNA_{species}_Length_plot.png') # save plot as a png
        plt.savefig(f'{output}/sRNA_{species}_Length_plot.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with _sRNA present for length plots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def sRNAcomboLengthPlot(workingPath, species, errFile):
    tagDirectory = f'{workingPath}data/{species}/tagDirs/'
    tagdirs = []  
    for eachFile in os.listdir(tagDirectory):
        # code this for the sRNA
        if '_sRNA' not in eachFile:
            continue
        elif '-r' in eachFile:
            continue
        else:
            tagdirs.append(tagDirectory+'/'+ eachFile)        
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        tagLength_start = pd.read_csv(tagdirs[1] + '/tagLengthDistribution.txt', sep='\t')
        tagLength_start = tagLength_start.rename(columns={tagLength_start.columns[0]: 'Length (nt)'})
        tagLength_start = tagLength_start.iloc[:,[0]]
        for f in sorted(tagdirs):
            file = f +'/tagLengthDistribution.txt'
            # read in file
            read_tagLength_file = pd.read_csv(file, sep='\t')
            # define variables for names
            column_name = f.split('/')[-1]
            # rename columns for merge
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[0]: 'Length (nt)'})
            read_tagLength_file = read_tagLength_file.rename(columns={read_tagLength_file.columns[1]: column_name})
            # cat all tagDir values together
            merged_frame = pd.merge(tagLength_start, read_tagLength_file, left_on='Length (nt)', right_on='Length (nt)', how='left')
            tagLength_start = merged_frame.copy()
        tagLength_start = tagLength_start.iloc[1: , :] # delete 0 cause not a read after trimming
        tagLength_start = tagLength_start.set_index('Length (nt)')
        tagLength_start.to_csv(output_dir + f'/{species}_csRNA_sum_TagLengthDist.txt', sep = '\t')
        combined_frames = tagLength_start.replace(['0', 0], np.nan) # make 0 to NaN 
        combined_frames_stacked = combined_frames.stack().reset_index() # stack and reset index
        combined_frames_stacked.columns = ['Length (nt)','Library','Fraction of Reads'] # rename columns
        ax = sns.lineplot(data=combined_frames_stacked, x="Length (nt)", y='Fraction of Reads', hue="Library", linewidth=2, alpha = 0.6)
        sns.set(rc={'figure.figsize':(20,15)}, font_scale=1)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f'{output}/sRNA_{species}_Combo_Length_plot.png') # save plot as a png
        plt.savefig(f'{output}/sRNA_{species}_Combo_Length_plot.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with not sRNA present for length plots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def ntPrefs(workingPath, species, errFile, output):
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        # complicated way of making a column of 0-999, but just in case someone wants wider windows in the nt freq output tool upstream
        nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep='\t')
        nt_freq_file_start = nt_freq_file_start.iloc[:,[0]]        
        for f in tagdirs:
            nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
            nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
            plot_me = nt_freq_file.set_index('Offset')
            plot_me_stack = plot_me.stack().reset_index() # stack and reset index
            axisName = f.split('/')[-1]
            plot_me_stack.columns = ['Distance from TSS', 'nt', axisName + ' - %'] # rename columns
            sns.set(rc={'figure.figsize':(20,15)}, font_scale=2)
            plt.figure()
            g = sns.lineplot(data = plot_me_stack, x = "Distance from TSS", y = axisName + ' - %', hue = "nt")
            plt.tight_layout()
            plt.savefig(f'{output}/{axisName}_nt_Preference.png') # save plot as a png
            plt.savefig(f'{output}/{axisName}_nt_Preference.svg') # and as an svg
            plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'sampeInfo.txt file appears to be empty\nCheck your fastq directory\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr at lin no. {exception_traceback.tb_lineno}\n')
        
def csRNAaPlot(workingPath, species, errFile):
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt', 'r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/' + eachLine
            if 'csRNA' not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        # select key parameter
        nucleotide_to_plot = 'A'
        selected_files = 'csRNA'        
        # complicated way of making a column of 0-999, but just in case if someone wants wider windows in the nt freq output tool upstream
        nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep = '\t')
        nt_freq_file_start = nt_freq_file_start[nt_freq_file_start.columns[:5]]
        nt_freq_file_start = nt_freq_file_start.set_index('Offset')
        nt_freq_file_start = nt_freq_file_start.stack().reset_index() # stack and reset index
        nt_freq_file_start = nt_freq_file_start.rename(columns={nt_freq_file_start.columns[1]: 'nt'})
        nt_freq_file_start = nt_freq_file_start.iloc[: , :2] # select both offset and nt info
        unused = []
        for f in tagdirs:
            if selected_files in f:
                nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
                nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
                plot_me = nt_freq_file.set_index('Offset')
                plot_me_stack = plot_me.stack().reset_index() # stack and reset index
                modName = f.split('/')[-1]
                plot_me_stack.columns = ['Distance from TSS','nt', modName] # rename columns        
                del plot_me_stack['Distance from TSS']
                del plot_me_stack['nt']        
                merged_frame = pd.merge(nt_freq_file_start, plot_me_stack, left_index=True, right_index=True, how='left')
                nt_freq_file_start = merged_frame.copy()
            else:
                # append unused values to a list
                unused.append(f)        
        AntPlot_frame = nt_freq_file_start[nt_freq_file_start["nt"].str.contains(nucleotide_to_plot)==True]
        del AntPlot_frame['nt']
        AntPlot_frame = AntPlot_frame.set_index('Offset')
        AntPlot_frame_stack = AntPlot_frame.stack().reset_index() # stack and reset index
        AntPlot_frame_stack.columns = ['Distance to TSS','Library', nucleotide_to_plot + ' [%]'] # rename columns
        # plot with legend next to plot
        sns.set(rc={'figure.figsize':(20,15)},font_scale=1)
        plt.figure()
        ax = sns.lineplot(data=AntPlot_frame_stack, x="Distance to TSS", y= nucleotide_to_plot + ' [%]', hue="Library", linewidth=2, alpha = 0.6)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f'{output}/csRNA_{species}_Aplot.png') # save plot as a png
        plt.savefig(f'{output}/csRNA_{species}_Aplot.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with csRNA present for Aplots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def csRNAcomboAPlot(workingPath, species, errFile):
    tagDirectory = f'{workingPath}data/{species}/tagDirs/'
    tagdirs = []
    for eachFile in os.listdir(tagDirectory):
        # code this for the csRNA
        if 'csRNA' not in eachFile:
            continue
        elif '-r' in eachFile:
            continue
        else:
            tagdirs.append(tagDirectory+'/'+ eachFile)    
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        # select key parameter
        nucleotide_to_plot = 'A'
        selected_files = 'csRNA'    
        # complicated way of making a column of 0-999 but just in case if someone wants wider windows in the nt freq output tool upstream
        nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep='\t')
        nt_freq_file_start = nt_freq_file_start[nt_freq_file_start.columns[:5]]
        nt_freq_file_start = nt_freq_file_start.set_index('Offset')
        nt_freq_file_start = nt_freq_file_start.stack().reset_index() # stack and reset index
        nt_freq_file_start = nt_freq_file_start.rename(columns={nt_freq_file_start.columns[1]: 'nt'})
        nt_freq_file_start = nt_freq_file_start.iloc[: , :2] # select both offset and nt info
        unused = []
        for f in tagdirs:
            if selected_files in f:
                nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
                nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
                plot_me = nt_freq_file.set_index('Offset')
                plot_me_stack = plot_me.stack().reset_index() # stack and reset index
                modName = f.split('/')[-1]
                plot_me_stack.columns = ['Distance from TSS','nt', modName] # rename columns
                del plot_me_stack['Distance from TSS']
                del plot_me_stack['nt']    
                merged_frame = pd.merge(nt_freq_file_start, plot_me_stack, left_index=True, right_index=True, how='left')
                nt_freq_file_start = merged_frame.copy()
            else: # append unused values to a list
                unused.append(f)    
        AntPlot_frame = nt_freq_file_start[nt_freq_file_start["nt"].str.contains(nucleotide_to_plot)==True]
        del AntPlot_frame['nt']
        AntPlot_frame = AntPlot_frame.set_index('Offset')
        AntPlot_frame_stack = AntPlot_frame.stack().reset_index() # stack and reset index
        AntPlot_frame_stack.columns = ['Distance to TSS','Library', nucleotide_to_plot + ' [%]'] # rename columns
        # plot with legend next to plot
        sns.set(rc={'figure.figsize':(20,15)},font_scale=1)
        plt.figure()
        ax = sns.lineplot(data=AntPlot_frame_stack, x="Distance to TSS", y= nucleotide_to_plot + ' [%]', hue="Library", linewidth=2, alpha = 0.6)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f'{output}/csRNA_{species}_Combo_Aplots.png') # save plot as a png
        plt.savefig(f'{output}/csRNA_{species}_Combo_Aplots.svg') # and as an svg
        plt.close()
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with csRNA present for Aplots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def sRNAaPlot(workingPath, species, errFile):
    tagdirs = []
    with open(f'{workingPath}data/{species}/fastq/sampleInfo.txt','r') as sampFile:
        for eachLine in sampFile:
            eachLine = eachLine.strip('\n')
            eachLine = eachLine.split('\t')[0]
            eachLine = '/'+eachLine
            if '_sRNA' not in eachLine:
                continue
            tagdirs.append(eachLine)
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        # select key parameter
        nucleotide_to_plot = 'A'
        selected_files = '_sRNA'        
        # complicated way of making a column of 0-999 but just in case if someone wants wider windows in the nt freq output tool upstream
        nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep='\t')
        nt_freq_file_start = nt_freq_file_start[nt_freq_file_start.columns[:5]]
        nt_freq_file_start = nt_freq_file_start.set_index('Offset')
        nt_freq_file_start = nt_freq_file_start.stack().reset_index() # stack and reset index
        nt_freq_file_start = nt_freq_file_start.rename(columns={nt_freq_file_start.columns[1]: 'nt'})
        nt_freq_file_start = nt_freq_file_start.iloc[: , :2] # select both offset and nt info
        unused = []
        for f in tagdirs:
            if selected_files in f:
                nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
                nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
                plot_me = nt_freq_file.set_index('Offset')
                plot_me_stack = plot_me.stack().reset_index() # stack and reset index
                modName = f.split('/')[-1]
                plot_me_stack.columns = ['Distance from TSS','nt', modName] # rename columns
                del plot_me_stack['Distance from TSS']
                del plot_me_stack['nt']        
                merged_frame = pd.merge(nt_freq_file_start, plot_me_stack, left_index=True, right_index=True, how='left')
                nt_freq_file_start = merged_frame.copy()
            else:
                # append unused values to a list
                unused.append(f)        
        AntPlot_frame = nt_freq_file_start[nt_freq_file_start["nt"].str.contains(nucleotide_to_plot)==True]
        del AntPlot_frame['nt']
        AntPlot_frame = AntPlot_frame.set_index('Offset')
        AntPlot_frame_stack = AntPlot_frame.stack().reset_index() # stack and reset index
        AntPlot_frame_stack.columns = ['Distance to TSS','Library', nucleotide_to_plot + ' [%]'] # rename columns
        # plot with legend next to plot
        sns.set(rc={'figure.figsize':(20,15)},font_scale=1)
        plt.figure()
        ax = sns.lineplot(data=AntPlot_frame_stack, x="Distance to TSS", y= nucleotide_to_plot + ' [%]', hue="Library", linewidth=2, alpha = 0.6)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f'{output}/sRNA_{species}_Aplot.png') # save plot as a png
        plt.savefig(f'{output}/sRNA_{species}_Aplot.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with _sRNA present for Aplots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
        
def sRNAcomboAPlot(workingPath, species, errFile):
    tagDirectory = f'{workingPath}data/{species}/tagDirs/'
    tagdirs = []
    for eachFile in os.listdir(tagDirectory):
        # code this for the cRNA
        if '_sRNA' not in eachFile:
            continue
        elif '-r' in eachFile:
            continue
        else:
            tagdirs.append(tagDirectory+'/'+ eachFile)    
    file = '/tagCountDistribution.txt' # declaration in case of FileNotFoundError
    try:
        # select key parameter
        nucleotide_to_plot = 'A'
        selected_files = '_sRNA'    
        # complicated way of making a column of 0-999 but just in case if someone wants wider windows in the nt freq output tool upstream
        nt_freq_file_start = pd.read_csv(tagdirs[1] + '/tagFreqUniq.txt', sep='\t')
        nt_freq_file_start = nt_freq_file_start[nt_freq_file_start.columns[:5]]
        nt_freq_file_start = nt_freq_file_start.set_index('Offset')
        nt_freq_file_start = nt_freq_file_start.stack().reset_index() # stack and reset index
        nt_freq_file_start = nt_freq_file_start.rename(columns={nt_freq_file_start.columns[1]: 'nt'})
        nt_freq_file_start = nt_freq_file_start.iloc[: , :2] # select both offset and nt info    
        unused = []
        for f in tagdirs:
            if selected_files in f:
                nt_freq_file = pd.read_csv(f + '/tagFreqUniq.txt', sep='\t')
                nt_freq_file = nt_freq_file[nt_freq_file.columns[:5]]
                plot_me = nt_freq_file.set_index('Offset')
                plot_me_stack = plot_me.stack().reset_index() # stack and reset index
                modName = f.split('/')[-1]
                plot_me_stack.columns = ['Distance from TSS','nt', modName] # rename columns
                del plot_me_stack['Distance from TSS']
                del plot_me_stack['nt']    
                merged_frame = pd.merge(nt_freq_file_start, plot_me_stack, left_index=True, right_index=True, how='left')
                nt_freq_file_start = merged_frame.copy()
            else:
                # append unused values to a list
                unused.append(f)    
        AntPlot_frame = nt_freq_file_start[nt_freq_file_start["nt"].str.contains(nucleotide_to_plot)==True]
        del AntPlot_frame['nt']
        AntPlot_frame = AntPlot_frame.set_index('Offset')
        AntPlot_frame_stack = AntPlot_frame.stack().reset_index() # stack and reset index
        AntPlot_frame_stack.columns = ['Distance to TSS','Library', nucleotide_to_plot + ' [%]'] # rename columns
        # plot with legend next to plot
        sns.set(rc={'figure.figsize':(20,15)},font_scale=1)
        plt.figure()
        ax = sns.lineplot(data=AntPlot_frame_stack, x="Distance to TSS", y= nucleotide_to_plot + ' [%]', hue="Library", linewidth=2, alpha = 0.6)
        sns.move_legend(ax, "upper right")
        plt.tight_layout()
        plt.savefig(f'{output}/sRNA_{species}_Combo_Aplots.png') # save plot as a png
        plt.savefig(f'{output}/sRNA_{species}_Combo_Aplots.svg') # and as an svg
        plt.close() 
    except(IndexError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No lines in samplenfo.txt with csRNA present for Aplots\nErr. at lin no. {exception_traceback.tb_lineno}\n')
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        errFile.write(f'No file found for {file}\nErr. at lin no. {exception_traceback.tb_lineno}\n')