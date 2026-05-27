import os # operation system functions
import sys # system functions
import subprocess # subproccess systm package. similar to os functions

# Functions for finding peaks, analyzing TSRs, etc.

def peakfinding(workingPath, species_lst):
    dataPath = f'{workingPath}data/'
    os.chdir(dataPath)
    for eachSpecies in species_lst:
        speciesPath = f'{dataPath}{eachSpecies}/'
        tagdir_path = f'{dataPath}{eachSpecies}/tagDirs/'
        os.chdir(tagdir_path)
        cmdLst = []
        fileSheet = []
        for tagDir in os.listdir():
            # define the regex used to find all csRNA-r,r2,-r etc.
            regMatch = re.search('csRNA*-r[0-9]',tagDir)
            if tagDir.endswith('_csRNA'):
                sampleName = tagDir.split('_csRNA')[0]
                totalRNA = f'{tagdir_path}{sampleName}_RNA'
                sRNA = f'{tagdir_path}{sampleName}_sRNA'
                output = f'{workingPath}files/TSR/{sampleName}-{ntag}'
                genomefa = f'{workingPath}genomes/{eachSpecies}/*.fa'
                gtf = f'{workingPath}genomes/{eachSpecies}/*.gtf'
                # now run findcsRNATSS.pl for each sample
                if os.path.exists(totalRNA):
                    cmd = f' findcsRNATSS.pl {tagdir_path}{tagDir} -o {output} -i {sRNA} -rna {totalRNA} -genome {genomefa} -gtf {gtf} -ntagThreshold {ntag}'
                    cmdLst.append(cmd)
                elif not os.path.exists(totalRNA):
                    cmd = f' findcsRNATSS.pl {tagdir_path}{tagDir} -o {output} -i {sRNA} -genome {genomefa} -gtf {gtf} -ntagThreshold {ntag}'
                    cmdLst.append(cmd)
                # create a string of values that contains the files required for each sample to be processed
                # i.e. csrna file \t srna file \t rna file
                fileLst = f'{tagdir_path}{tagDir}\t{sRNA}\t{totalRNA}'
                fileSheet.append(fileLst)                
            elif (regMatch != None):
                sampleRep = regMatch.string.split('_csRNA')[-1]
                sampleName = regMatch.string.split('_csRNA')[0]
                totalRNA = f'{tagdir_path}{sampleName}_RNA{sampleRep}'
                sRNA = f'{tagdir_path}{sampleName}_sRNA{sampleRep}'
                output = f'{workingPath}files/TSR/{sampleName}{sampleRep}-{ntag}'
                genomefa = f'{workingPath}genomes/{eachSpecies}/*.fa'
                gtf = f'{workingPath}genomes/{eachSpecies}/*.gtf'
                # now run find csRNATSS.pl for each sample
                if os.path.exists(totalRNA):
                    cmd = f' findcsRNATSS.pl {tagdir_path}{tagDir} -o {output} -i {sRNA} -rna {totalRNA} -genome {genomefa} -gtf {gtf} -ntagThreshold {ntag}'
                    cmdLst.append(cmd)
                elif not os.path.exists(totalRNA):
                    cmd = f' findcsRNATSS.pl {tagdir_path}{tagDir} -o {output} -i {sRNA} -genome {genomefa} -gtf {gtf} -ntagThreshold {ntag}'
                    cmdLst.append(cmd)
                # create a string of values that contains the files required for each sample to be processed
                # i.e. csrna file \t srna file \t rna file
                fileLst = f'{tagdir_path}{tagDir}\t{sRNA}\t{totalRNA}'
                fileSheet.append(fileLst)
            else:
                continue
        if len(cmdLst) > 0:
            print('\nRunning PeakFinding commands:')
            procPool = multiprocessing.Pool(2)
            procPool.map(genericFunc, cmdLst)
            # close the pool
            procPool.close()
            procPool.join()        
        # now create a pandas dataframe and save it to a file location
        peakfindingFileDF = pd.DataFrame(fileSheet)
        peakfindingFileDF.to_csv(f'{workingPath}files/{eachSpecies}_Peakfinding_dirs.txt',sep='\t', index=None)
        
def TSRstabilityCount(workingPath, errFile):
    TSRPath = f'{workingPath}files/TSR/'
    os.chdir(TSRPath)
    stability_dict = {"species":[],"stableTSRs":[],"unstableTSRs":[]}
    for TSRFile in os.listdir():
        if TSRFile.endswith(f'-{ntag}.tss.txt'):
            try:
                TSR_frame = pd.read_csv(TSRFile, sep='\t')
                file_UDI = TSRFile.split('-10')[0]
                species = TSRFile.split('_')[0] + '_' + TSRFile.split('_')[1].split('-')[0]
                genomefa = f'{workingPath}genomes/{species}/genome.fa'+ species + '/genome.fa'        
                TSR_frame = TSR_frame[['#tssClusterID', 'chr', 'start', 'end', 'strand','Stable/Unstable']]
                TSR_frame_stable = TSR_frame[~TSR_frame['Stable/Unstable'].str.contains('unstable')==True]
                TSR_frame_stable.to_csv('tsrs_' + file_UDI + '_stable.txt', sep = '\t', index=False) #save for iTSS counts
                number_stableTSRs = len(TSR_frame_stable)        
                TSR_frame_unstable = TSR_frame[TSR_frame['Stable/Unstable'].str.contains('unstable')==True]
                TSR_frame_unstable.to_csv('tsrs_' + file_UDI + '_unstable.txt', sep = '\t', index=False)
                number_UNstableTSRs = len(TSR_frame_unstable)        
                stability_dict["species"].append(file_UDI)
                stability_dict["stableTSRs"].append(number_stableTSRs)
                stability_dict["unstableTSRs"].append(number_UNstableTSRs)                
            except(FileNotFoundError):
                exception_type, exception_object, exception_traceback = sys.exc_info()
                splitObjects = str(exception_object).split("'")
                errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')    
    stability_dictTSR_frame = pd.DataFrame(stability_dict)
    stability_dictTSR_frame = stability_dictTSR_frame.sort_values(by=['species'])
    stability_dictTSR_frame.to_csv(f'../Summary_{species}_TSR_stability_counts.tsv', sep = '\t')
    
def TSSfromReads(workingPath, errFile):
    os.chdir(iTSSPath)
    for eachSpecies in species_lst:
        speciesPath = f'{dataPath}{eachSpecies}/'
        speciesFastqPath = f'{dataPath}{eachSpecies}/fastq/'
        infoFile_path = speciesPath +'sampleInfoFile.txt'
        tagdir_path = f'{dataPath}{eachSpecies}/tagDirs/'
        for tagDir in os.listdir(tagdir_path):
            if tagDir.endswith('_csRNA'):
                try:
                    # name out all the files I will want to use for the program
                    csRNAfile = f'{tagdir_path}{tagDir}'
                    sRNAfile = tagDir.replace('csRNA','sRNA')
                    RNAfile = tagDir.replace('csRNA','RNA')
                    outputfile = 'iTSS_'+tagDir.split('_csRNA')[0]+'.txt'        
                    # run the command for getTSSfromReads.pl
                    cmd = f'getTSSfromReads.pl -d {csRNAfile} -dinput {tagdir_path}{sRNAfile} -min 7 > {outputfile}'
                    print(cmd)
                    proc = subprocess.Popen([cmd],shell=True).wait()                
                except(FileNotFoundError):
                    exception_type, exception_object, exception_traceback = sys.exc_info()
                    splitObjects = str(exception_object).split("'")
                    errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                    continue
                    
def plusTSSstats(workingPath, errFile):
    iTSSPath = f'{workingPath}files/iTSS/'
    os.chdir(iTSSPath)
    rmCmd = 'rm *coBoundBy*.txt'
    os.system(rmCmd) # remove the previous run attempts if any
    TSS_stability_dict = {"species":[],"stableTSSs":[],"unstableTSSs":[]}
    for eachSpecies in species_lst:
        speciesPath = f'{dataPath}{eachSpecies}/'
        speciesFastqPath = f'{dataPath}{eachSpecies}/fastq/'
        infoFile_path = speciesPath +'sampleInfoFile.txt'
        tagdir_path = f'{dataPath}{eachSpecies}/tagDirs/'
        for iTSS_file in os.listdir():
            if eachSpecies in iTSS_file:
                if iTSS_file.startswith('iTSS_'):
                    if iTSS_file.endswith('.txt'): # gets both stable & unstable
                        try:
                            UDI_file = iTSS_file.split('iTSS_')[1].split('.txt')[0]
                            prefix_stable = 'iTSS_' + UDI_file + '_stable'
                            prefix_UNstable = 'iTSS_' + UDI_file + '_unstable'
                            species_tss = iTSS_file.split('iTSS_')[1].split('.txt')[0]    
                            stable_tsr = '../TSR/tsrs_' + species_tss + '_stable.txt'
                            UNstable_tsr = '../TSR/tsrs_' + species_tss + '_unstable.txt'    
                            cmd1 = f'mergePeaks {iTSS_file} {stable_tsr} -strand -cobound 1 -prefix {prefix_stable}'    
                            os.system(cmd1)    
                            cmd2 = f'mergePeaks {iTSS_file} {UNstable_tsr} -strand -cobound 1 -prefix {prefix_UNstable}'
                            os.system(cmd2)    
                            stable_TSS = pd.read_csv(prefix_stable + '.coBoundBy1.txt', sep='\t')
                            number_stableTSRs = len(stable_TSS)-1    
                            UNstable_TSS = pd.read_csv(prefix_UNstable + '.coBoundBy1.txt', sep='\t')
                            number_UNstableTSRs = len(UNstable_TSS)-1    
                            TSS_stability_dict["species"].append(UDI_file)
                            TSS_stability_dict["stableTSSs"].append(number_stableTSRs)
                            TSS_stability_dict["unstableTSSs"].append(number_UNstableTSRs)                        
                        except(FileNotFoundError):
                            exception_type, exception_object, exception_traceback = sys.exc_info()
                            splitObjects = str(exception_object).split("'")
                            errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                            continue                        
            rmCmd = 'rm *stable.coBoundBy0.txt'
            os.system(rmCmd) # remove the cobound files i don't wish to keep    
            stability_dict_frame = pd.DataFrame(TSS_stability_dict)
            stability_dict_frame = stability_dict_frame.sort_values(by=['species'])
            stability_dict_frame.to_csv(f'../{species}_Summary_TSS_Stability_counts.tsv', sep = '\t')
            
def tagDirStats(workingPath, species_lst, errFile):
    for species in species_lst:
        tagDirsStats_dict = {"File":[],"UniquePositions":[],"TotalTags":[],"tagsPerBP":[],"tagPosinAnalysis":[],"averageTagsPerPosition":[],"averageFragmentGCcontent":[]}
        # jump down that species path
        speciesPath = f'{workingPath}data/{species}/'
        tagPath = f'{speciesPath}tagDirs/'
        os.chdir(tagPath)
        for tagDir in os.listdir():
            if os.path.exists(tagPath + tagDir + '/tagInfo.txt'):
                try:
                    tagInfo_frame = pd.read_csv(tagPath + tagDir + '/tagInfo.txt',  sep='\t')
                    UniquePositions = (tagInfo_frame.loc[0][1])
                    TotalTags = (tagInfo_frame.loc[0][2])
                    tagsPerBP = (tagInfo_frame.loc[3][0]).split('=')[1]
                    tagPosinAnalysis = str(float(tagInfo_frame.loc[0][2])/float(tagInfo_frame.loc[4][0].split('=')[-1]))
                    averageTagsPerPosition = (tagInfo_frame.loc[4][0]).split('=')[1]
                    averageFragmentGCcontent = (tagInfo_frame.loc[8][0]).split('=')[1]
                    tagDirsStats_dict["File"].append(tagDir)
                    tagDirsStats_dict["UniquePositions"].append(UniquePositions)
                    tagDirsStats_dict["TotalTags"].append(TotalTags)
                    tagDirsStats_dict["tagsPerBP"].append(tagsPerBP)
                    tagDirsStats_dict["tagPosinAnalysis"].append(tagPosinAnalysis)
                    tagDirsStats_dict["averageTagsPerPosition"].append(averageTagsPerPosition)
                    tagDirsStats_dict["averageFragmentGCcontent"].append(averageFragmentGCcontent)                    
                except(FileNotFoundError):
                    exception_type, exception_object, exception_traceback = sys.exc_info()
                    splitObjects = str(exception_object).split("'")
                    errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                    continue
                except(ZeroDivisionError):
                    exception_type, exception_object, exception_traceback = sys.exc_info()
                    errFile.write(f'Something wrong with tagDir {tagDir}\nErr. at lin no. {exception_traceback.tb_lineno}\n')
                    continue    
        tagDirsStats_dict_frame = pd.DataFrame(tagDirsStats_dict)
        tagDirsStats_dict_frame = tagDirsStats_dict_frame.sort_values(by=['File'])
        tagDirsStats_dict_frame.to_csv(f'{workingPath}files/{species}_summary_tagDirsStats.tsv', sep = '\t')
        
def TSRannotation(workingPath, errFile):
    TSRFile = f'{workingPath}files/TSR/'
    os.chdir(TSRFile)
    TSR_dict = {"Library":[],"TSRs w. annotation":[],"tss":[],"first Exon": [], "single Exon": [],"tssAntisense":[],"exon":[],"other":[]}    
    for tsr_file in os.listdir(f'{workingPath}files/TSR/'):
        if tsr_file.endswith('-10.tss.txt'):
            try:
                tsr_frame = pd.read_csv(f'{workingPath}files/TSR/' + tsr_file,  sep='\t')
                dataset = tsr_file.split('-10')[0]
                TSR_dict["Library"].append(dataset)
                tsr_frame_counts = pd.DataFrame(tsr_frame['annotation'].value_counts()) # count annotations from the tsr file    
                # create a frame to merge counts into: this is to avoid errors if one annotations is NOT found
                start_vals = {'traits': ['other','tss','singleExon', 'tssAntisense', 'otherExon', 'firstExon', 'otherExonBidirectional'],
                              'values':['0','0','0','0','0','0','0']}
                startframe = pd.DataFrame(data=start_vals)
                startframe = startframe.set_index('traits')
                sum_frame = startframe.join(tsr_frame_counts)
                sum_frame = sum_frame.astype(np.float64) # convert NaN etc
                sum_frame['sum'] = sum_frame['values'] + sum_frame['count'] # so 0 + 0 = 0 but the order is maintained
                sum_frame = sum_frame.reset_index()    
                other = sum_frame.loc[sum_frame['traits'] == 'other', 'sum'].values[0]
                TSR_dict["other"].append(other)    
                tssAntisense = sum_frame.loc[sum_frame['traits'] == 'tssAntisense', 'sum'].values[0]
                TSR_dict["tssAntisense"].append(tssAntisense)    
                tss = sum_frame.loc[sum_frame['traits'] == 'tss', 'sum'].values[0]
                TSR_dict["tss"].append(tss)    
                singleExon = sum_frame.loc[sum_frame['traits'] == 'singleExon', 'sum'].values[0]
                TSR_dict["single Exon"].append(singleExon)    
                firstExon = sum_frame.loc[sum_frame['traits'] == 'firstExon', 'sum'].values[0]
                TSR_dict["first Exon"].append(firstExon)    
                otherExon = sum_frame.loc[sum_frame['traits'] == 'otherExon', 'sum'].values[0]
                otherExonBidirectional = sum_frame.loc[sum_frame['traits'] == 'otherExonBidirectional', 'sum'].values[0]
                exon = otherExon + otherExonBidirectional
                TSR_dict["exon"].append(exon)    
                TSRs = other + tssAntisense + tss + singleExon + otherExon + otherExonBidirectional + firstExon
                TSR_dict["TSRs w. annotation"].append(TSRs)                
            except(FileNotFoundError):
                exception_type, exception_object, exception_traceback = sys.exc_info()
                splitObjects = str(exception_object).split("'")
                errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                continue
    TSR_annotations_frame = pd.DataFrame(TSR_dict)
    TSR_annotations_frame.to_csv(f'../{species}_summary_TSRannotations.tsv', sep = '\t')
    
def TSRstabAnno(workingPath, errFile):
    TSRFile = f'{workingPath}files/TSR/'
    os.chdir(TSRFile)
    for tsr_file in os.listdir():
        if tsr_file.startswith('tsrs_'):
            if tsr_file.endswith('stable.txt'):
                try:
                    tsr_frame = pd.read_csv(tsr_file,  sep='\t')
                    tsr_frame = tsr_frame[['#tssClusterID']]    
                    # load original .tss file (which are tsrs)
                    all_tsrs_file = tsr_file.split('tsrs_')[-1].split('_stable.txt')[0] + '-10.tss.txt'
                    all_tsrs_frame = pd.read_csv(all_tsrs_file,  sep='\t')    
                    new_frame = pd.merge(all_tsrs_frame,tsr_frame , left_on='#tssClusterID', right_on='#tssClusterID', how='inner')    
                    outputname = tsr_file.split('.txt')[0] + '_anno.txt'    
                    new_frame.to_csv(outputname, sep = '\t')                    
                except(FileNotFoundError):
                    exception_type, exception_object, exception_traceback = sys.exc_info()
                    splitObjects = str(exception_object).split("'")
                    errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                    continue    
    TSR_dict = {"Library":[],"TSRs w. annotation":[],"tss":[],"first Exon": [], "single Exon": [],"tssAntisense":[],"exon":[],"other":[]}
    for tsr_file in os.listdir():
        if tsr_file.startswith('tsrs_'):
            if tsr_file.endswith('_stable_anno.txt'):
                try:
                    tsr_frame = pd.read_csv(tsr_file,  sep='\t')    
                    dataset = tsr_file.split('.txt')[0].split('tsrs_')[1]
                    TSR_dict["Library"].append(dataset)    
                    tsr_frame_counts = pd.DataFrame(tsr_frame['annotation'].value_counts()) # count annotations from the tsr file
                    # create a frame to merge counts into: this is to avoid errors if one annotations is NOT found
                    start_vals = {'traits': ['other','tss','singleExon', 'tssAntisense', 'otherExon', 'firstExon', 'otherExonBidirectional'],
                                  'values':['0','0','0','0','0','0','0']}
                    startframe = pd.DataFrame(data=start_vals)
                    startframe = startframe.set_index('traits')
                    sum_frame = startframe.join(tsr_frame_counts)
                    sum_frame = sum_frame.astype(np.float64) # convert NaN etc
                    sum_frame['sum'] = sum_frame['values'] + sum_frame['count'] # so 0 + 0 = 0 but the order is maintained
                    sum_frame = sum_frame.reset_index()    
                    other = sum_frame.loc[sum_frame['traits'] == 'other', 'sum'].values[0]
                    TSR_dict["other"].append(other)    
                    tssAntisense = sum_frame.loc[sum_frame['traits'] == 'tssAntisense', 'sum'].values[0]
                    TSR_dict["tssAntisense"].append(tssAntisense)    
                    tss = sum_frame.loc[sum_frame['traits'] == 'tss', 'sum'].values[0]
                    TSR_dict["tss"].append(tss)    
                    singleExon = sum_frame.loc[sum_frame['traits'] == 'singleExon', 'sum'].values[0]
                    TSR_dict["single Exon"].append(singleExon)    
                    firstExon = sum_frame.loc[sum_frame['traits'] == 'firstExon', 'sum'].values[0]
                    TSR_dict["first Exon"].append(firstExon)    
                    otherExon = sum_frame.loc[sum_frame['traits'] == 'otherExon', 'sum'].values[0]
                    otherExonBidirectional = sum_frame.loc[sum_frame['traits'] == 'otherExonBidirectional', 'sum'].values[0]
                    exon = otherExon + otherExonBidirectional
                    TSR_dict["exon"].append(exon)    
                    TSRs = other + tssAntisense + tss + singleExon + otherExon + otherExonBidirectional + firstExon
                    TSR_dict["TSRs w. annotation"].append(TSRs)    
                except(FileNotFoundError):
                    exception_type, exception_object, exception_traceback = sys.exc_info()
                    splitObjects = str(exception_object).split("'")
                    errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                    continue    
    stableTSR_annotations_frame = pd.DataFrame(TSR_dict)    
    for tsr_file in os.listdir():
        if tsr_file.startswith('tsrs_'):
            if tsr_file.endswith('_unstable_anno.txt'):
                try:
                    tsr_frame = pd.read_csv(tsr_file,  sep='\t')    
                    dataset = tsr_file.split('.txt')[0].split('tsrs_')[1]
                    TSR_dict["Library"].append(dataset)    
                    tsr_frame_counts = pd.DataFrame(tsr_frame['annotation'].value_counts()) # count annotations from the tsr file    
                    # create a frame to merge counts into: this is to avoid errors if one annotations is NOT found
                    start_vals = {'traits': ['other','tss','singleExon', 'tssAntisense', 'otherExon', 'firstExon', 'otherExonBidirectional'],
                                  'values':['0','0','0','0','0','0','0']}
                    startframe = pd.DataFrame(data=start_vals)
                    startframe = startframe.set_index('traits')
                    sum_frame = startframe.join(tsr_frame_counts)
                    sum_frame = sum_frame.astype(np.float64) # convert NaN etc
                    sum_frame['sum'] = sum_frame['values'] + sum_frame['count'] # so 0 + 0 = 0 but the order is maintained
                    sum_frame = sum_frame.reset_index()
                    other = sum_frame.loc[sum_frame['traits'] == 'other', 'sum'].values[0]
                    TSR_dict["other"].append(other)    
                    tssAntisense = sum_frame.loc[sum_frame['traits'] == 'tssAntisense', 'sum'].values[0]
                    TSR_dict["tssAntisense"].append(tssAntisense)    
                    tss = sum_frame.loc[sum_frame['traits'] == 'tss', 'sum'].values[0]
                    TSR_dict["tss"].append(tss)    
                    singleExon = sum_frame.loc[sum_frame['traits'] == 'singleExon', 'sum'].values[0]
                    TSR_dict["single Exon"].append(singleExon)    
                    firstExon = sum_frame.loc[sum_frame['traits'] == 'firstExon', 'sum'].values[0]
                    TSR_dict["first Exon"].append(firstExon)    
                    otherExon = sum_frame.loc[sum_frame['traits'] == 'otherExon', 'sum'].values[0]
                    otherExonBidirectional = sum_frame.loc[sum_frame['traits'] == 'otherExonBidirectional', 'sum'].values[0]
                    exon = otherExon + otherExonBidirectional
                    TSR_dict["exon"].append(exon)    
                    TSRs = other + tssAntisense + tss + singleExon + otherExon + otherExonBidirectional + firstExon
                    TSR_dict["TSRs w. annotation"].append(TSRs)                
                except(FileNotFoundError):
                    exception_type, exception_object, exception_traceback = sys.exc_info()
                    splitObjects = str(exception_object).split("'")
                    errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                    continue    
    unstableTSR_annotations_frame = pd.DataFrame(TSR_dict)    
    # summarize
    concat_frame = pd.concat([TSR_annotations_frame, stableTSR_annotations_frame, unstableTSR_annotations_frame], axis=0).copy()
    concat_frame.to_csv(f'../{species}_summary_TSR_stability_and_annotations.tsv', sep = '\t')
    
def TSRsummary(workingPath, errFile):
    TSRFile = f'{workingPath}files/TSR/'
    os.chdir(TSRFile)
    tsr_dict = {"Library":[],"putative TSRs":[], "valid TSRs":[],"Fraction Promoter-Distal TSS clusters":[],
                "Fraction of bidirectional TSS clusters":[],"Fraction of stable transcript TSS clusters":[],
                "SS":[],"SU":[], "US":[],"S":[],"UU":[],"U":[]}
    for stats_file in os.listdir():
        if stats_file.endswith('10.stats.txt'):
            try:
                tsr_file = pd.read_csv(stats_file, sep=' ', header=None,)
                print(stats_file)
                # check to see if the RNA tagdir exists for the given exp.
                tagDirName = stats_file.split('-10')[0] + '_RNA'
                if os.path.exists(f'{tagPath}{tagDirName}/tagInfo.txt'):            
                    # get values to populate tha dataframe
                    library = (tsr_file.loc[0][5]).split('/')[8].split('-10')[0].split('10')[0]
                    putative_TSRs = float((tsr_file.loc[4][3]).split('clusters\t')[1])
                    valid_TSRs = float((tsr_file.loc[5][2]).split('clusters\t')[1])
                    distal_TSRs = float((tsr_file.loc[9][4]).split('%')[0])
                    stable_TSRs = float((tsr_file.loc[10][6]).split('%')[0])
                    bidir_TSRs = float((tsr_file.loc[11][5]).split('%')[0])
                    SS_TSRs = float((tsr_file.loc[13][2]).split('%)')[0].split('(')[1])
                    SU_TSRs = float((tsr_file.loc[14][2]).split('%)')[0].split('(')[1])
                    US_TSRs = float((tsr_file.loc[16][2]).split('%)')[0].split('(')[1])
                    S_TSRs = float((tsr_file.loc[15][2]).split('%)')[0].split('(')[1])
                    UU_TSRs = float((tsr_file.loc[17][2]).split('%)')[0].split('(')[1])
                    U_TSRs = float((tsr_file.loc[18][2]).split('%)')[0].split('(')[1])
                    # append to library
                    tsr_dict["Library"].append(library)
                    tsr_dict["putative TSRs"].append(putative_TSRs)
                    tsr_dict["valid TSRs"].append(valid_TSRs)
                    tsr_dict["Fraction Promoter-Distal TSS clusters"].append(distal_TSRs)
                    tsr_dict["Fraction of stable transcript TSS clusters"].append(stable_TSRs)
                    tsr_dict["Fraction of bidirectional TSS clusters"].append(bidir_TSRs)
                    tsr_dict["SS"].append(SS_TSRs)
                    tsr_dict["SU"].append(SU_TSRs)
                    tsr_dict["S"].append(S_TSRs)
                    tsr_dict["UU"].append(UU_TSRs)
                    tsr_dict["US"].append(US_TSRs)
                    tsr_dict["U"].append(U_TSRs)
                else:
                    # get values to populate tha dataframe
                    library = (tsr_file.loc[0][5]).split('/')[8].split('-10')[0].split('10')[0]
                    putative_TSRs = float((tsr_file.loc[3][3]).split('clusters\t')[1])
                    valid_TSRs = float((tsr_file.loc[4][2]).split('clusters\t')[1])
                    distal_TSRs = float((tsr_file.loc[8][4]).split('%')[0])
                    stable_TSRs = float((tsr_file.loc[9][6]).split('%')[0])
                    bidir_TSRs = float((tsr_file.loc[10][5]).split('%')[0])
                    SS_TSRs = float((tsr_file.loc[12][2]).split('%)')[0].split('(')[1])
                    SU_TSRs = float((tsr_file.loc[13][2]).split('%)')[0].split('(')[1])
                    US_TSRs = float((tsr_file.loc[15][2]).split('%)')[0].split('(')[1])
                    S_TSRs = float((tsr_file.loc[14][2]).split('%)')[0].split('(')[1])
                    UU_TSRs = float((tsr_file.loc[16][2]).split('%)')[0].split('(')[1])
                    U_TSRs = float((tsr_file.loc[17][2]).split('%)')[0].split('(')[1])
                    # append to library
                    tsr_dict["Library"].append(library)
                    tsr_dict["putative TSRs"].append(putative_TSRs)
                    tsr_dict["valid TSRs"].append(valid_TSRs)
                    tsr_dict["Fraction Promoter-Distal TSS clusters"].append(distal_TSRs)
                    tsr_dict["Fraction of stable transcript TSS clusters"].append(stable_TSRs)
                    tsr_dict["Fraction of bidirectional TSS clusters"].append(bidir_TSRs)
                    tsr_dict["SS"].append(SS_TSRs)
                    tsr_dict["SU"].append(SU_TSRs)
                    tsr_dict["S"].append(S_TSRs)
                    tsr_dict["UU"].append(UU_TSRs)
                    tsr_dict["US"].append(US_TSRs)
                    tsr_dict["U"].append(U_TSRs)
            except(FileNotFoundError):
                exception_type, exception_object, exception_traceback = sys.exc_info()
                splitObjects = str(exception_object).split("'")
                errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                continue
            except(IndexError):
                exception_type, exception_object, exception_traceback = sys.exc_info()
                errFile.write(f'misformated stats file on {stats_file} giving error:\nKeyError\t:\t{exception_object}')
                continue
            except(KeyError):
                exception_type, exception_object, exception_traceback = sys.exc_info()
                errFile.write(f'misformated stats file on {stats_file} giving error:\nKeyError\t:\t{exception_object}')
                continue
    TSR_stats_frame = pd.DataFrame(tsr_dict)
    TSR_stats_frame.to_csv(f'{species}_summaryTSRstats.tsv', sep = '\t')    
    summaryStats_frame = pd.merge(TSR_stats_frame, stability_dict_frame, left_on='Library', right_on='species', how='outer')
    summaryStats_frame = pd.merge(summaryStats_frame, stability_dictTSR_frame, left_on='Library', right_on='species', how='outer')
    summaryStats_frame = pd.merge(summaryStats_frame, TSR_annotations_frame, left_on='Library', right_on='Library', how='outer')    
    summaryStats_frame = summaryStats_frame[['Library', 'putative TSRs', 'valid TSRs', 'stableTSRs',
           'unstableTSRs','stableTSSs', 'unstableTSSs',
           'Fraction Promoter-Distal TSS clusters',
           'Fraction of bidirectional TSS clusters',
           'Fraction of stable transcript TSS clusters', 'S','SS', 'SU', 'U','US',
           'UU',  'TSRs w. annotation', 'tss', 'first Exon',
           'single Exon', 'tssAntisense', 'exon', 'other']]    
    summaryStats_frame = summaryStats_frame.sort_values(by=['Library'])
    summaryStats_frame.to_csv(f'../{species}_summary_allTSRsDetails.tsv', sep = '\t')
    
def mappingStats(workingPath, errFile):
    TSRFile = f'{workingPath}files/TSR/'
    mappingPath = f'{workingPath}files/mappingStats/'
    os.chdir(f'{TSRFile}../')
    mappingStats_dict = {"Library":[],"Reads":[], "Adapter reads":[],
                         "Aligned 0 times":[],"Aligned 1 time":[],
                         "Aligned >1 times":[], "Adapters %":[],
                         "Aligned 0 times %":[],"Aligned 1 time %":[],
                         "Aligned >1 times %":[],"Alignment rate":[]}
    for mapping_file in os.listdir(mappingPath):
        if mapping_file.endswith('_mappingstats.txt'):
            try:
                mapping_frame = pd.read_csv(mappingPath + mapping_file,  sep='\t')
                # write alternate ways to perform pulling stats from csRNA, sRNA, and totalRNA.
                if ('_totalRNA' in mapping_file) or ('ChIP' in mapping_file):                
                    library = mapping_file.split('.fastq')[0]
                    reads = (mapping_frame.loc[0][0]).split(' ')[2]
                    aligned_0 = (mapping_frame.loc[10][0]).split(' (')[0].split(' ')[-1]
                    aligned_0percent = (mapping_frame.loc[10][0]).split('(')[1].split(')')[0]                
                    aligned_1 = (mapping_frame.loc[11][0]).split(' (')[0].split(' ')[-1]
                    aligned_1percent = (mapping_frame.loc[11][0]).split('(')[1].split(')')[0]    
                    aligned_more = (mapping_frame.loc[12][0]).split(' (')[0].split(' ')[-1]
                    aligned_morePercent = (mapping_frame.loc[12][0]).split('(')[1].split(')')[0]            
                    rate = (mapping_frame.loc[13][0]).split(' ')[0]        
                    # also read out adapter dimers
                    species = mapping_file.split('_')[0] + '_' + mapping_file.split('_')[1].split('-')[0]                
                    if ('_totalRNA' in mapping_file) or ('_ChIP' in mapping_file):
                        seqType = 'totalRNA'
                    elif '_sRNA' in mapping_file:
                        seqType = 'sRNA'
                    elif '_csRNA' in mapping_file:
                        seqType = 'csRNA'
                    adapter_dimers_reads = 'NA'
                    adapters_percent = 'NA'                
                    mappingStats_dict["Library"].append(library)
                    mappingStats_dict["Reads"].append(reads)
                    mappingStats_dict["Adapter reads"].append(adapter_dimers_reads)     
                    mappingStats_dict["Aligned 0 times"].append(aligned_0)
                    mappingStats_dict["Aligned 1 time"].append(aligned_1)
                    mappingStats_dict["Aligned >1 times"].append(aligned_more)
                    mappingStats_dict["Adapters %"].append(adapters_percent)
                    mappingStats_dict["Aligned 0 times %"].append(aligned_0percent)
                    mappingStats_dict["Aligned 1 time %"].append(aligned_1percent)
                    mappingStats_dict["Aligned >1 times %"].append(aligned_morePercent)
                    mappingStats_dict["Alignment rate"].append(rate)    
                elif ('_csRNA' in mapping_file) or ('_sRNA' in mapping_file):       
                    library = mapping_file.split('.fastq')[0]
                    reads = (mapping_frame.loc[0][0]).split(' ')[0]
                    aligned_0 = (mapping_frame.loc[2][0]).split(' (')[0].split(' ')[-1]
                    aligned_0percent = (mapping_frame.loc[2][0]).split('(')[1].split(')')[0]        
                    aligned_1 = (mapping_frame.loc[3][0]).split(' (')[0].split(' ')[-1]
                    aligned_1percent = (mapping_frame.loc[3][0]).split('(')[1].split(')')[0]        
                    aligned_more = (mapping_frame.loc[4][0]).split(' (')[0].split(' ')[-1]
                    aligned_morePercent = (mapping_frame.loc[4][0]).split('(')[1].split(')')[0]        
                    rate = (mapping_frame.loc[5][0]).split(' ')[0]
                    # also read out adapter dimers
                    species = mapping_file.split('_')[0] + '_' + mapping_file.split('_')[1].split('-')[0]
                    if '_totalRNA' in mapping_file:
                        seqType = 'totalRNA'
                    elif '_sRNA' in mapping_file:
                        seqType = 'sRNA'
                    elif '_csRNA' in mapping_file:
                        seqType = 'csRNA'
                    trimmed_lengths_file = dataPath + species + f'/fastq/{seqType}/' + mapping_file.split('_mappingstats.txt')[0] + '.fastq.gz.lengths'
                    trimmed_lengths_frame = pd.read_csv(trimmed_lengths_file,  sep='\t')
                    adapter_dimers_reads = trimmed_lengths_frame.loc[0][1]
                    adapters_percent = round(float((trimmed_lengths_frame.loc[0][2]).split('%')[0]),2)                    
                    mappingStats_dict["Library"].append(library)
                    mappingStats_dict["Reads"].append(reads)
                    mappingStats_dict["Adapter reads"].append(adapter_dimers_reads)     
                    mappingStats_dict["Aligned 0 times"].append(aligned_0)
                    mappingStats_dict["Aligned 1 time"].append(aligned_1)
                    mappingStats_dict["Aligned >1 times"].append(aligned_more)
                    mappingStats_dict["Adapters %"].append(adapters_percent)
                    mappingStats_dict["Aligned 0 times %"].append(aligned_0percent)
                    mappingStats_dict["Aligned 1 time %"].append(aligned_1percent)
                    mappingStats_dict["Aligned >1 times %"].append(aligned_morePercent)
                    mappingStats_dict["Alignment rate"].append(rate)
            except(FileNotFoundError):
                exception_type, exception_object, exception_traceback = sys.exc_info()
                splitObjects = str(exception_object).split("'")
                errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
                continue    
    mappingStats_dict_frame = pd.DataFrame(mappingStats_dict)
    mappingStats_dict_frame = mappingStats_dict_frame.sort_values(by=['Library'])
    mappingStats_dict_frame.to_csv(f'{species}_summary_mappingStats.tsv', sep = '\t',index=False)
    
def tagDirStats(workingPath, errFile):
    # Now that we have each sheet seperately populated with all data, combine the common lines between them
    try:
        os.chdir(f'{workingPath}files')
        tagDirStats = pd.read_csv(f'{species}_summary_mappingStats.tsv', sep = '\t')
        mappingStats = pd.read_csv(f'{species}_summary_tagDirsStats.tsv', sep = '\t')        
        with open(f'{species}_summary_mappingAndtagDirstats.tsv','w') as sumFile:
            headerStr = 'Name\tReads\tAdapter Reads\tAligned 0 Times\tAligned 1 time\tAligned>1 times\tAdapters %\tAligned 0 times %\tAligned 1 times %\tAligned>1 times%\tAlignment rate\tName again\tUnique Positions\tTotal Tags\tTags Per BP\tTag Positions in Analysis\tAverage Tags Per Position\tAverage Fragment GC Content\n'
            sumFile.write(headerStr)
            for MSindex,MSrow in tagDirStats.iterrows():
                # split the library name by the -r(1|2) value and re-append the removed value to the name. use this to search the other df
                if 'r1' in MSrow['Library']:
                    file = re.split(r"-r[1|2|3|4|5]",MSrow['Library'])[0]
                    file = file + '-r1'
                elif 'r2' in MSrow['Library']:
                    file = re.split(r"-r[1|2|3|4|5]",MSrow['Library'])[0]
                    file = file + '-r2'
                elif 'r3' in MSrow['Library']:
                    file = re.split(r"-r[1|2|3|4|5]",MSrow['Library'])[0]
                    file = file + '-r3'
                elif 'r4' in MSrow['Library']:
                    file = re.split(r"-r[1|2|3|4|5]",MSrow['Library'])[0]
                    file = file + '-r3'
                elif 'r5' in MSrow['Library']:
                    file = re.split(r"-r[1|2|3|4|5]",MSrow['Library'])[0]
                    file = file + '-r3'        
                for TDSindex, TDSrow in mappingStats.iterrows():
                    if file == TDSrow['File']:
                        MSrow = MSrow.to_string().split('\n')
                        TDSrow = TDSrow.to_string().split('\n')
                        newRow = []
                        count = 0
                        for val in MSrow:
                            if count >= 1:
                                newRow.append(val.split(' ')[-1])
                            else:
                                count+=1
                        count = 0
                        for val in TDSrow:
                            if count >=1:
                                newRow.append(val.split(' ')[-1])
                            else:
                                count+=1
                        rowStr = "\t".join(newRow)
                        rowStr = file+'\t'+rowStr
                        sumFile.write(rowStr+'\n')
        sumFile.close()
    except(FileNotFoundError):
        exception_type, exception_object, exception_traceback = sys.exc_info()
        splitObjects = str(exception_object).split("'")
        errFile.write(f'{splitObjects[0]}\t:\t{splitObjects[1]}\n')
        print(f'missing either {species}_summary_mappingStats.tsv or {species}_summary_tagDirsStats.tsv')
        
def genericFunc(cmd):
    print(cmd)
    proc = subprocess.Popen([cmd],shell=True, stdout=subprocess.PIPE)
    proc.communicate()