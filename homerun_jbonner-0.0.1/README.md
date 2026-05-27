Homerun
Supplementary Software for Homer
Duttke Lab


INTRODUCTION

Homerun automates several functions of Homer, STAR, and HISAT2 for ease of use and time efficiency.
Resulting files are neatly organized and can then be utilized by Homerun for quality control and statistical analysis.

A full Homerun project would involve calling the software for 3 steps:
1. Tag directory creation
2. Quality control analysis
3. Statistical analysis

Tag directories must be created before quality control can be conducted on them, and the quality of the tag directories
should generally be confirmed before performing further analysis.

DEPENDENCIES

STAR, HISAT2, HOMER, pandas

USING HOMERUN

The primary Homerun command uses this format:

python Homerun.py -p [species] -f [fastq directory] -g [genome directory] -m [mode] -t [step]

species:		  the species being analyzed
fastq directory:  the directory where the fastq files to be analyzed are located - only required for tagdir step
genome directory: the directory where the appropriate species' genome file is located - only required for tagdir step
mode:			  "star" or "hisat"
step:			  "tagdir", "qc", or "stat"

Optional fields:
-w [working path]: The parent directory that will contain all outputs of the Homerun project
				   Default value: current directory when Homerun is called
-n [n-tag count]:  n-tag number
				   Default value: 7

An example of a command as a user might write it on the command line:

python ../homerun/Homerun.py -p Homo_sapiens -f ../fastq/project_1/ -g ../genome/Homo_sapiens/ -n 10 -m star -t tagdir

The user-provided fastq and genome files are copied into the new directories that Homerun makes for the project.
This preserves an archive of all input data, but users should remain aware of this when space is a scarce resource.


TROUBLESHOOTING

If tag directories fail during a STAR run but no useful error is logged, try running again with more RAM.

Notes for Jake below this point:

TO DO LIST:
ls *gz -> .tsv
TSS and TSR in tagdir step
exclude median plots from QC
optional command: check tsv
					3' adapter
					-t all
					
					
length plots good for showing enrichment (show %sRNA 20-24nt, save that stat)
A-plot % difference from mean, save that stat (K562 A-like)
frac vs pos plot

stats pdf:
,
total reads (before trimming)
% adapter
% unique/nonaligned
, tag dir
Position in analysis
reads per bp
, TSR
TSR count
% stable
% bidirectional
,
A-plot/K562
% sRNA
frac vs pos plot
enrichment difference
, layout
		A-plot	length
sRNA	plot	plot
csRNA	plot	plot
	PCA, cluster (subsample 20,000)
	
checkTSVqc