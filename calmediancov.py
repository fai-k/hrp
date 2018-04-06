#Calculate median coverage of each bam file
#Required text file: input-namelocat.txt which contains "isolate name" (tab or space) "file location" in the same line
#Required: pysam 0.11.2.2, pysamstats 1.0.1from __future__ import print_function
import numpy as np
import pysam
import pysamstats as pss
import csv
import argparse

#Paths and Parameters#
parser = argparse.ArgumentParser(description='Calculate median coverage, output to the same folder as input files')
parser.add_argument('-r', '--ref', type=str, metavar='<reffile>', help='reference file as fasta')
parser.add_argument('-i', '--infile', type=str, metavar='<inputfile>', help='input text file containing "isolate name"<tab>"file location" one isolate per line')
args = parser.parse_args()
reffile = args.ref
inputfile = args.infile
core_win = [['Pf3D7_01_v3',92901,575900],
	['Pf3D7_02_v3',105801,862500],
	['Pf3D7_03_v3',70631,1003060],
	['Pf3D7_09_v3',79101,1473560],
	['Pf3D7_10_v3',68971,1571815],
	['Pf3D7_11_v3',110001,2003320],
	['Pf3D7_13_v3',74414,2791900]]
###

samples=[]
with open(inputfile,'r') as f:
    for line in csv.reader(f, delimiter='\t'):
        samples.append(line)

cov_median=[]
for i in range(0,len(samples)):
	bam = pysam.AlignmentFile(samples[i][1],'rb')
	#Set reference#
	for j in range(0,len(core_win)):
		if j == 0:
			rcoverage = pss.load_coverage_binned(bam, reffile, chrom = core_win[j][0], start = core_win[j][1], end = core_win[j][2])
		else:
			rcoverage = np.append(rcoverage,pss.load_coverage_binned(bam, reffile, chrom = core_win[j][0], start = core_win[j][1], end = core_win[j][2]))
	cov=[]	
	for k in range(0,len(rcoverage)):
		if rcoverage[k][2] >= 20:
			cov.append(rcoverage[k][3])
	cov_median.append(np.median(cov)*100/300)
	print('Done: %s' % samples[i][0])
with open('%s.CovMed' % inputfile,'w') as f:
	for p in range(0,len(samples)):
		f.write('%s\t%s\t%.3f\n' % (samples[p][0], samples[p][1], cov_median[p]))

#Median coverage is listed as the third column in "inputfile".CovMed with the same folder as input file
