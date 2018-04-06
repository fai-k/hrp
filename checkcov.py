#Check coverage for suspected HRP2/3 deletion
#Required text file: input-namelocat.txt which contains "isolate name" (tab or space) "file location" in the same line
#Required: pysam 0.11.2.2, pysamstats 1.0.1
from __future__ import print_function
import numpy as np
import pysam
import pysamstats as pss
import csv
import argparse

#Paths and Parameters#
parser = argparse.ArgumentParser(description='Check normalized coverage of exons of HRP2, HRP3 and neighboring genes')
parser.add_argument('-r', '--ref', type=str, metavar='<reffile>', help='reference file as fasta')
parser.add_argument('-i', '--infile', type=str, metavar='<inputfile>', help='input text file containing "isolate name"<tab>"file location<tab>"median coverage" one isolate per line')
parser.add_argument('-o', '--outfold', type=str, metavar='<outfolder>', help='output folder for HRPcov_"isolate name".txt')
args = parser.parse_args()
reffile = args.ref
inputfile = args.infile
outfolder = args.outfold
cutoff = 0.1	#the normalized coverage cut point for flagging an exon as suspected deletion
cutoffHRP = 0.3 	#the normalized coverage cut point for flagging HRP exon as suspected deletion
exons = [['Pf3D7_08_v3',1358314,1359142,'PF3D7_0831600','Chr8_CLAG8'],
	['Pf3D7_08_v3',1359287,1359859,'PF3D7_0831600','Chr8_CLAG8'],
	['Pf3D7_08_v3',1361249,1362015,'PF3D7_0831600','Chr8_CLAG8'],
	['Pf3D7_08_v3',1362197,1363618,'PF3D7_0831600','Chr8_CLAG8'],
	['Pf3D7_08_v3',1365467,1367506,'PF3D7_0831700','Chr8_HSP70x'],
	['Pf3D7_08_v3',1374236,1375084,'PF3D7_0831800','Chr8_HRP2'],
	['Pf3D7_08_v3',1383349,1384212,'PF3D7_0832000','Chr8_stevor'],
	['Pf3D7_08_v3',1392407,1393135,'PF3D7_0832200','Chr8_PHISTa-like'],
	['Pf3D7_13_v3',2833204,2834322,'PF3D7_1372000','Chr13_PHISTa'],
	['Pf3D7_13_v3',2837313,2839058,'PF3D7_1372100','Chr13_PHISTb'],
	['Pf3D7_13_v3',2840727,2841485,'PF3D7_1372200','Chr13_HRP3'],
	['Pf3D7_13_v3',2846002,2846538,'PF3D7_1372300','Chr13_PHIST'],
	['Pf3D7_13_v3',2850891,2853482,'PF3D7_1372400','Chr13_ACS4']]
###

samples=[]
with open(inputfile,'r') as f:
    for line in csv.reader(f, delimiter='\t'):
        samples.append(line)

suspectdel=[]
for i in range(0,len(samples)):
	bam = pysam.AlignmentFile(samples[i][1],'rb')
	#Coverage of HRP and neighboring gene#
	exoncov=[]
	for k in range(0,len(exons)):
		ecoverage = pss.load_coverage(bam, reffile, pad = True, truncate = True, chrom = exons[k][0], start = exons[k][1], end = exons[k][2])
		exonmean = np.mean(ecoverage.reads_all)
		exoncov.append([exonmean, exonmean/float(samples[i][2])])
	with open('%s/HRPcov_%s.txt' % (outfolder, samples[i][0]),'w') as f:
		f.write('chromosome\tstart\tend\tgeneID\tname\tcoverage\tnormalizedcov\n')
		for m in range(0,len(exons)):
			f.write('%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n' % (exons[m][0], exons[m][1], exons[m][2], exons[m][3], exons[m][4], exoncov[m][0], exoncov[m][1]))
	#Flag suspected deletion#
	for n in range(0,len(exoncov)):
		if exoncov[n][1] < cutoff:
			suspectdel.append(samples[i])
			break
		else:
			if ((n == 5) or (n == 10)) and (exoncov[n][1] < cutoffHRP):
				suspectdel.append(samples[i])
				break
	print('Done: %s\n' % samples[i][0])
with open('%s/suspectdel.txt' % outfolder,'w') as f:
	for p in range(0,len(suspectdel)):
		f.write('%s\t%s\t%.3f\n' % (suspectdel[p][0], suspectdel[p][1], suspectdel[p][2]))


#Normalized coverage of exons in HRP2 and HRP3 neigboring genes are listed in HRPcov_IsolateName.txt
#Isolate names of the suspected HRP2 or HRP3 deletion are listed in suspectdel.txt
