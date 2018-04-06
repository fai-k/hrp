#Scan for telomere split reads
#Required text file: input-namelocat.txt which contains "isolate name" (tab or space) "file location" in the same line
#Required: pysam 0.11.2.2, pysamstats 1.0.1
from __future__ import print_function
import numpy as np
import pysam
import pysamstats as pss
import csv
import argparse
import re

#Paths and Parameters#
parser = argparse.ArgumentParser(description='Extract out-of-place soft-clipped telomere reads')
parser.add_argument('-i', '--infile', type=str, metavar='<inputfile>', help='input text file containing "isolate name"<tab>"file location" one isolate per line')
parser.add_argument('-o', '--outfold', type=str, metavar='<outfolder>', help='output folder for "isolate name"-telo.bam')
args = parser.parse_args()
inputfile = args.infile
outfolder = args.outfold
forward=re.compile('A{0,1}G{2,3}T{2,4}C{0,1}')
reverse=re.compile('A{2,4}C{2,3}T{0,1}G{0,1}')
notelo = [['Pf3D7_01_v3',29510,614893],
	['Pf3D7_02_v3',25232,923648],
	['Pf3D7_03_v3',36965,1038254],
	['Pf3D7_04_v3',28706,1180226],
	['Pf3D7_05_v3',20929,1342964],
	['Pf3D7_06_v3',653,1382627],
	['Pf3D7_07_v3',20307,1426234],
	['Pf3D7_08_v3',21361,1443449],
	['Pf3D7_09_v3',20080,1503336],
	['Pf3D7_10_v3',28490,1649948],
	['Pf3D7_11_v3',24160,2035886],
	['Pf3D7_12_v3',16973,2248962],
	['Pf3D7_13_v3',21364,2892340],
	['Pf3D7_14_v3',1393,3291501]]
###

samples=[]
with open(inputfile,'r') as f:
    for line in csv.reader(f, delimiter='\t'):
        samples.append(line)

for i in range(0,len(samples)):
	bam = pysam.AlignmentFile(samples[i][1],'rb')
	get = pysam.AlignmentFile('%s/%s-telo.bam' % (outfolder, samples[i][0]),'wb', template=bam)
	for j in range(0,len(notelo)):
		for read in bam.fetch(notelo[j][0],notelo[j][1],notelo[j][2]):
			cig = read.cigartuples
			if cig != None :
				start = 0
				end = 0
				for k in range(0,len(cig)):
					start = end
					end = end + cig[k][1]
					if cig[k][0] == 4 :
						countfor = forward.findall(read.query_sequence,start,end)
						countrev = reverse.findall(read.query_sequence,start,end)
						if (len(countfor) > 2) or (len(countrev) > 2):
							get.write(read)
							break
	bam.close()
	get.close()
	print('Done: %s' % samples[i][0])
#All reads with telomere repeat match > 2 in softclipped region is collected in "isolate name"-telo.bam
