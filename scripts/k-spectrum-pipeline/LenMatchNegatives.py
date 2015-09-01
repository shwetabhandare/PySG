'''
LenMatchNegatives.py
Requires Python 2.7 and fasta.py in the same directory

Creates length-matched negative sequences for each sequence in a list of positive transcripts

Takes two command line arguments:
Positive File Name (to match to, in fasta format)		Negative File Name (to draw sequences from, in fasta format)

Outputs file called NON + Positive File Name
'''
import fasta
import sys
import random
import re
import os

#Prepare and open files
PositiveFileName = sys.argv[1]
NegativeFileName = sys.argv[2]
PosFileBaseName = os.path.basename(PositiveFileName);
print PosFileBaseName
OutputFileName = "NON" + PosFileBaseName 
print OutputFileName;
PositiveFile = open(PositiveFileName, "r")
NegativeFile = open(NegativeFileName, "r")
OutputFile = open(OutputFileName , "w")

PosLengths = []
for record in fasta.fasta_itr(PositiveFileName):						#Remember lengths of positive sequences
	seqLength = len(record.sequence)
	PosLengths.append(seqLength)

NegSequences = []
NegHeaders = []
for record in fasta.fasta_itr(NegativeFileName):						#Create parallel arrays for negative sequence/headers
	sequence = record.sequence
	header = record.header
	
	NegSequences.append(sequence)
	NegHeaders.append(header)

indexArr = []
#Iterate through lengths in positive length array, picking a random negative sequence 
#and picking LENGTH nucleotides at random start point
for number in PosLengths:										
	gotLength = False
	seqLength = 0
	index = 0
	allNucs = ""
	while(gotLength != True):											#Pick random negative sequence, ensuring it hasn't been picked before
		index = random.randint(0, (len(NegSequences)-1))				
		if index in indexArr:
			while index in indexArr:
				index = random.randint(0, (len(NegSequences)-1))
		indexArr.append(index)
		allNucs = NegSequences[index]
		seqLength = len(allNucs)
		if seqLength > number:											#Ensure that negative sequence has enough nucleotides
			gotLength = True
		else:
			print "Did not get length";
	startPoint = random.randint(0, seqLength - (number+1))
	counter = 0
	selectedNucs = ""
	while (counter < number):											#Select LENGTH nucleotides
		selectedNucs = selectedNucs + allNucs[startPoint + counter]
		counter = counter + 1

	OutputFile.write(">")												#Write to output file
	OutputFile.write(NegHeaders[index])
	OutputFile.write("\n")
	OutputFile.write(selectedNucs)
	OutputFile.write("\n")

PositiveFile.close()
NegativeFile.close()
OutputFile.close()
