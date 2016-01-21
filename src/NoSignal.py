import fasta
import sys
import random
import re
import os
import yaml
import atgcDistribution

indexArr = []

def CreateNegDict(NegativeFileName):
	NegSequences = []
	NegHeaders = []
	NegativeFile = open(NegativeFileName, "r")
	for record in fasta.fasta_itr(NegativeFileName):						#Create parallel arrays for negative sequence/headers
		sequence = record.sequence
		header = record.header
		
		NegSequences.append(sequence)
		NegHeaders.append(header)
	return NegSequences, NegHeaders;

def GetSequenceToShuffle(NegSequences, targetSeqLength):
	gotLength = False;
	index = 0;
	allNucs = ""

	while(gotLength != True):											#Pick random negative sequence, ensuring it hasn't been picked before
		index = random.randint(0, (len(NegSequences)-1))				
		if index in indexArr:
			while index in indexArr:
				index = random.randint(0, (len(NegSequences)-1))
		indexArr.append(index)
		allNucs = NegSequences[index]
		seqLength = len(allNucs)
		if seqLength >= targetSeqLength:											#Ensure that negative sequence has enough nucleotides
			gotLength = True

	return allNucs, index;

def WriteSeqToFile(OutputFile, NegHeaders, index, selectedNucs):
	OutputFile.write(">")												#Write to output file
	OutputFile.write(NegHeaders[index])
	OutputFile.write("\n")
	OutputFile.write(selectedNucs)
	OutputFile.write("\n")

def GenerateNoSignalSequences(NegSequences, NegHeaders, NumSeqsToGenerate, 
	                           SeqLength, OutputFileName):

	
	NumNegSequences = len(NegSequences)
	OutputFile = open(OutputFileName , "w")

	for i in range(0, NumSeqsToGenerate):
		allNucs, index = GetSequenceToShuffle(NegSequences, SeqLength)
		allNucsLength = len(allNucs)
		#print "All Nucs Length: " + str(allNucsLength)
		startPoint = random.randint(0, allNucsLength - 1)
		#print "Start Point: " + str(startPoint)
		if startPoint + SeqLength < allNucsLength:
			endPoint = startPoint + SeqLength;
		else:
			startPoint = 0;
			endPoint = startPoint + SeqLength;

		selectedNucs = ''.join(random.sample(allNucs[startPoint:endPoint], (endPoint - startPoint)))
		WriteSeqToFile(OutputFile, NegHeaders, index, selectedNucs)		

	OutputFile.close()

def GetConf(configFile):
	f = open(configFile);
	confMap = yaml.load(f)
	f.close();
	return confMap;

if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	confMap = GetConf(configFile)


	NegativeFileName = confMap["sequence"]["fastaFile"]
	NumSeqsToGenerate = int(confMap["sequence"]["numSeq"])
	SeqLength = int(confMap["sequence"]["seqLen"])
	OutFileName = confMap["sequence"]["outFastaFile"]

	print "ATGC Distribution for: ", NegativeFileName;
	atgcDistribution.printATGCDistribution(NegativeFileName);
	NegSequences, NegHeaders = CreateNegDict(NegativeFileName);
	GenerateNoSignalSequences(NegSequences, NegHeaders, 
		                       NumSeqsToGenerate, SeqLength, OutFileName)
	print "ATGC Distribution for: ", OutFileName;
	atgcDistribution.printATGCDistribution(OutFileName);
