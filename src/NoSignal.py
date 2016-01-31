import fasta
import sys
import random
import re
import os
import yaml
import atgcDistribution
import SeqGenUtils

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



def CreateNoSignalSequences(NegSequences, NegHeaders, NumSeqsToGenerate, 
	                           SeqLength, OutputFileName):
	NumNegSequences = len(NegSequences)
	NegSeqDict = dict();

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

		selectedNucs = allNucs[startPoint:endPoint];
		l = list(str(selectedNucs))
		random.shuffle(l)
		shuffledNucs = ''.join(l);
		NegSeqDict[NegHeaders[index]] = shuffledNucs

	print type(NegSeqDict)
	return NegSeqDict;



def GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, SeqLength, OutFileName):
	NegSequences, NegHeaders = CreateNegDict(NegativeFileName);
	return CreateNoSignalSequences(NegSequences, NegHeaders, NumSeqsToGenerate, SeqLength, OutFileName)


if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)


	NegativeFileName = confMap["sequence"]["fastaFile"]
	NumSeqsToGenerate = int(confMap["sequence"]["numSeq"])
	SeqLength = int(confMap["sequence"]["seqLen"])
	OutFileName = confMap["sequence"]["outFastaFile"]

	NegSeqDict = GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, SeqLength, OutFileName);
	SeqGenUtils.WriteSeqDictToFile(NegSeqDict, OutFileName);
