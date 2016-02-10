import fasta
import sys
import random
import re
import os
import yaml
import atgcDistribution
import SeqGenUtils
from random import choice
from numpy import *
import numpy as np

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
	return NegSeqDict;

def GenerateNoSignalWithDirichlet(confMap, NumSeqsToGenerate, SeqLength):
	seqBackGroundDict = confMap["sequence"]["nosignal"]["seqBackGround"]
	seqDistList = SeqGenUtils.GetDirichletDistribution(seqBackGroundDict, 
			                                  NumSeqsToGenerate);
	NegSeqDict = SeqGenUtils.GenerateNoSignalFromDirichlet(seqDistList, seqBackGroundDict, 
			SeqLength);	
	return NegSeqDict


def GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, SeqLength, OutFileName):
	NegSequences, NegHeaders = CreateNegDict(NegativeFileName);
	print SeqLength
	return CreateNoSignalSequences(NegSequences, NegHeaders, NumSeqsToGenerate, SeqLength, OutFileName)

def CreateNoSignalDict(confMap):
	NegativeFileName = "";
	NegSeqDict  = dict();
	SeqInfo  = dict();
		
	SeqInfo = SeqGenUtils.GetSequenceInfoFromConfMap(confMap);

	OutFileName = SeqGenUtils.GetNoSignalOutFileName(confMap);
	NegativeFileName = SeqInfo['inputName'];

	print SeqInfo, OutFileName;

	if NegativeFileName != "":
		NegSeqDict = GenerateNoSignalSequences(SeqInfo['inputName'], SeqInfo['numSeq'], 
			                                    SeqInfo['seqLen'], OutFileName);
	else:
	 	NegSeqDict = GenerateNoSignalWithDirichlet(confMap, SeqInfo['numSeq'], SeqInfo['seqLen']);	

	return NegSeqDict, OutFileName;

def CreateNoSignalFastaFile(configFile):
	confMap = SeqGenUtils.GetConf(configFile)
	NegSeqDict = CreateNoSignalDict(confMap);
	OutFileName = SeqGenUtils.GetNoSignalOutFileName(confMap);
	SeqGenUtils.WriteSeqDictToFile(NegSeqDict, OutFileName);	

if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]

