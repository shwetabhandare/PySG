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


def GetRandomNucleotide(items): 
    return choice("".join(x * y for x, y in items))

def GenerateNoSignalFromDirichlet(seqDistDirichletList, seqBackGroundDict,
	SeqLength):
	seq = ""
	keys = list(seqBackGroundDict.keys())
	NegSeqDict = dict()
	for seqNum, atgcDistribution in  enumerate(seqDistDirichletList):
		a = atgcDistribution
		b = list()
		b = [(keys[0], a[0] * 100), (keys[1], a[1] * 100), (keys[2], a[2] * 100), 
		     (keys[3], a[3] * 100)]
		print b
		for count in range(SeqLength):
			seq += GetRandomNucleotide(b)
		key = "NoSignal_" + str(seqNum);
		NegSeqDict[key] = seq;

	for key, value in NegSeqDict.iteritems():
		print key, value
	return NegSeqDict;


		


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
	return NegSeqDict;

def GetDirichletDistribution(seqBackGroundDict, NumSeqsToGenerate):
	s = np.random.dirichlet((seqBackGroundDict["A"] * 100, seqBackGroundDict["C"] * 100,
		                      seqBackGroundDict["T"] * 100, seqBackGroundDict["G"] * 100), 
	                         NumSeqsToGenerate);
	return s;

def GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, SeqLength, OutFileName):
	NegSequences, NegHeaders = CreateNegDict(NegativeFileName);
	return CreateNoSignalSequences(NegSequences, NegHeaders, NumSeqsToGenerate, SeqLength, OutFileName)


if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)
	NegativeFileName = "";
	NegSeqDict  = dict();

	NumSeqsToGenerate = int(confMap["sequence"]["numSeq"])

	if confMap["sequence"].get("fastaFile"):
		NegativeFileName = confMap["sequence"]["fastaFile"]
		

	SeqLength = int(confMap["sequence"]["seqLen"])
	OutFileName = confMap["sequence"]["outFastaFile"]

	if NegativeFileName != "":
		NegSeqDict = GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, SeqLength, OutFileName);
	else:
		seqBackGroundDict = confMap["sequence"]["seqBackGround"]
		seqDistList = GetDirichletDistribution(seqBackGroundDict, 
			                                    NumSeqsToGenerate);
		NegSeqDict = GenerateNoSignalFromDirichlet(seqDistList, seqBackGroundDict, 
			SeqLength);	

	SeqGenUtils.WriteSeqDictToFile(NegSeqDict, OutFileName);
