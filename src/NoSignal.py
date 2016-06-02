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
import subprocess

from altschulEriksonDinuclShuffle import dinuclShuffle


indexArr = []
uShuffle = "/projects/bhandare/workspace/PySG/src/ushuffle/ushuffle"

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
		#l = list(str(selectedNucs))
		#random.shuffle(l)
		#shuffledNucs = ''.join(l);
		NegSeqDict[NegHeaders[index]] = selectedNucs;
	return NegSeqDict;

def GenerateNoSignalWithDirichlet(confMap, alpha, NumSeqsToGenerate, SeqLength, signalFlag):
	seqBackGroundDict = confMap["sequence"]["nosignal"]["seqBackGround"]
	seqDistList = SeqGenUtils.GetDirichletDistribution(seqBackGroundDict, alpha, NumSeqsToGenerate);
	NegSeqDict = SeqGenUtils.GenerateNoSignalFromDirichlet(seqDistList, seqBackGroundDict, 
			SeqLength, signalFlag);	
	return NegSeqDict


def GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, SeqLength, OutFileName):
	NegSequences, NegHeaders = CreateNegDict(NegativeFileName);
	return CreateNoSignalSequences(NegSequences, NegHeaders, NumSeqsToGenerate, SeqLength, OutFileName)

def GetShuffledSequence(sequence):
	k = 3;
	shuffled_seq = subprocess.check_output([uShuffle, "-s", sequence, "-k", str(k)])
	shuffled_seq= shuffled_seq[:-1]
	return shuffled_seq;


def ShuffleToCreateNoSignalSequences(PosSeqDict, configFile):
	NegSeqDict = dict();
	for seq_id, sequence in PosSeqDict.iteritems():
		shuffledSeq = GetShuffledSequence(sequence)
		#shuffledSeq = dinuclShuffle(sequence);
		#print "Original : ", sequence;
		#print "Shuffled: ", shuffledSeq;
		NegSeqDict[seq_id] = shuffledSeq;

	confMap = SeqGenUtils.GetConf(configFile)
	NoSignalFileName = confMap["sequence"]["nosignal"]["outNoSignalFastaFile"]
	SeqGenUtils.WriteSeqDictToFile(NegSeqDict, NoSignalFileName);

	return NegSeqDict;

def CreateNoSignalDict(confMap, signalFlag):
	NegativeFileName = "";
	NegSeqDict  = dict();
	SeqInfo  = dict();
		
	SeqInfo = SeqGenUtils.GetSequenceInfoFromConfMap(confMap);

	OutFileName = SeqGenUtils.GetNoSignalOutFileName(confMap);
	NegativeFileName = SeqInfo['inputName'];

	#print SeqInfo, OutFileName

	if NegativeFileName != "":
		NegSeqDict = GenerateNoSignalSequences(NegativeFileName, SeqInfo['numSeq'], 
			                                    SeqInfo['seqLen'], OutFileName);
	else:
	 	NegSeqDict = GenerateNoSignalWithDirichlet(confMap, SeqInfo['alpha'], 
		                                           SeqInfo['numSeq'], SeqInfo['seqLen'], signalFlag);	
	return NegSeqDict;

def CreateNoSignalFastaFile(configFile):
	confMap = SeqGenUtils.GetConf(configFile)
	NegSeqDict = CreateNoSignalDict(confMap, False);
	OutFileName = SeqGenUtils.GetNoSignalOutFileName(confMap);


if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	CreateNoSignalFastaFile(configFile)
