from random import choice
import yaml
import os, fnmatch
from numpy import *
import numpy as np

def weightedchoice(items): # this doesn't require the numbers to add up to 100
	return choice("".join(x * y for x, y in items))


def GetConf(configFile):
	f = open(configFile);
	confMap = yaml.load(f)
	f.close();
	return confMap;

def GetNoSignalOutFileName(confMap):
	return confMap["sequence"]["nosignal"]["outFastaFile"]

def GetSequenceInfoFromConfMap(confMap):

	SeqInfoDict = dict();
	NegativeFileName = ""

	NumNoSignalSeq = int(confMap["sequence"]["nosignal"]["numSeq"])
	if confMap["sequence"]["nosignal"].get("fastaFile"):
		NegativeFileName = confMap["sequence"]["nosignal"]["fastaFile"]	
	if confMap["sequence"]["nosignal"].get("seqLen"):
		SeqLength = int(confMap["sequence"]["nosignal"]["seqLen"])
	else:
		SeqLength = 200;

	SeqInfoDict['seqLen'] = SeqLength;
	SeqInfoDict['inputName'] = NegativeFileName;
	SeqInfoDict['numSeq'] = NumNoSignalSeq;


	return SeqInfoDict;

def GetRandomNucleotide(items): 
    return choice("".join(x * y for x, y in items))

def WriteSeqDictToFile(NegSeqDict, OutputFileName):
	OutputFile = open(OutputFileName , "w")
	for key, value in NegSeqDict.iteritems():
		OutputFile.write(">")												#Write to output file
		OutputFile.write(key)
		OutputFile.write("\n")
		OutputFile.write(value)
		OutputFile.write("\n")
	OutputFile.close()

def GetDirichletDistribution(seqBackGroundDict, NumSeqsToGenerate):
	s = np.random.dirichlet((seqBackGroundDict["A"] * 100, seqBackGroundDict["C"] * 100,
		                      seqBackGroundDict["T"] * 100, seqBackGroundDict["G"] * 100), 
	                         NumSeqsToGenerate);
	return s;

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
		#print b
		count = 0;
		seq = ""
		for count in range(SeqLength):
			seq += GetRandomNucleotide(b)
		key = "NoSignal_" + str(seqNum);
		NegSeqDict[key] = seq;

	sorted(NegSeqDict, key=lambda key: NegSeqDict[key])

	#for key, value in NegSeqDict.iteritems():
	#	print key, len(value)
	return NegSeqDict;

def findFiles (path, filter):
   for root, dirs, files in os.walk(path):
      for file in fnmatch.filter(files, filter):
         yield os.path.join(root, file)
#seq = GetRandomSequence(50);
#print seq;
