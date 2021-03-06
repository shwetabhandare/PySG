from random import choice
import yaml
import os, fnmatch
from numpy import *
import numpy as np
import fasta

np.set_printoptions(precision=2)

def weightedchoice(items): # this doesn't require the numbers to add up to 100
	return choice("".join(x * y for x, y in items))

def createFastaFileFromKmers(kmerList, faPrefix):
   counter = 0;
   new_lines = "";

   for kmer in kmerList:
      add_line = ">" + faPrefix + str(counter) + "\n"
      new_lines = new_lines + add_line
      new_lines = new_lines + kmer  + "\n";
      counter = counter + 1

   return new_lines;

def writeFastaLinesToFile(kmerFastaLines, outFile):
	f1 = open(outFile, "w")
	for line in kmerFastaLines:
		f1.write(line)
	f1.close();
def ChangeUsToTs(seq_dict):
	for header, sequence in seq_dict.iteritems():
		sequence = sequence.upper().replace('U', 'T');
		seq_dict[header] = sequence;
	return seq_dict;
	
def fasta_read(file_name):
   """read the sequence from a file in fasta format"""
   seq_dict = dict()
   for record in fasta.fasta_itr(file_name):
      seq_dict[record.header] = record.sequence;
   return seq_dict;

def GetConf(configFile):
	f = open(configFile);
	confMap = yaml.load(f)
	f.close();
	return confMap;

def GetNoSignalOutFileName(confMap):
	return confMap["sequence"]["nosignal"]["outNoSignalFastaFile"]

def GetSequenceInfoFromConfMap(confMap):

	SeqInfoDict = dict();
	NegativeFileName = ""
	alpha = 0;

	NumNoSignalSeq = int(confMap["sequence"]["nosignal"]["numSeq"])
	if confMap["sequence"]["nosignal"].get("fastaFile"):
		NegativeFileName = confMap["sequence"]["nosignal"]["fastaFile"]	
	else:
		alpha = int(confMap["sequence"]["nosignal"]["alpha"])
	if confMap["sequence"]["nosignal"].get("seqLen"):
		SeqLength = int(confMap["sequence"]["nosignal"]["seqLen"])
	else:
		SeqLength = 200;

	SeqInfoDict['seqLen'] = SeqLength;
	SeqInfoDict['inputName'] = NegativeFileName;
	SeqInfoDict['numSeq'] = NumNoSignalSeq	;
	SeqInfoDict['alpha'] = alpha;

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

def GetDirichletDistribution(seqBackGroundDict, scaleFactor, NumSeqsToGenerate):
	utr_dist = []
	for key, value in seqBackGroundDict.iteritems():
		utr_dist.append(value)
	alpha = scaleFactor * np.array(utr_dist)
	s = np.random.dirichlet(alpha, NumSeqsToGenerate);
	return s;

def GenerateNoSignalFromDirichlet(seqDistDirichletList, seqBackGroundDict,
	SeqLength, signalFlag):
	seq = ""
	keys = []
	for key, value in seqBackGroundDict.iteritems():
		keys.append(key)
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
		if signalFlag == False:
			key = "NoSignal_" + str(seqNum);
		else:
			key = "Signal_" + str(seqNum);
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
