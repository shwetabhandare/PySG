import sys
import NoSignal
import yaml
import random
import SeqGenUtils
import TAMO_Motif

def GetKmersToEmbed(type, numSeqsWithSignal, confMap):
	kmers = list();
	generateKmer = False;
	motifBackGround = ""

	if type == "pwm" or type == "motif":
		if confMap['sequence']['signal'].get('motifBackGround'):
			motifBackGround = confMap["sequence"]["signal"]["motifBackGround"]

	if type == 'pwm':
		motiffile = confMap["sequence"]['signal']["pwmFile"]
		motif = TAMO_Motif.Make_PWM_Motif(motiffile, motifBackGround)
		generateKmer = True;
	elif type == 'motif':
		textMotif  = confMap["sequence"]['signal']["textMotif"]
		motif = TAMO_Motif.Make_Text_Motif(textMotif)
		generateKmer = True;
	else:
		kmerToEmbed = confMap["sequence"]['signal']["kmer"]

	for i in range(numSeqsWithSignal):
		if generateKmer:
			kmerToEmbed = motif.random_kmer();

		print "Appending kmer: ", kmerToEmbed;
		kmers.append(kmerToEmbed);
	return kmers;

def EmbedMotif(SeqDict, kmerList, numSeqsWithSignal, locationFromStart):
	generateKmer = False;
	motifBackGround = ""
	print "Num Sequences with Signal", str(numSeqsWithSignal)
	keysToReplace = random.sample(SeqDict, numSeqsWithSignal)
	print len(keysToReplace)
	for idx, key in enumerate(keysToReplace):
		kmerToEmbed = kmerList[idx];

		value = SeqDict[key]
		endIndex = locationFromStart + len(kmerToEmbed)
		if (endIndex > len(value)):
			shiftStart = endIndex - len(value);
			locationFromStart = locationFromStart - shiftStart
		newValue = value[:locationFromStart] + kmerToEmbed + value[endIndex:]
		SeqDict[key] = newValue;
	return SeqDict;


def PrintConfMap(confMap):
	for doc in confMap:
		for k,v in doc.items(): 
			print k, "->", v 
			print "\n",

def GetMotifType(confMap):
	kmerToEmbed = ""
	if confMap['sequence']['signal'].get('kmer'):
		kmerToEmbed = confMap["sequence"]['signal']["kmer"]
		motifType = "kmer"
	elif confMap['sequence']['signal'].get('pwmFile'):
		motifFile = confMap['sequence']['signal'].get('pwmFile');
		motifType = "pwm"
	elif confMap['sequence']['signal'].get('textMotif'):
		textMotif = confMap['sequence']['signal'].get('textMotif');
		motifType = "motif"
	return motifType;

def GetMotifLocation(confMap, SeqLength):
	if confMap['sequence']['signal'].get('locationFromEnd'):
		locationFromEnd = int(confMap["sequence"]["motif"]["locationFromEnd"])
		locationFromStart = SeqLength - locationFromEnd;
	elif confMap['sequence']['signal'].get('locationFromStart'):
		locationFromStart = int(confMap["sequence"]['motif']["locationFromStart"])
	else:
		locationFromStart = random.randint(0, SeqLength)
	return locationFromStart;

def GetSignalSeqInfo(confMap, SeqLength):
	SignalSeqInfo = dict();

	if confMap['sequence'].get('signal'):
		motifType = GetMotifType(confMap);
		locationFromStart = GetMotifLocation(confMap, SeqLength);
		numSeqsWithSignal = int(confMap["sequence"]["signal"]["seqWithSignal"])
		SignalSeqInfo['motifType'] = motifType;
		SignalSeqInfo['locationFromStart'] = locationFromStart;
		SignalSeqInfo['seqWithSignal'] = numSeqsWithSignal;
	return SignalSeqInfo;


def CreateFastaWithSignal(configFile):
	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)
	SeqLength = confMap["sequence"]["nosignal"]["seqLen"]
	SignalSeqInfo = GetSignalSeqInfo(confMap, SeqLength);
	if len(SignalSeqInfo) == 0:
		print "Did not find any signal data in the yaml file: ", configFile;
	else:
		kmerList = GetKmersToEmbed(SignalSeqInfo['motifType'], SignalSeqInfo['seqWithSignal'], confMap);
		SeqDict = NoSignal.CreateNoSignalDict(confMap);

		# Embed motif into the sequences.
		SeqDict = EmbedMotif(SeqDict, kmerList, SignalSeqInfo['seqWithSignal'], SignalSeqInfo['locationFromStart']);

		NegSeqDict = NoSignal.CreateNoSignalDict(confMap);

		NoSignalOutFileName = confMap["sequence"]["nosignal"]["outNoSignalFastaFile"]
		MotifOutFileName = confMap["sequence"]["signal"]["outSignalFile"]
		# Write sequences to file.
		SeqGenUtils.WriteSeqDictToFile(SeqDict, MotifOutFileName);
		SeqGenUtils.WriteSeqDictToFile(NegSeqDict, NoSignalOutFileName);

if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	CreateFastaWithSignal(configFile);
