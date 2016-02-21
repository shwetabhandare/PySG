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
		if confMap['sequence']['motif'].get('motifBackGround'):
			motifBackGround = confMap["sequence"]["motif"]["motifBackGround"]

	if type == 'pwm':
		motiffile = confMap["sequence"]['motif']["pwmFile"]
		motif = TAMO_Motif.Make_PWM_Motif(motiffile, motifBackGround)
		generateKmer = True;
	elif type == 'motif':
		textMotif  = confMap["sequence"]['motif']["textMotif"]
		motif = TAMO_Motif.Make_Text_Motif(textMotif)
		generateKmer = True;
	else:
		kmerToEmbed = confMap["sequence"]['motif']["kmer"]

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
		print newValue
		SeqDict[key] = newValue;
	return SeqDict;


def PrintConfMap(confMap):
	for doc in confMap:
		for k,v in doc.items(): 
			print k, "->", v 
			print "\n",

def GetMotifType(confMap):
	if confMap['sequence']['motif'].get('kmer'):
		kmerToEmbed = confMap["sequence"]['motif']["kmer"]
		motifType = "kmer"
	elif confMap['sequence']['motif'].get('pwmFile'):
		motifFile = confMap['sequence']['motif'].get('pwmFile');
		motifType = "pwm"
	elif confMap['sequence']['motif'].get('textMotif'):
		textMotif = confMap['sequence']['motif'].get('textMotif');
		motifType = "motif"
	return motifType;

def GetMotifLocation(confMap, SeqLength):
	if confMap['sequence']['motif'].get('locationFromEnd'):
		locationFromEnd = int(confMap["sequence"]["motif"]["locationFromEnd"])
		locationFromStart = SeqLength - locationFromEnd;
	elif confMap['sequence']['motif'].get('locationFromStart'):
		locationFromStart = int(confMap["sequence"]['motif']["locationFromStart"])
	else:
		locationFromStart = random.randint(0, SeqLength - len(kmerToEmbed))
	return locationFromStart;

if __name__ == "__main__":
	import sys

	kmerToEmbed = ""
	type = "None";

	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)
	SeqLength = confMap["sequence"]["nosignal"]["seqLen"]

	if confMap['sequence'].get('motif'):
		motifType = GetMotifType(confMap);
		locationFromStart = GetMotifLocation(confMap, SeqLength);
		numSeqsWithSignal = int(confMap["sequence"]["motif"]["seqWithSignal"])
	else:
		print "Unable to generate sequences with motif"

	kmerList = GetKmersToEmbed(motifType, numSeqsWithSignal, confMap);


	SeqDict = NoSignal.CreateNoSignalDict(confMap);

	# Embed motif into the sequences.
	SeqDict = EmbedMotif(SeqDict, kmerList, numSeqsWithSignal, locationFromStart);

	NegSeqDict = NoSignal.CreateNoSignalDict(confMap);

	NoSignalOutFileName = confMap["sequence"]["nosignal"]["outNoSignalFastaFile"]
	MotifOutFileName = confMap["sequence"]["motif"]["outSignalFile"]


	# Write sequences to file.
	SeqGenUtils.WriteSeqDictToFile(SeqDict, MotifOutFileName);
	SeqGenUtils.WriteSeqDictToFile(NegSeqDict, NoSignalOutFileName);