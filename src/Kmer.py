import sys
import NoSignal
import yaml
import random
import SeqGenUtils
import TAMO_Motif

def EmbedMotif(SeqDict, type, confMap, numSeqsWithSignal, locationFromStart):
	generateKmer = False;
	motifBackGround = ""

	if type == "pwm" or type == "motif":
		#Check if background is provided.
		if confMap['sequence'].get('motifBackGround'):
			motifBackGround = confMap["sequence"]["motifBackGround"]
			print motifBackGround;

	if type == 'pwm':
		motiffile = confMap["sequence"]["pwmFile"]
		motif = TAMO_Motif.Make_PWM_Motif(motiffile, motifBackGround)
		generateKmer = True;
	elif type == 'motif':
		textMotif  = confMap["sequence"]["textMotif"]
		motif = TAMO_Motif.Make_Text_Motif(textMotif)
		generateKmer = True;
	else:
		kmerToEmbed = confMap["sequence"]["kmer"]


	keysToReplace = random.sample(SeqDict, numSeqsWithSignal)
	for key in keysToReplace:
		if generateKmer:
			kmerToEmbed = motif.random_kmer();

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

def GetMotifLocation(confMap):
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
	CreateNoSignalFastaFile(configFile)

	if confMap['sequence'].get('motif'):
		motifType = GetMotifType(confMap);
		locationFromStart = GetMotifLocation(confMap);
		numSeqsWithSignal = int(confMap["sequence"]["motif"]["seqWithSignal"])
	else:
		print "Unable to generate sequences with motif"

	# Embed motif into the sequences.
	SeqDict = EmbedMotif(SeqDict, motifType, confMap, numSeqsWithSignal, locationFromStart)

	# Write sequences to file.
	SeqGenUtils.WriteSeqDictToFile(SeqDict, OutFileName);
