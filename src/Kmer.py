import sys
import NoSignal
import yaml
import random
import SeqGenUtils
import TAMO_Motif

def AddSignalFromStart(SeqDict, type, confMap, numSeqsWithSignal, locationFromStart):
	generateKmer = False;
	if type == 'pwm':
		motiffile = confMap["sequence"]["pwmFile"]
		motif = TAMO_Motif.Make_PWM_Motif(motiffile)
		generateKmer = True;
	elif type == 'motif':
		motiffile = confMap["sequence"]["pwmFile"]
		motif = TAMO_Motif.Make_PWM_Motif(motiffile)
		generateKmer = True;
	else:
		kmerToEmbed = confMap["sequence"]["kmer"]

	keysToReplace = random.sample(SeqDict, numSeqsWithSignal)
	for key in keysToReplace:
		if generateKmer:
			kmerToEmbed = motif.random_kmer();
			print "k-mer to embed ", kmerToEmbed;

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
	
if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)
	kmerToEmbed = ""
	type = "None";

	NegativeFileName = confMap["sequence"]["fastaFile"]
	NumSeqsToGenerate = int(confMap["sequence"]["numSeq"])
	SeqLength = int(confMap["sequence"]["seqLen"])
	OutFileName = confMap["sequence"]["outFastaFile"]

	if confMap['sequence'].get('kmer'):
		kmerToEmbed = confMap["sequence"]["kmer"]
		type = "kmer"
	elif confMap['sequence'].get('pwmFile'):
		motifFile = confMap['sequence'].get('pwmFile');
		type = "pwm"
	elif confMap['sequence'].get('textMotif'):
		textMotif = confMap['sequence'].get('kmer');
		type = "motif"
	
	if confMap['sequence'].get('locationFromEnd'):
		locationFromEnd = int(confMap["sequence"]["locationFromEnd"])
		locationFromStart = SeqLength - locationFromEnd;
	elif confMap['sequence'].get('locationFromStart'):
		locationFromStart = int(confMap["sequence"]["locationFromStart"])
	else:
		locationFromStart = random.randint(0, SeqLength - len(kmerToEmbed))

	print "Location to embed signal : " , str(locationFromStart)
	numSeqsWithSignal = int(confMap["sequence"]["seqWithSignal"])

	SeqDict = NoSignal.GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, 
	          SeqLength, OutFileName);

	SeqDict = AddSignalFromStart(SeqDict, type, confMap, numSeqsWithSignal, locationFromStart)
	SeqGenUtils.WriteSeqDictToFile(SeqDict, OutFileName);
