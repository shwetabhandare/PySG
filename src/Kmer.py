import sys
import NoSignal
import yaml
import random
import SeqGenUtils

def AddSignalFromStart(SeqDict, kmerToEmbed, numSeqsWithSignal, locationFromStart):
	keysToReplace = random.sample(SeqDict, numSeqsWithSignal)
	for key in keysToReplace:
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

	NegativeFileName = confMap["sequence"]["fastaFile"]
	NumSeqsToGenerate = int(confMap["sequence"]["numSeq"])
	SeqLength = int(confMap["sequence"]["seqLen"])
	OutFileName = confMap["sequence"]["outFastaFile"]
	kmerToEmbed = confMap["sequence"]["kmer"]
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

	SeqDict = AddSignalFromStart(SeqDict, kmerToEmbed, numSeqsWithSignal, locationFromStart)
	SeqGenUtils.WriteSeqDictToFile(SeqDict, OutFileName);
