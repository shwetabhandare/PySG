import sys
import NoSignal
import yaml
import random
import SeqGenUtils

def AddSignal(SeqDict, kmerToEmbed, numSeqsWithSignal, locationFromStart):
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


	
if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)

	NegativeFileName = confMap["sequence"]["fastaFile"]
	NumSeqsToGenerate = int(confMap["sequence"]["numSeq"])
	SeqLength = int(confMap["sequence"]["seqLen"])
	OutFileName = confMap["sequence"]["outFastaFile"]
	kmerToEmbed = confMap["sequence"]["kmer"]
	numSeqsWithSignal = int(confMap["sequence"]["seqWithSignal"])
	locationFromStart = int(confMap["sequence"]["locationFromStart"])

	SeqDict = NoSignal.GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, 
	          SeqLength, OutFileName);

	SeqDict = AddSignal(SeqDict, kmerToEmbed, numSeqsWithSignal, locationFromStart)
	SeqGenUtils.WriteSeqDictToFile(SeqDict, OutFileName);
