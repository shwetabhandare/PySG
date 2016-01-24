import sys
import NoSignal
import yaml
import random
import SeqGenUtils

def EmbedKmer(SeqDict, kmerToEmbed):
	for key, value in SeqDict.iteritems():
		location = random.randint(0, len(value) - 1 - len(kmerToEmbed));

	
if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)

	NegativeFileName = confMap["sequence"]["fastaFile"]
	NumSeqsToGenerate = int(confMap["sequence"]["numSeq"])
	SeqLength = int(confMap["sequence"]["seqLen"])
	OutFileName = confMap["sequence"]["outFastaFile"]
	kmerToEmbed = confMap["sequence"]["kmer"]

	SeqDict = NoSignal.GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, SeqLength, OutFileName);

	SeqDict = EmbedKmer(SeqDict, kmerToEmbed)
