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
	
if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	confMap = SeqGenUtils.GetConf(configFile)
	kmerToEmbed = ""
	type = "None";
	NegativeFileName = "";
	NegSeqDict  = dict();	

	if confMap["sequence"].get("fastaFile"):
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
		textMotif = confMap['sequence'].get('textMotif');
		type = "motif"
	
	if confMap['sequence'].get('locationFromEnd'):
		locationFromEnd = int(confMap["sequence"]["locationFromEnd"])
		locationFromStart = SeqLength - locationFromEnd;
	elif confMap['sequence'].get('locationFromStart'):
		locationFromStart = int(confMap["sequence"]["locationFromStart"])
	else:
		locationFromStart = random.randint(0, SeqLength - len(kmerToEmbed))

	numSeqsWithSignal = int(confMap["sequence"]["seqWithSignal"])

	if NegativeFileName != "":

		# Generate sequences that contain no signal.
		SeqDict = NoSignal.GenerateNoSignalSequences(NegativeFileName, NumSeqsToGenerate, 
	   	       SeqLength, OutFileName);
	else:
		SeqDict = NoSignal.GenerateNoSignalWithDirichlet(confMap, NumSeqsToGenerate, SeqLength);

	# Embed motif into the sequences.
	SeqDict = EmbedMotif(SeqDict, type, confMap, numSeqsWithSignal, locationFromStart)

	# Write sequences to file.
	SeqGenUtils.WriteSeqDictToFile(SeqDict, OutFileName);
