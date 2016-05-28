import sys
import NoSignal
import yaml
import random
import SeqGenUtils
import TAMO_Motif
import os

def GetSequencesWithSignal(SignalSeqInfo):
	signalPercent = SignalSeqInfo['seqWithSignal']
	locationFromStart = SignalSeqInfo['locationFromStart'];
	totalSeq = SignalSeqInfo['numSeq'];
	numSeqsWithSignal = (signalPercent * totalSeq)/100

	return numSeqsWithSignal;

def CreateLengthMatchingKmers(kmers, numSeqsWithSignal):

	while(len(kmers) < numSeqsWithSignal):
		randomIndex = random.randint(0, len(kmers) - 1)
		kmers.append(kmers[randomIndex])
		#print "Appending Kmer: ", kmers[randomIndex], " from index: ", randomIndex;
	return kmers;

def GetKmersFromStructureFile(structureFile, numSeqsWithSignal):
	structureDict = SeqGenUtils.fasta_read(structureFile);
	structureDict = SeqGenUtils.ChangeUsToTs(structureDict);
	originalKmers = structureDict.values();
	
	updatedKmers = CreateLengthMatchingKmers(originalKmers, numSeqsWithSignal)
	return updatedKmers;




def GetKmersToEmbed(SignalSeqInfo, confMap):
	kmers = list();
	generateKmer = False;
	motifBackGround = ""

	type = SignalSeqInfo['motifType'];
	#print "Motif Type: ", type;

	if type == "pwm" or type == "motif":
		if confMap['sequence']['signal'].get('seqBackGround'):
			motifBackGround = confMap["sequence"]["signal"]["seqBackGround"]

	if type == 'pwm':
		motiffile = confMap["sequence"]['signal']["pwmFile"]
		motif = TAMO_Motif.Make_PWM_Motif(motiffile, motifBackGround)
		generateKmer = True;
	elif type == 'pfm':
		motiffile = confMap["sequence"]['signal']["pfmFile"]
		motif = TAMO_Motif.Make_PFM_Motif(motiffile, motifBackGround)
		generateKmer = True;
	elif type == 'motif':
		textMotif  = confMap["sequence"]['signal']["textMotif"]
		motif = TAMO_Motif.Make_Text_Motif(textMotif)
		generateKmer = True;
	elif type == "structure":
		structureFile = confMap["sequence"]['signal']['structureFile']
		numSeqsWithSignal = GetSequencesWithSignal(SignalSeqInfo);
		kmers = GetKmersFromStructureFile(structureFile, numSeqsWithSignal)
	else:
		kmerToEmbed = confMap["sequence"]['signal']["kmer"]

	numSeqsWithSignal = GetSequencesWithSignal(SignalSeqInfo)

	if type != "structure":
		for i in range(numSeqsWithSignal):
			if generateKmer:
				kmerToEmbed = motif.emit();

			#print "Appending kmer: ", kmerToEmbed;
			kmers.append(kmerToEmbed);
	return kmers;

def EmbedMotif(SeqDict, kmerList, SignalSeqInfo):
	generateKmer = False;
	motifBackGround = ""
	embeddedKmerDict = dict();

	numSeqsWithSignal = GetSequencesWithSignal(SignalSeqInfo);
	locationFromStart = SignalSeqInfo['locationFromStart'];

	#print "Num Sequences with Signal", str(numSeqsWithSignal)
	keysToReplace = random.sample(SeqDict, numSeqsWithSignal)
	#print len(keysToReplace)
	#print len(kmerList)

	for idx, key in enumerate(keysToReplace):
		kmerToEmbed = kmerList[idx];

		value = SeqDict[key]
		endIndex = locationFromStart + len(kmerToEmbed)
		if (endIndex > len(value)):
			shiftStart = endIndex - len(value);
			locationFromStart = locationFromStart - shiftStart
		newValue = value[:locationFromStart] + kmerToEmbed + value[endIndex:]
		SeqDict[key] = newValue;
		embeddedKmerDict[key] = [kmerToEmbed, locationFromStart]
	return SeqDict, embeddedKmerDict;


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
	elif confMap['sequence']['signal'].get('pfmFile'):
		motifFile = confMap['sequence']['signal'].get('pfmFile');
		motifType = "pfm"		
		#print "MOTIF FILE:", motifFile
	elif confMap['sequence']['signal'].get('textMotif'):
		textMotif = confMap['sequence']['signal'].get('textMotif');
		motifType = "motif"
	elif confMap['sequence']['signal'].get('structureFile'):
		textMotif = confMap['sequence']['signal'].get('structureFile');
		motifType = "structure"		
	return motifType;

def GetMotifLocation(confMap, SeqLength):
	if confMap['sequence']['signal'].get('locationFromEnd'):
		locationFromEnd = int(confMap["sequence"]["motif"]["locationFromEnd"])
		locationFromStart = SeqLength - locationFromEnd;
	elif confMap['sequence']['signal'].get('locationFromStart'):
		locationFromStart = int(confMap["sequence"]['signal']["locationFromStart"])
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

	SignalSeqInfo['numSeq'] = int(confMap['sequence']['nosignal']['numSeq'])


	return SignalSeqInfo;

def WriteKmersToFile(embeddedKmerDict, OutputFileName):
	OutputFile = open(OutputFileName , "w")
	for key, value in embeddedKmerDict.iteritems():
  		OutputFile.write(key);
  		OutputFile.write(",");
  		OutputFile.write(value[0]);
  		OutputFile.write(",");
  		OutputFile.write(str(value[1]));
		OutputFile.write("\n")
	OutputFile.close()


def CreateFastaWithSignal(configFile):
	confMap = SeqGenUtils.GetConf(configFile)
	SeqLength = confMap["sequence"]["nosignal"]["seqLen"]
	SignalSeqInfo = GetSignalSeqInfo(confMap, SeqLength);
	SeqDict = dict();
	
	if len(SignalSeqInfo) == 0:
		print "Did not find any signal data in the yaml file: ", configFile;
	else:
		# Get the list of k-mers to embed in sequences.
		kmerList = GetKmersToEmbed(SignalSeqInfo, confMap);

		# Get no signal sequences to embed the kmers into.
		SeqDict = NoSignal.CreateNoSignalDict(confMap, True);

		# Embed motif into the sequences.
		SeqDict, embeddedKmerDict = EmbedMotif(SeqDict, kmerList, SignalSeqInfo);


		MotifOutFileName = confMap["sequence"]["signal"]["outSignalFile"]

		KmersFileName = os.path.splitext(MotifOutFileName)[0] + ".kmers"
		WriteKmersToFile(embeddedKmerDict, KmersFileName);

		# Write sequences to file.
		SeqGenUtils.WriteSeqDictToFile(SeqDict, MotifOutFileName);
	return SeqDict;

if __name__ == "__main__":
	import sys
	configFile = sys.argv[1]
	CreateFastaWithSignal(configFile);

