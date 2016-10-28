import TAMO_Motif
import csv
import SeqGenUtils
import parseKspectrum
import os;



def GetTopXKmers(topKmerFile, N):
	topKmerList = [];
	with open(topKmerFile, 'rb') as csvfile:
		reader = csv.reader(csvfile);
		kmersFound = 0;
		for row in reader:
			if kmersFound < N:
				topKmerList.append(row[1])
				kmersFound = kmersFound + 1;
	#print "Top-",N, " kmer list: ", topKmerList;
	return topKmerList;

def GetMotifForSignal(signalType, signalFile):

	if signalType == "PFM":
		motif = TAMO_Motif.Make_PFM_Motif(signalFile);
	else: #PWM
		motif = TAMO_Motif.Make_PWM_Motif(signalFile);
	return motif;

def matchKmerListWithBestKmers(sequence, kmerList, bestScoreKmerList):
	bestKmerList = [i[1] for i in bestScoreKmerList];
	numKmersInBestList = 0;
	matchedKmerList = list()

	for predictedKmer in kmerList:
		if predictedKmer in bestKmerList:
			#print "Found ", predictedKmer, " in bestKmerList"
			numKmersInBestList = numKmersInBestList + 1;
			matchedKmerList.append(predictedKmer)

	return numKmersInBestList, matchedKmerList;

def FindBestScoreSeqForMotif(motif, seqDict, kmerList):
	totalMatchKmers = 0;
	totalKmersSet = set();
	for header, sequence in seqDict.iteritems():
		bestScoreKmerList = motif.bestseqs(0.90);
		#print "Finding predicted kmers for sequence: ", header;
		matchKmers, matchedKmerList = matchKmerListWithBestKmers(sequence, kmerList, bestScoreKmerList)
		totalMatchKmers = totalMatchKmers + matchKmers;
		tempSet = set(matchedKmerList)
		totalKmersSet.update(tempSet)
		#print "Found ", matchKmers, " for sequence: ", header;

		#print type(bestScoreKmers), len(bestScoreKmers)
		#print bestScoreKmers;
	print "Total Match Kmers : ", totalMatchKmers;
	print "Predicted kmers matched: ", totalKmersSet;

def addFastaHeaderToKmers(kmerDict, faPrefix):
   counter = 0;
   new_lines = "";

   for kmer, score in kmerDict.iteritems():
      add_line = ">" + faPrefix + str(counter) + "\n"
      new_lines = new_lines + add_line
      new_lines = new_lines + kmer  + "\n";
      counter = counter + 1

   return new_lines;


if __name__ == "__main__":
	import sys

	featureFile = sys.argv[1]
	sequenceFile = sys.argv[2]
	numKmers = int(sys.argv[3])

	kmerDict = parseKspectrum.FindKspectrumKmers(featureFile, numKmers)
	HuR_ReString = '[^-](\d+\.\d+)\,HuR_([ATGC]+)'
	TTP_ReString = '[^-](\d+\.\d+)\,TTP_([ATGC]+)'

	HuRKmerDict = parseKspectrum.FindRBPSpecificKmers(featureFile, HuR_ReString, numKmers);
	TTPKmerDict = parseKspectrum.FindRBPSpecificKmers(featureFile, TTP_ReString, numKmers);


	kmerFastaLines = addFastaHeaderToKmers(kmerDict, "FeatureKmer")
	HuRFastaLines = addFastaHeaderToKmers(HuRKmerDict, "HuR_Specific")
	TTPFastaLines = addFastaHeaderToKmers(TTPKmerDict, "TTP_Specific")

	filename_prefix = os.path.splitext(sequenceFile)[0];
	filename_ext = os.path.splitext(sequenceFile)[1];
	HuR_Filename = filename_prefix + "_HuR" + filename_ext;
	TTP_Filename = filename_prefix + "_TTP" + filename_ext;

	SeqGenUtils.writeFastaLinesToFile(kmerFastaLines, sequenceFile)
	SeqGenUtils.writeFastaLinesToFile(HuRFastaLines, HuR_Filename)
	SeqGenUtils.writeFastaLinesToFile(TTPFastaLines, TTP_Filename)
	# topKmerFile = sys.argv[1]
	# signalType = sys.argv[2] ## PWM or PFM
	# signalFile = sys.argv[3]
	# seqFile = sys.argv[4]
	# topKmers = int(sys.argv[5])

	# kmerList = GetTopXKmers(topKmerFile, topKmers);
	# motif = GetMotifForSignal(signalType, signalFile);
	# seqDict  = SeqGenUtils.fasta_read(seqFile);
	# FindBestScoreSeqForMotif(motif, seqDict, kmerList)
