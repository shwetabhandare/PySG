import TAMO_Motif
import csv
import SeqGenUtils



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

if __name__ == "__main__":
	import sys
	topKmerFile = sys.argv[1]
	signalType = sys.argv[2] ## PWM or PFM
	signalFile = sys.argv[3]
	seqFile = sys.argv[4]
	topKmers = int(sys.argv[5])

	kmerList = GetTopXKmers(topKmerFile, topKmers);
	motif = GetMotifForSignal(signalType, signalFile);
	seqDict  = SeqGenUtils.fasta_read(seqFile);
	FindBestScoreSeqForMotif(motif, seqDict, kmerList)