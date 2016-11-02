import sys
import parseDreme
import parseRealKmers
import SeqGenUtils
import compareKmerCommon


def isRealKmerInREString(realKmer, kmerREString):
	kmerREString = kmerREString[1:-1]
	predictedKmers = kmerREString.split("|")
	matching = [s for s in predictedKmers if realKmer in s]
	print "Predicted K-mers: ", predictedKmers, ", Matching Kmer: ", matching

	return matching, len(predictedKmers)

def compareRealAndPredicted(realKmerDict, seqDict, pssmList):
	numFN = numTP = numFP = 0;

	for seqid, seq in seqDict.iteritems():
		realKmer, realStart, realEnd = compareKmerCommon.getRealKmerDetails(realKmerDict, seqid);
		kmerREString = compareKmerCommon.getKmerFromPSSM(pssmList, seq);

		print "===== Seq ID: ", seqid, ": EMBEDDED K-MER: ", realKmer, ", PREDICTED: ", kmerREString
		matchingKmers, numPredictedKmers = isRealKmerInREString(realKmer, kmerREString)
		numMatchingKmers = len(matchingKmers)

		if numMatchingKmers == 0: #no predicted k-mers matched real kmer.
			numFN = numFN + 1;
			numFP = numFP + numPredictedKmers;
		else:
			numTP = numTP + 1; # should this be incremented by numMatchingKmers?
			numFP = numFP + (numPredictedKmers - numMatchingKmers);
		print "TP: ", numTP, ", FP: ", numFP, ", FN: ", numFN
	return numTP, numFP, numFN;


if __name__ == "__main__":
	import sys
	seqFile = sys.argv[2]
	dremeFile = sys.argv[1]
	#realCsvFile = sys.argv[3]

	#realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);
	seqDict  = SeqGenUtils.fasta_read(seqFile);
	pssmList = parseDreme.getPSSMListFromDremeFile(dremeFile)
	for seqid, seq in seqDict.iteritems():
		bestKmer = compareKmerCommon.getKmerFromPSSM(pssmList, seq)

	#numTP, numFP, numFN = compareRealAndPredicted(realKmerDict, seqDict, pssmList)
