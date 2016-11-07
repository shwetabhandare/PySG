import sys
import parseDreme
import parseRealKmers
import SeqGenUtils
import compareKmerCommon
import compareKmers


def isRealKmerInREString(realKmer, kmerREString):
	kmerREString = kmerREString[1:-1]
	predictedKmers = kmerREString.split("|")
	matching = [s for s in predictedKmers if realKmer in s]
	#print "Predicted K-mers: ", predictedKmers, ", Matching Kmer: ", matching

	return matching, len(predictedKmers)

def compareRealAndPredicted(realKmerDict, seqDict, pssmList, positive=True):
	numFN = numTP = numFP = numTN = 0;

	for seqid, seq in seqDict.iteritems():
		realKmer, realStart, realEnd = compareKmerCommon.getRealKmerDetails(realKmerDict, seqid);
		kmerREString = compareKmerCommon.getKmerFromPSSM(pssmList, seq);

		#print "===== Seq ID: ", seqid, ": EMBEDDED K-MER: ", realKmer, ", PREDICTED: ", kmerREString
		matchingKmers, numPredictedKmers = isRealKmerInREString(realKmer, kmerREString)
		numMatchingKmers = len(matchingKmers)

		if positive:
			if numMatchingKmers == 0: #no predicted k-mers matched real kmer.
				numFN = numFN + 1;
			else:
				numTP = numTP + 1; # should this be incremented by numMatchingKmers?
		else:
			if numMatchingKmers == 0: #no predicted k-mers matched real kmer.
				numTN = numTN + 1;
			else:
				numFP = numFP + 1; # should this be incremented by numMatchingKmers?
		#print "TP: ", numTP, ", FP: ", numFP, ", FN: ", numFN, ", numTN: ", numTN

	return numTP, numFP, numFN, numTN;


def computeSequenceBasedDREMEResults(dremeFile, realCsvFile, posSeqFile, negSeqFile):
	totalPosTP = totalPosFP = totalPosFN = totalPosTN = 0;
	totalNegTP = totalNegFP = totalNegFN = totalNegTN = 0;

	posSeqDict  = SeqGenUtils.fasta_read(posSeqFile);
	negSeqDict  = SeqGenUtils.fasta_read(negSeqFile);
	realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);
	pssmList = parseDreme.getPSSMListFromDremeFile(dremeFile)

	numPosTP, numPosFP, numPosFN, numPosTN = compareRealAndPredicted(realKmerDict, posSeqDict, pssmList, positive=True)
	numNegTP, numNegFP, numNegFN, numNegTN = compareRealAndPredicted(realKmerDict, negSeqDict, pssmList, False)

	#print "Positive: TP: ", numPosTP, ", FP: ", numPosFP, ", FN: ", numPosFN, ", TN: ", numPosTN
	#print "Negative: TP: ", numNegTP, ", FP: ", numNegFP, ", FN: ", numNegFN, ", TN: ", numNegTN

	totalPos = len(posSeqDict)
	totalNeg = len(negSeqDict)

	sensitivity, ppv = compareKmers.GetSensitivityAndPPV((numPosTP + numNegTP) , (numPosFP + numNegFP), (numPosFN + numNegFN))
	accuracy = compareKmers.GetAccuracy( (numPosTP + numNegTP), (numPosTN + numNegTN),  (totalPos + totalNeg) )
	specificity = compareKmers.GetSpecificity( (numPosFP + numNegFP), totalNeg);

	#print "Senitivity: ", sensitivity, ", PPV: ", ppv, ", Accuracy: ", accuracy, ", Specificity: ", specificity;
	return sensitivity, ppv;	

if __name__ == "__main__":
	import sys
	posSeqFile = sys.argv[1]
	negSeqFile = sys.argv[2]
	dremeFile = sys.argv[3]
	realCsvFile = sys.argv[4]

	computeSequenceBasedDREMEResults(predictedFile, realCsvFile, posSeqFile, negSeqFile)

