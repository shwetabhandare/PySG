import sys
import parseDreme
import parseRealKmers
import SeqGenUtils
import compareKmerCommon
import compareKmers


def getPredictedKmerList(pssmList, seq):
	kmerREString = compareKmerCommon.getKmerFromPSSM(pssmList, seq);
	kmerREString = kmerREString[1:-1]
	predictedKmerList = kmerREString.split("|")

	return predictedKmerList;

def compareRealAndPredicted(realKmerDict, seqDict, pssmList, positive=True):
	totalFN = totalTP = totalFP = totalTN = 0;
	kmerEmbedded = True;

	for seqid, seq in seqDict.iteritems():
		realKmer, realStart, realEnd = compareKmerCommon.getRealKmerDetails(realKmerDict, seqid);
		if realKmer == "" and realStart == 0 and realEnd == 0:
			#print "Sequence ID: ", seqid , " should be treated like a negative."
			kmerEmbedded = False;

		#print "===== Seq ID: ", seqid, ": EMBEDDED K-MER: ", realKmer, ", PREDICTED: ", kmerREString
		predictedKmerList = getPredictedKmerList(pssmList, seq);
		kmerFoundInSeq = compareKmerCommon.isRealKmerInPredictedKmerList(realKmer, predictedKmerList)
		numTP, numFP, numFN, numTN = compareKmerCommon.getNumbers(kmerFoundInSeq, kmerEmbedded, positive);
		totalTP = totalTP + numTP;
		totalFP = totalFP + numFP;
		totalFN = totalFN + numFN;
		totalTN = totalTN + numTN;


	return totalTP, totalFP, totalFN, numTN;


def computeSequenceBasedDREMEResults(dremeFile, realCsvFile, posSeqFile, negSeqFile):
	totalPosTP = totalPosFP = totalPosFN = totalPosTN = 0;
	totalNegTP = totalNegFP = totalNegFN = totalNegTN = 0;

	posSeqDict  = SeqGenUtils.fasta_read(posSeqFile);
	negSeqDict  = SeqGenUtils.fasta_read(negSeqFile);
	realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);
	pssmList = parseDreme.getPSSMListFromDremeFile(dremeFile)

	numPosTP, numPosFP, numPosFN, numPosTN = compareRealAndPredicted(realKmerDict, posSeqDict, pssmList, positive=True)
	numNegTP, numNegFP, numNegFN, numNegTN = compareRealAndPredicted(realKmerDict, negSeqDict, pssmList, False)

	print "Positive: TP: ", numPosTP, ", FP: ", numPosFP, ", FN: ", numPosFN, ", TN: ", numPosTN
	print "Negative: TP: ", numNegTP, ", FP: ", numNegFP, ", FN: ", numNegFN, ", TN: ", numNegTN

	totalPos = len(posSeqDict)
	totalNeg = len(negSeqDict)

	sensitivity, ppv = compareKmers.GetSensitivityAndPPV((numPosTP + numNegTP) , (numPosFP + numNegFP), (numPosFN + numNegFN))
	accuracy = compareKmers.GetAccuracy( (numPosTP + numNegTP), (numPosTN + numNegTN),  (totalPos + totalNeg) )
	specificity = compareKmers.GetSpecificity( (numPosFP + numNegFP), totalNeg);

	print "Senitivity: ", sensitivity, ", PPV: ", ppv, ", Accuracy: ", accuracy, ", Specificity: ", specificity;
	return sensitivity, ppv;	

if __name__ == "__main__":
	import sys
	posSeqFile = sys.argv[1]
	negSeqFile = sys.argv[2]
	dremeFile = sys.argv[3]
	realCsvFile = sys.argv[4]

	computeSequenceBasedDREMEResults(dremeFile, realCsvFile, posSeqFile, negSeqFile)

