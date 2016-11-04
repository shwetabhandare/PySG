import sys
import parseDreme
import parseRealKmers
import SeqGenUtils
import compareKmerCommon
import compareKmers

def isBinder(kmer, kmerScore, cutOff):
	if kmerScore < cutOff:
		#print "FP: ", kmer, ",", kmerScore
		return False;
	return True;

def isNonBinder(kmer, kmerScore, cutOff):
	if kmerScore > cutOff:
		return False;

	return True;

def getNumbersForNegativeSet(outStr):
	tupleList = outStr.split("|")
	numTP = numFP = numFN = 0;

	for tuples in tupleList:
		kmer, kmerScore = tuples.split(",")
		kmerScore = float(kmerScore)

	return numTP, numFP, numFN

def getNumbersForMatchedKmers(outStr, thresholdPercent, positive=True):
	tupleList = outStr.split("|")
	totalTP = totalFP = totalFN = totalTN = 0;

	for tuples in tupleList:
		kmer, kmerScore, maxMotifScore = tuples.split(",")
		threshold = thresholdPercent * float(maxMotifScore);
		kmerScore = float(kmerScore)
		if positive:
			if isBinder(kmer, kmerScore, threshold):
				totalTP = totalTP + 1;
				break;
			else:
				continue;
		else:
			if isNonBinder(kmer, kmerScore, threshold):
				totalTN = 1;
				continue;
			else:
				totalFP = totalFP + 1;
				break;

	# None of the predicted motifs scored higher than 0.6 * max score.
	if positive and totalTP == 0:
		totalFN = 1;

	#We found more than one motif that didn't meet the cut-off.
	if not positive and totalFP == 0:
		totalTN = 1;

	return totalTP, totalFP, totalFN, totalTN;

def computeNumbersForSequence(dremeFile, seqFile, thresholdPercent, positive=True):
	seqDict  = SeqGenUtils.fasta_read(seqFile);
	totalExamples = len(seqDict)
	pssmList = parseDreme.getPSSMListFromDremeFile(dremeFile)
	totalTP = totalFP = totalFN = totalTN = 0;

	for seqid, seq in seqDict.iteritems():
		outStr = compareKmerCommon.getKmerAndScoreFromPSSM(pssmList, seq)
		numTP, numFP, numFN, numTN = getNumbersForMatchedKmers(outStr, thresholdPercent, positive)
		totalTP = totalTP + numTP
		totalFN = totalFN + numFN
		totalFP = totalFP + numFP
		totalTN = totalTN + numTN

	print "Total TP: ", totalTP, ", Total FP: ", totalFP, ", Total FN: ", totalFN, ", Total TN:", totalTN
	return totalExamples, totalTP, totalFP, totalFN, totalTN;

if __name__ == "__main__":
	import sys

	dremeFile = sys.argv[1]
	posSeqFile = sys.argv[2]
	negSeqFile = sys.argv[3]
	thresholdPercent = float(sys.argv[4])

	totalPos, numPosTP, numPosFP, numPosFN, numPosTN = computeNumbersForSequence(dremeFile, posSeqFile, thresholdPercent)
	totalNeg, numNegTP, numNegFP, numNegFN, numNegTN = computeNumbersForSequence(dremeFile, negSeqFile, thresholdPercent, False)
	sensitivity, ppv = compareKmers.GetSensitivityAndPPV((numPosTP + numNegTP) , (numPosFP + numNegFP), (numPosFN + numNegFN))
	accuracy = compareKmers.GetAccuracy( (numPosTP + numNegTP), (numPosTN + numNegTN),  (totalPos + totalNeg) )

	print "Senitivity: ", sensitivity, ", PPV: ", ppv, ", Accuracy: ", accuracy;



