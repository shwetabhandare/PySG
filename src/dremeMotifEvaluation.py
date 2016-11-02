import sys
import parseDreme
import parseRealKmers
import SeqGenUtils
import compareKmerCommon

def getNumbersForKmerAndScore(kmer, kmerScore, cutOff, positive=True):
	numTP = numFP = numFN = 0;
	if positive:
		if kmerScore < 0.0:
			print "FP: ", kmer, ",", kmerScore
			numFP = numFP + 1
		elif kmerScore < cutOff:
			print "FN: ", kmer, ",", kmerScore
			numFN = numFN + 1
		else:
			print "TP: ", kmer, ",", kmerScore
			numTP = numTP + 1
	else:
		if kmerScore < 0.0:
			numTP = numTP + 1
		elif kmerScore > cutOff:
			numFP = numFP + 1
		else:
			numFN = numFN + 1
	return  numTP, numFP, numFN;

def getNumbersForNegativeSet(outStr):
	tupleList = outStr.split("|")
	numTP = numFP = numFN = 0;

	for tuples in tupleList:
		kmer, kmerScore = tuples.split(",")
		kmerScore = float(kmerScore)

	return numTP, numFP, numFN
def getNumbersForMatchedKmers(outStr):
	tupleList = outStr.split("|")
	totalTP = totalFP = totalFN = 0;


	for tuples in tupleList:
		kmer, kmerScore = tuples.split(",")
		kmerScore = float(kmerScore)
		numTP, numFP, numFN = getNumbersForKmerAndScore(kmer, kmerScore, 3.0)
		totalTP = totalTP + numTP;
		totalFP = totalFP + numFP;
		totalFN = totalFN + numFN;

	return totalTP, totalFP, totalFN;

def getHighScoringKmer(dremeFile, seqFile):
	seqDict  = SeqGenUtils.fasta_read(seqFile);
	pssmList = parseDreme.getPSSMListFromDremeFile(dremeFile)
	totalTP = totalFP = totalFN = 0;

	for seqid, seq in seqDict.iteritems():
		outStr = compareKmerCommon.getKmerAndScoreFromPSSM(pssmList, seq)
		numTP, numFP, numFN = getNumbersForMatchedKmers(outStr)
		#numTP, numFP, numFN = getNumbersForNegativeSet(outStr)
		totalTP = totalTP + numTP
		totalFN = totalFN + numFN
		totalFP = totalFP + numFP

	print "Total TP: ", totalTP, ", Total FP: ", totalFP, ", Total FN: ", totalFN
	return totalTP, totalFP, totalFN;

if __name__ == "__main__":
	import sys
	dremeFile = sys.argv[1]
	seqFile = sys.argv[2]
	getHighScoringKmer(dremeFile, seqFile)

