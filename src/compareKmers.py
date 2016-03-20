#import re, ahocorasick,random,time
import parseDreme
import parseRealKmers
import fasta
import SeqGenUtils
import re
import TAMO_Motif;
#from itertools import imap

def getPredictedKmerDict(dremeResultFile):
	return getPredictedDremeMotifs(dremeResultFile);

def getReadKmerDict(realKmersCsvFile):
	return GetRealKmerDict(realKmersCsvFile)


def getStartIndexAndUpdateNumbers(predictedStart, predictedEnd, realStart, realEnd, predKmer):
	numFP = 0;
	numFN = 0;

	if predictedStart > realStart:
		if predictedStart > realEnd:
			numFP = numFP + len(predKmer)
			numFN = numFN + (realEnd - realStart)
			#There is no k-mer to compare, so we set startIndex to realEnd
			# so that we never compare the k-mer. 
			startIndex = realEnd;
		else:
			startIndex = predictedStart;
			numFN = predictedStart - realStart;
	else:
		#Predicted k-mer starts before the real k-mer.
		if predictedEnd < realStart:
			#Predicted k-mer ends before the real k-mer starts,
			# There is no real k-mer to compare.
			startIndex = realStart;
			numFP = numFP + len(predKmer)
			# We can't update the numFN here as there might be more k-mers
			# that are after the real start.
		else:
			#Predicted k-mer ends after the real-kmer starts.
			startIndex = realStart;
			numFP = numFP + (realStart - predictedStart);

	return startIndex, numFP, numFN;

def getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer):
	numTP = 0;
	numFP = 0;

	predKmerIndex = 0;
	#print "Start Index: ", str(startIndex), ", End Index: ", str(endIndex)
	for index in range(startIndex, endIndex):
		if seq[index] == predKmer[predKmerIndex]:
			numTP = numTP + 1;
		else:
			numFP = numFP + 1;
		predKmerIndex = predKmerIndex + 1;
	#print "After KMER: TP: ", str(numTP), ", FP: " , str(numFP)

	return numTP, numFP;

def getNumbersAfterAllKmers(realStart, realEnd, predictedStart, predictedEnd):
	numFN = 0;
	numFP = 0;

	if realStart < realEnd:
		numFN = numFN + (realEnd - realStart)

	#if realEnd < predictedEnd:
	#	numFP = numFP + (predictedEnd - predictedStart)
	#	print "REST OF PREDICTED KMER: FP: " , str(numFP)

	return numFP, numFN;

def getKmerLengthFromREString(kmerREString):
	kmers = kmerREString.split('|')
	kmerStr = "".join([str(x) for x in kmers] )
	kmerStr = kmerStr[1:-1]
	return len(kmerStr);

def getPredictedInfo(m):
	return int(m.start()), int(m.end()), m.group(1)

def getNumbersForSeq(kmerREString, realStart, realEnd, seq):
	predictedStart = predictedEnd = 0;
	notFoundFP = 0;
	notFoundFN = 0;
	numTP = 0;
	numFN = 0;
	numFP = 0;

	for m in re.finditer(kmerREString, seq):

		predictedStart, predictedEnd, predKmer  = getPredictedInfo(m);

		#print "Predicted Kmer: ", predKmer, ", Predicted Start: ", str(predictedStart), ", End: ", str(predictedEnd)
		startIndex, startFP, startFN = getStartIndexAndUpdateNumbers(predictedStart, predictedEnd, realStart, realEnd, predKmer);

		numFN = numFN + startFN;
		numFP = numFP + startFP;

		#print "Based on Start Location: FP: " , str(startFP), ", FN: ", str(startFN)
		endIndex = getEndIndex(predictedEnd, realEnd)	

		# This can happen if the predicted kmer ends before the real kmer starts.
		if endIndex < startIndex:
			endIndex = startIndex;
	
		kmerTP, kmerFP = getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		
		numTP = numTP + kmerTP
		numFP = numFP + kmerFP

		# Update realStart Index after k-mer is parsed.
		if endIndex >= realStart:
			realStart = endIndex;
		
		
	# We did not find the predicted k-mer in the sequence.
	if predictedStart == 0 and predictedEnd == 0:
		notFoundFP = notFoundFP + getKmerLengthFromREString(kmerREString)

	# We did not find any predicted k-mers in the sequence. 
	if realStart < realEnd:
		numFN = numFN + (realEnd - realStart);

	numFP = numFP + notFoundFP
	numFN = numFN + notFoundFN

	return numTP, numFP, numFN;

def getEndIndex(predictedEnd, realEnd):
	if predictedEnd < realEnd:
		endIndex = predictedEnd
	else:
		endIndex = realEnd;

	return endIndex;

def getRealKmerDetails(realKmerDict, seqid):
	if seqid in realKmerDict:
		return realKmerDict[seqid][0], int(realKmerDict[seqid][1]), (int(realKmerDict[seqid][1]) + len(realKmerDict[seqid][0]))
	else:
		return "", 0, 0

def getKmerFromPSSM(pssmList, seq):
	kmerReToSearchFor = "("
	for pssmLines in pssmList:
		bestSeqKmerTuple = TAMO_Motif.GetKmerFromMotifFromPSSM(pssmLines, seq);
		kmer = bestSeqKmerTuple[1]
		kmerReToSearchFor = kmerReToSearchFor + kmer
		kmerReToSearchFor = kmerReToSearchFor + '|'
	kmerReToSearchFor = kmerReToSearchFor[:-1]
	kmerReToSearchFor = kmerReToSearchFor + ")"
	#print "PSSM RE: ", kmerReToSearchFor
	return kmerReToSearchFor;

def getKmerRE(predictedMotifs, seq):
	kmerReToSearchFor = "("
	for motif in predictedMotifs:
		bestSeqKmerTuple = TAMO_Motif.GetKmerFromTextMotifForSeq(motif, seq)
		kmer = bestSeqKmerTuple[1]
		kmerReToSearchFor = kmerReToSearchFor + kmer
		kmerReToSearchFor = kmerReToSearchFor + '|'
	kmerReToSearchFor = kmerReToSearchFor[:-1]
	kmerReToSearchFor = kmerReToSearchFor + ")"
	return kmerReToSearchFor;

def findKmerInSeqFromPSSM(realKmerDict, pssmList, posFile, negFile):
	PosSeqDict = SeqGenUtils.fasta_read(posFile);
	NegSeqDict = SeqGenUtils.fasta_read(negFile);
	numTotalFN = 0;
	numTotalTP = 0;
	numTotalFP = 0;	
	for seqid, seq in PosSeqDict.iteritems():
		realKmer, realStart, realEnd = getRealKmerDetails(realKmerDict, seqid);
		if realKmer == "" and realStart == 0 and realEnd == 0:
			print "SeqID: ", seqid , " does not contain kmer."
		kmerREString = getKmerFromPSSM(pssmList, seq)

		numTP, numFP, numFN = getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		#print seqid, ": TP: ", str(numTP), ", FP: ", str(numFP), ", FN: ", str(numFN)
		numTotalTP = numTotalTP + numTP;
		numTotalFP = numTotalFP + numFP;
		numTotalFN = numTotalFN + numFN;

	print "Pos File: TP: ", str(numTotalTP), ", FP: ", numTotalFP, ", FN: ", numTotalFN
	return numTotalTP, numTotalFP, numTotalFN;

def findKmerInSeq(realKmerDict, predictedMotifs, posFile, negFile):
	PosSeqDict = SeqGenUtils.fasta_read(posFile);
	NegSeqDict = SeqGenUtils.fasta_read(negFile);
	numTotalFN = 0;
	numTotalTP = 0;
	numTotalFP = 0;
	for seqid, seq in PosSeqDict.iteritems():
		realKmer, realStart, realEnd = getRealKmerDetails(realKmerDict, seqid);
		if realKmer == "" and realStart == 0 and realEnd == 0:
			print "SeqID: ", seqid , " does not contain kmer."
		kmerREString = getKmerRE(predictedMotifs, seq)
		numTP, numFP, numFN = getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		#print seqid, ": TP: ", str(numTP), ", FP: ", str(numFP), ", FN: ", str(numFN)
		numTotalTP = numTotalTP + numTP;
		numTotalFP = numTotalFP + numFP;
		numTotalFN = numTotalFN + numFN;
	
	print "Pos File: TP: ", str(numTotalTP), ", FP: ", numTotalFP, ", FN: ", numTotalFN
	return numTotalTP, numTotalFP, numTotalFN;
	

def CompareKmers(realCsvFile, predictedDremeFile, posFile, negFile, textMotif=0):
	realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);
	numTP = 0
	numFP = 0
	numFN = 0

	if textMotif == 1:
		predictedKmers = parseDreme.getPredictedDremeMotifs(predictedDremeFile);
		numTP, numFP, numFN = findKmerInSeq(realKmerDict, predictedKmers, posFile, negFile)
	else:
		pssmList = parseDreme.getPSSMListFromDremeFile(predictedDremeFile);
		numTP, numFP, numFN = findKmerInSeqFromPSSM(realKmerDict, pssmList, posFile, negFile)

	return numTP, numFP, numFN;

if __name__ == "__main__":
	import sys
	realCsvFile = sys.argv[1]
	predictedDremeFile= sys.argv[2]
	posFile = sys.argv[3]
	negFile = sys.argv[4]
	textMotif = int(sys.argv[5])

	numTP, numFP, numFN = CompareKmers(realCsvFile, predictedDremeFile, posFile, negFile, textMotif)
