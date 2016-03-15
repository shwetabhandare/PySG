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
			startIndex = realEnd;
			numFP = numFP + len(predKmer)
		else:
			#Predicted k-mer ends after the real-kmer starts.
			startIndex = realStart;
			numFP = numFP + (realStart - predictedStart);

	return startIndex, numFP, numFN;

def getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer):
	numTP = 0;
	numFP = 0;

	predKmerIndex = 0;
	print "Start Index: ", str(startIndex), ", End Index: ", str(endIndex)
	for index in range(startIndex, endIndex):
		if seq[index] == predKmer[predKmerIndex]:
			numTP = numTP + 1;
		else:
			numFP = numFP + 1;
		predKmerIndex = predKmerIndex + 1;
	print "After KMER: TP: ", str(numTP), ", FP: " , str(numFP)

	return numTP, numFP;

def getNumbersAfterAllKmers(realStart, realEnd, predictedStart, predictedEnd):
	numFN = 0;
	numFP = 0;

	if realStart < realEnd:
		numFN = numFN + (realEnd - realStart)
		print "REST OF KMER: FN: ", str(numFN)

	#if realEnd < predictedEnd:
	#	numFP = numFP + (predictedEnd - predictedStart)
	#	print "REST OF PREDICTED KMER: FP: " , str(numFP)

	return numFP, numFN;

def getKmerLengthFromREString(kmerREString):
	kmers = kmerREString.split('|')
	kmerStr = "".join([str(x) for x in kmers] )
	kmerStr = kmerStr[1:-1]
	return len(kmerStr);

def getNumbersForSeq(kmerREString, realStart, realEnd, seq):
	predictedStart = predictedEnd = 0;
	notFoundFP = 0;
	notFoundFN = 0;
	numTP = 0;
	numFN = 0;
	numFP = 0;

	print "KMER RE STRING: ", kmerREString;
	for m in re.finditer(kmerREString, seq):

		predictedStart = int(m.start())
		predictedEnd = int(m.end())
		predKmer = m.group(1)

		print "Predicted Kmer: ", predKmer, ", Predicted Start: ", str(predictedStart), ", End: ", str(predictedEnd)
		startIndex, startFP, startFN = getStartIndexAndUpdateNumbers(predictedStart, predictedEnd, realStart, realEnd, predKmer);

		numFN = numFN + startFN;
		numFP = numFP + startFP;

		print "Based on Start Location: FP: " , str(startFP), ", FN: ", str(startFN)
		endIndex = getEndIndex(predictedEnd, realEnd)	
	

		kmerTP, kmerFP = getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		print "After Kmer comparison: TP: ", str(kmerTP), ", FP: ", str(kmerFP);
		numTP = numTP + kmerTP

		#No more k-mers found, and you want to handle the case for the rest of the real kmer
		if endIndex == realEnd:
			numFPAfterKmer, numFNAfterKmer = getNumbersAfterAllKmers(realStart, realEnd, predictedStart, predictedEnd)
			numFP = numFP + numFPAfterKmer
			numFN = numFN + numFNAfterKmer

		# Update realStart Index after k-mer is parsed.
		if endIndex > realStart:
			realStart = endIndex;


	# We did not find the predicted k-mer in the sequence.
	if predictedStart == 0 and predictedEnd == 0:
		notFoundFP = notFoundFP + getKmerLengthFromREString(kmerREString)
		notFoundFN = notFoundFN + (realEnd - realStart)

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
	return realKmerDict[seqid][0], int(realKmerDict[seqid][1]), (int(realKmerDict[seqid][1]) + len(realKmerDict[seqid][0]))

def getKmerRE(predictedMotifs, seq):
	print predictedMotifs
	kmerReToSearchFor = "("
	for motif in predictedMotifs:
		bestSeqKmerTuple = TAMO_Motif.GetKmerForSeq(motif, seq)
		kmer = bestSeqKmerTuple[1]
		kmerReToSearchFor = kmerReToSearchFor + kmer
		kmerReToSearchFor = kmerReToSearchFor + '|'
	kmerReToSearchFor = kmerReToSearchFor[:-1]
	kmerReToSearchFor = kmerReToSearchFor + ")"
	return kmerReToSearchFor;

def findKmerInSeq(realKmerDict, predictedMotifs, posFile, negFile):
	PosSeqDict = SeqGenUtils.fasta_read(posFile);
	NegSeqDict = SeqGenUtils.fasta_read(negFile);
	numFN = 0;
	numTP = 0;
	numFP = 0;
	for seqid, seq in PosSeqDict.iteritems():
		realKmer, realStart, realEnd = getRealKmerDetails(realKmerDict, seqid);
		print "Real KMER: ", realKmer;
		kmerREString = getKmerRE(predictedMotifs, seq)
		numTP, numFP, numFN = getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		print "SeqID" , seqid, "Numbers: TP: ", str(numTP), ", FP: ", str(numFP), ", FN: ", str(numFN)
		
			

if __name__ == "__main__":
	import sys
	realCsvFile = sys.argv[1]
	predictedDremeFile= sys.argv[2]
	posFile = sys.argv[3]
	negFile = sys.argv[4]

	realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);
	predictedKmers = parseDreme.getPredictedDremeMotifs(predictedDremeFile);

	#print "Predicted: ", predictedKmers, "\nReal: ", realKmerDict.keys()
	findKmerInSeq(realKmerDict, predictedKmers, posFile, negFile)
