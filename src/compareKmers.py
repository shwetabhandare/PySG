#import re, ahocorasick,random,time
from parseDreme import *
from parseRealKmers import *
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

	if realEnd < predictedEnd:
		numFP = numFP + (predictedEnd - realEnd)
		print "REST OF PREDICTED KMER: FP: " , str(numFP)

	return numFP, numFN;

def getNumbersForSeq(realKmer, realStart, realEnd, predKmer, seq):
	predictedStart = predictedEnd = 0;
	notFoundFP = 0;
	numTP = 0;
	numFN = 0;
	numFP = 0;

	for m in re.finditer(predKmer, seq):
		predictedStart = int(m.start())
		predictedEnd = int(m.end())

		print "Predicted Start: ", str(predictedStart), ", End: ", str(predictedEnd)
		startIndex, startFP, startFN = getStartIndexAndUpdateNumbers(predictedStart, predictedEnd, realStart, realEnd, predKmer);

		numFN = numFN + startFN;
		numFP = numFP + startFP;

		print "Based on Start Location: FP: " , str(startFP), ", FN: ", str(startFN)
		endIndex = getEndIndex(predictedEnd, realEnd)	
	
		kmerTP, kmerFP = getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		print "After Kmer comparison: TP: ", str(kmerTP), ", FP: ", str(kmerFP);
		numTP = numTP + kmerTP
		
		# Update realStart Index after k-mer is parsed.
		if endIndex > realStart:
			realStart = endIndex;


	# We did not find the predicted k-mer in the sequence.
	if predictedStart == 0 and predictedEnd == 0:
		notFoundFP = notFoundFP + len(predKmer)

	numFPAfterKmer, numFNAfterKmer = getNumbersAfterAllKmers(realStart, realEnd, predictedStart, predictedEnd)	

	numFP = numFP + notFoundFP + numFPAfterKmer;
	numFN = numFN + numFNAfterKmer;

	print "Total numbers: TP: ", str(numTP), ", FP: ", str(numFP), ", FN: ", str(numFN)

	return numTP, numFP, numFN;

def getEndIndex(predictedEnd, realEnd):
	if predictedEnd < realEnd:
		endIndex = predictedEnd
	else:
		endIndex = realEnd;

	return endIndex;

def getRealKmerDetails(realKmerDict, seqid):
	return realKmerDict[seqid][0], int(realKmerDict[seqid][1]), (int(realKmerDict[seqid][1]) + len(realKmerDict[seqid][0]))

def findKmerInSeq(realKmerDict, predictedMotifs, posFile, negFile):
	PosSeqDict = SeqGenUtils.fasta_read(posFile);
	NegSeqDict = SeqGenUtils.fasta_read(negFile);
	numFN = 0;
	numTP = 0;
	numFP = 0;
	for seqid, seq in PosSeqDict.iteritems():
		realKmer, realStart, realEnd = getRealKmerDetails(realKmerDict, seqid);
		for motif in predictedMotifs:
			kmer = TAMO_Motif.GetKmerForSeq(motif, seq)
			numTP, numFP, numFN = getNumbersForSeq(realKmer, realStart, realEnd, kmer, seq);
			
		print "End of KMER: TP: ", str(numTP), ", FP: " , str(numFP), ", FN: ", str(numFN)
		
			
def compareRealAndPredictedKmers(realKmerDict, predictedKmerDict, posFile, negFile):
	PosSeqDict = SeqGenUtils.fasta_read(posFile);
	NegSeqDict = SeqGenUtils.fasta_read(negFile);
	numFN = 0;
	numTP = 0;
	numFP = 0;
	for seqid, seq in PosSeqDict.iteritems():
		realStart = int(realKmerDict[seqid][1])
		realEnd = int(realKmerDict[seqid][1]) + len(realKmerDict[seqid][0])
		print "***** SeqID: ", seqid, ", KMER" , realKmerDict[seqid][0],  ", Real Start:", str(realStart), ", Real End: ", str(realEnd)
		for predKmer, value in predictedKmerDict.iteritems():
			print "Predicted KMER: ", predKmer
			predictedStart = predictedEnd = 0;
			for m in re.finditer(predKmer, seq):
				predictedStart = int(m.start())
				predictedEnd = int(m.end())
				print "Predicted Start: ", str(predictedStart), ", End: ", str(predictedEnd)

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
						#Predicted k-mer ends before the real k-mer starts
						startIndex = predictedEnd;
						numFP = numFP + len(predKmer)
					else:
						#Predicted k-mer ends after the real-kmer starts.
						startIndex = realStart;
						numFP = numFP + (realStart - predictedStart);
				if predictedEnd < realEnd:
					endIndex = predictedEnd
				else:
					endIndex = realEnd;

				print "Based on Start Location: TP: ", str(numTP), ", FP: " , str(numFP), ", FN: ", str(numFN)
				# Now go through the predicted kmer, and compare with real kmer
				predKmerIndex = 0;
				print "Start Index: ", str(startIndex), ", End Index: ", str(endIndex)
				for index in range(startIndex, endIndex):
					print "Seq Index: ", seq[index], ", Pred index: ", predKmer[predKmerIndex]
					if seq[index] == predKmer[predKmerIndex]:
						numTP = numTP + 1;
					else:
						numFP = numFP + 1;
					predKmerIndex = predKmerIndex + 1;
				print "After KMER: TP: ", str(numTP), ", FP: " , str(numFP), ", FN: ", str(numFN)

				if endIndex > realStart:
					realStart = endIndex;

			if predictedStart == 0 and predictedEnd == 0:
				# We did not find the predicted k-mer in the sequence.
				numFP = numFP + len(predKmer)
				print "NOT FOUND KMER: TP: ", str(numTP), ", FP: " , str(numFP), ", FN: ", str(numFN)
			#Handle the case where we have nucleotides left after going 
			#through the predicted k-mers.	
		if realStart < realEnd:
			numFN = numFN + (realEnd - realStart)
			print "REST OF KMER: TP: ", str(numTP), ", FP: " , str(numFP), ", FN: ", str(numFN)

		if realEnd < predictedEnd:
			numFP = numFP + (predictedEnd - realEnd)
			print "REST OF PREDICTED KMER: TP: ", str(numTP), ", FP: " , str(numFP), ", FN: ", str(numFN)

		print "End of KMER: TP: ", str(numTP), ", FP: " , str(numFP), ", FN: ", str(numFN)

def aho():
	# search N words from dict
	N=3

	#file from http://norvig.com/big.txt
	with open("big.txt","r") as f:
		text = f.read()

		#words = set(re.findall('[atgc]+', text.lower())) 
		foundKmerDict = FindDremeKmers("/projects/bhandare/workspace/PySG/src/Test/dreme.txt");
		#search_words = random.sample([w for w in words],N)
		search_words = foundKmerDict.keys();

		A = ahocorasick.Automaton()
		for i,w in enumerate(search_words):
			A.add_word(w, (i, w))

		A.make_automaton()
		#test time for ahocorasic
		start = time.time()
		print("ah matches",sum(1 for i in A.iter(text))) 
		print("aho done in ", time.time() - start)


		exp = re.compile('|'.join(search_words))
		#test time for re
		start = time.time()
		m = exp.findall(text)
		print("re matches",sum(1 for _ in m))
		print("re done in ",time.time()-start)

if __name__ == "__main__":
	import sys
	realCsvFile = sys.argv[1]
	predictedDreme = sys.argv[2]
	posFile = sys.argv[3]
	negFile = sys.argv[4]

	realKmerDict = getReadKmerDict(realCsvFile);
	predictedKmerDict = getPredictedKmerDict(predictedDreme);

	#print "Predicted: ", predictedKmerDict.keys(), "\nReal: ", realKmerDict.keys()
	compareRealAndPredictedKmers(realKmerDict, predictedKmerDict, posFile, negFile)
