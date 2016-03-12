#import re, ahocorasick,random,time
from parseDreme import *
from parseRealKmers import *
import fasta
import SeqGenUtils
import re
#from itertools import imap

def getPredictedKmerDict(dremeResultFile):
	return FindDremeKmers(dremeResultFile);

def getReadKmerDict(realKmersCsvFile):
	return GetRealKmerDict(realKmersCsvFile)


def compareKmers(realKmerDict, predictedKmerDict, posFile, negFile):
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
	compareKmers(realKmerDict, predictedKmerDict, posFile, negFile)
