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

	for seqid, seq in PosSeqDict.iteritems():
		numTN = 0;
		numTP = 0;
		numFP = 0;
		realStart = int(realKmerDict[seqid][1])
		realEnd = int(realKmerDict[seqid][1]) + len(realKmerDict[seqid][0])
		for predKmer, value in predictedKmerDict.iteritems():
			realKmerIndex = 0;
			for m in re.finditer(predKmer, seq):
				print "Seq ID: ", seqid, ", : ", predKmer, " at: ", m.start(), m.end(), ", real kmer at: ", realStart, ", ", realEnd
				if m.start() > realStart:
					startIndex = realStart;
					numTN = m.start() - realStart;
				else:
					startIndex = m.start();
					numTP = m.start() - startIndex;

				if m.end() < realEnd:
					endIndex = m.end()
				else:
					endIndex = realEnd;

				# Now go through the predicted kmer, and compare with real kmer
				predKmerIndex = 0;
				for index in startIndex, endIndex:
					if seq[index] == predKmer[predKmerIndex]:
						numTP = numTP + 1;
					else:
						numFP = numFP + 1;
					predKmerIndex = predKmerIndex + 1;

				realStart = endIndex;

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
