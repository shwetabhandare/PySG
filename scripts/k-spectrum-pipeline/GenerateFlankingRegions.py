import yaml
import sys
import csv
import fasta
from fasta import *
import operator
from itertools import *
import math

def ReadConfigFile():
	confMap = {}
	configFile = sys.argv[1]
	f = open(configFile)
	confMap = yaml.load(f)
	f.close()
	return confMap;

def fasta_read(file_name):
	"""read the sequence from a file in fasta format"""
	seq_dict = dict()
	for record in fasta.fasta_itr(file_name):
		header = record.header
		seq_dict[header] = record.sequence;

						
	return seq_dict;

def CreateKmerDict(kmerFile):
	with open(kmerFile, 'rb') as csvFile:
		reader = csv.reader(csvFile, delimiter=',');
		numFeatures = 0;
		kmerDict = dict()
		for row in reader:
			if numFeatures <= 20:
				print row[0]
				featureScore = float(row[0])
				featureKmer = row[1]
				kmerDict[featureKmer] = featureScore;
				numFeatures = numFeatures + 1;
	return kmerDict;

def CreateFlankingRegion(sequence, kmer, flankStart, kmerIndex, flankEnd):
	l = list();

	flankStartRegion = sequence[flankStart:kmerIndex]
	flankEndRegion = sequence[kmerIndex + len(kmer):flankEnd]
	flankingRegion = sequence[flankStart:flankEnd]
	flankingRegionConcatStr = flankStartRegion + flankEndRegion
	
	l = [flankStartRegion, flankEndRegion, flankingRegionConcatStr, flankingRegion]
	return l;

def GetFlankingData(kmer, seq, kmerIndex):

	if kmerIndex > 0:
		flankStart = kmerIndex - 15 
		flankEnd   = kmerIndex + len(kmer) + 15
		seqLen = len(seq)
		if flankStart > 0 and flankEnd < seqLen:
			return CreateFlankingRegion(seq, kmer, flankStart, kmerIndex, flankEnd)										
		else:
			#print "Found a boundary condition "
			if flankStart <= 0:
				flankStart = 0
			if flankEnd > seqLen:
				flankEnd = seqLen
			return CreateFlankingRegion(seq, kmer, flankStart, kmerIndex, flankEnd)


def GetKmerCounts(kmerDict, seqDict):
	kmerCountDict = dict()

	for header, seq in seqDict.iteritems():
		for kmer, value in kmerDict.iteritems():
			#print "Looking for kmer: ",kmer, " in sequence: ", seq
			# Return the number of (non-overlapping) occurrences of substring sub in string s[start:end].
			kmerCount = seq.count(kmer);
			if kmer in kmerCountDict.keys():
				kmerCountDict[kmer] = kmerCountDict[kmer] + kmerCount;
			else:
				kmerCountDict[kmer] = kmerCount;	

	return kmerCountDict;

	
def GetSortedKmerDict(kmerCountDict):

	sorted_kmerCountDict = [ (v,k) for k,v in kmerCountDict.iteritems() ]
	sorted_kmerCountDict.sort(reverse=True)

	return sorted_kmerCountDict;

def WriteTopLogLikelihoodKmers(topLikelyKmers, topKmersFile):
	topKmersCsvWriter = csv.writer(open(topKmersFile, 'w'));
	for item in topLikelyKmers:
		row = list()
		row = [item[0], item[1]];
		print row
		topKmersCsvWriter.writerow(row);


def CreateKmerCountAndFlankingDict(kmerDict, seqDict, flankingRegionFile, topKmersFile):
	kmerCountDict = dict()
	kmerFlankingDict = dict()

	csvWriter = csv.writer(open(flankingRegionFile, 'w'));
	topKmersCsvWriter = csv.writer(open(topKmersFile, 'w'));

	for header, seq in seqDict.iteritems():
		for kmer, value in kmerDict.iteritems():
			print "Looking for kmer: ",kmer, " in sequence: ", seq
			kmerCount = seq.count(kmer);
			if kmer in kmerCountDict.keys():
				kmerCountDict[kmer] = kmerCountDict[kmer] + kmerCount;
			else:
				kmerCountDkict[kmer] = kmerCount;

			kmerIndex = seq.find(kmer)
			if kmerIndex > 0:
				l = GetFlankingData(kmer, seq, kmerIndex);
				kmerFlankingDict[(header, kmer)] = l
				row = list()
				row = [header, kmer, l[0], l[1]];
				csvWriter.writerow(row);
	#for key, value in kmerInSeqDict.iteritems():
		#print key, value;
	sorted_kmerCountDict = [ (v,k) for k,v in kmerCountDict.iteritems() ]
	sorted_kmerCountDict.sort(reverse=True)
	for key, value in sorted_kmerCountDict:
		row = list()
		row = [key, value];
		print row
		topKmersCsvWriter.writerow(row);


def CreateFlankingRegions(posSeq, negSeq, kmerFile, topKmersFile):

	kmerDict = CreateKmerDict(kmerFile)
	
	posSeqDict = fasta_read(posSeq)
	negSeqDict = fasta_read(negSeq)

	#CreateKmerCountAndFlankingDict(kmerDict, seqDict, flankingRegionFile, topKmersFile);
	posKmerCountDict = GetKmerCounts(kmerDict, posSeqDict);
	negKmerCountDict = GetKmerCounts(kmerDict, negSeqDict);

	for k, v in negKmerCountDict.iteritems():
		print k, v;

	sortedPosKmerCountList = GetSortedKmerDict(posKmerCountDict)
	sortedNegKmerCountList = GetSortedKmerDict(negKmerCountDict)

	LogLikeHood = dict();

	for pos_item in sortedPosKmerCountList:
		for neg_item in sortedNegKmerCountList:
			kmer = pos_item[1]
			if kmer == neg_item[1]:
				positiveCount = pos_item[0]
				negativeCount = neg_item[0]
				print positiveCount, negativeCount

				if positiveCount == 0:
					LogLikeHood[kmer] = 0;
				elif negativeCount == 0:
					print "Could not find ", kmer, " in negative set."
					LogLikeHood[kmer] = math.log10(positiveCount);
				else:
					LogLikeHood[kmer] = math.log10( float(positiveCount)/ float(negativeCount))

	sortedLogLikeHood =  [ (v,k) for k,v in LogLikeHood.iteritems() ]
	sortedLogLikeHood.sort(reverse=True)

	WriteTopLogLikelihoodKmers(sortedLogLikeHood, topKmersFile);

if __name__ == '__main__':
	import sys;
	posSeq = sys.argv[1]
	negSeq = sys.argv[2]
	kmerFile = sys.argv[3]
	topKmersFile = sys.argv[4]
	CreateFlankingRegions(posSeq, negSeq, kmerFile, topKmersFile);
