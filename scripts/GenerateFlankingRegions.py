import yaml
import sys
import csv
import fasta
from fasta import *
import operator

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
			if numFeatures <= 50:
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



def CreateKmerCountAndFlankingDict(kmerDict, seqDict, flankingRegionFile, topKmersFile):
	kmerInSeqDict = dict()
	kmerCountDict = dict()
	kmerFlankingDict = dict()

	csvWriter = csv.writer(open(flankingRegionFile, 'w'));
	topKmersCsvWriter = csv.writer(open(topKmersFile, 'w'));

	for header, seq in seqDict.iteritems():
		for kmer, value in kmerDict.iteritems():
			kmerCount = seq.count(kmer);
			kmerInSeqDict[(header, kmer)] = kmerCount;
			if kmer in kmerCountDict.keys():
				kmerCountDict[kmer] = kmerCountDict[kmer] + kmerCount;
			else:
				kmerCountDict[kmer] = kmerCount;

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
		topKmersCsvWriter.writerow(row);


def CreateFlankingRegions():
	confMap = ReadConfigFile();
	seqFile = confMap['input']['dataFile']
	kmerFile = confMap['input']['featureKmers']
	flankingRegionFile = confMap['output']['flankingRegionFile']
	topKmersFile = confMap['output']['topKmersFile']
	kmerDict = CreateKmerDict(kmerFile)
	seqDict = fasta_read(seqFile)
	CreateKmerCountAndFlankingDict(kmerDict, seqDict, flankingRegionFile, topKmersFile);



if __name__ == '__main__':
	CreateFlankingRegions();
