from __future__ import division
import sys
import re
import parseRealKmers

def read_file(dreme_file):
	with open (dreme_file, "r") as myfile:
		data=myfile.readlines()
	return data;

def findKmers(file_contents, maxKmers=None):
	kmerDict = dict();
	pattern2 = re.compile('[^-](\d+\.\d+)\,([HuR_|TTP_|ATGC]+)');
	kmerCount = 0;
	for match2 in pattern2.finditer(file_contents):
		kmer_score = match2.group(1);
		kmer = match2.group(2);

		# Take all the positive k-mers.
		if maxKmers == None:
			maxKmers = sys.maxint;
		if float(kmer_score) > 0.0 and kmerCount < maxKmers:
		#if float(kmer_score) > 0.0:
			#print "%s:%s" %(match2.group(1), match2.group(2))
			kmerDict[kmer] = kmer_score;
			kmerCount = kmerCount + 1;

	return kmerDict;

def findSpecificKmers(file_contents, ReString, maxKmers=None):
	kmerDict = dict();
	pattern2 = re.compile(ReString);
	kmerCount = 0;
	for match2 in pattern2.finditer(file_contents):
		kmer_score = match2.group(1);
		kmer = match2.group(2);

		# Take all the positive k-mers.
		if maxKmers == None:
			maxKmers = sys.maxint;
		if float(kmer_score) > 0.0 and kmerCount < maxKmers:
		#if float(kmer_score) > 0.0:
			#print "%s:%s" %(match2.group(1), match2.group(2))
			kmerDict[kmer] = kmer_score;
			kmerCount = kmerCount + 1;

	return kmerDict;


def FindKspectrumKmers(kspectrumFeatureFile, maxKmers=None):
	fileContents = str(read_file(kspectrumFeatureFile))
	return findKmers(fileContents, maxKmers);

def FindRBPSpecificKmers(kspectrumFeatureFile, ReString, maxKmers=None):
	fileContents = str(read_file(kspectrumFeatureFile))
	return findSpecificKmers(fileContents, ReString, maxKmers);

def GetUniqueRealKmers(realDict):
	realKmers = list();
	for seqid, value in realDict.iteritems():
		realKmers.append(value[0]);

	result = set (realKmers);
	return list(result);

def GetNumPredictedKmersFoundInReal(predictedDict, realDict):

	uniqueRealKmers = GetUniqueRealKmers(realDict);
	realKmersFound = list()

	numPredictedKmersFound = 0;
	for predictedKmer, value in predictedDict.iteritems():
		#print "Looking for predicted kmer: ", predictedKmer, ", real kmers: ", uniqueRealKmers
		matching = [s for s in uniqueRealKmers if predictedKmer in s]
		if len(matching) > 0:
			numPredictedKmersFound = numPredictedKmersFound + 1;
			realKmersFound.append(predictedKmer)

	print "Number of Predicted K-mers found in real-set:", numPredictedKmersFound
	return numPredictedKmersFound, realKmersFound;

def ComputeNucleotidesMatchedBetweenRealAndEmbeddedKmers(realKmers, kmersPredicted):
	numMatched = 0;
	numMismatched = 0;

	totalCompared = 0;
	for kmer in realKmers:
		totalCompared = totalCompared + len(kmer)
	print totalCompared;
	for predictedKmer in kmersPredicted:
		for realKmer in realKmers:
			u = zip(predictedKmer, realKmer)
			#print "Comparing predicted: ", predictedKmer, " with real: ", realKmer;
			for i, j in u:
				if i == j:
					#print i, '---', j
					numMatched = numMatched + 1;
				else:
					#print i, 'xxx', j
					numMismatched = numMismatched + 1;


	print "Matched percentage: ", numMatched
	print "Mis-matched percentage: ", numMismatched
	return numMatched, numMismatched;




if __name__ == "__main__":
	import sys
	
	kmerDict25 = FindKspectrumKmers(sys.argv[1], 25);
	kmerDict50 = FindKspectrumKmers(sys.argv[1], 50);
	kmerDict100 = FindKspectrumKmers(sys.argv[1], 100);

	kmer25Unique = set(kmerDict25)
	kmer50Unique = set(kmerDict50) - set(kmerDict25)
	kmer100Unique = set(kmerDict100) - set(kmerDict50)

	realDict = parseRealKmers.GetRealKmerDict(sys.argv[2]);
	realKmers = [x[0] for x in realDict.values()]

	#ComputeNucleotidesMatchedBetweenRealAndEmbeddedKmers(realKmers, kmer25Unique)
	#ComputeNucleotidesMatchedBetweenRealAndEmbeddedKmers(realKmers, kmer50Unique)
	#ComputeNucleotidesMatchedBetweenRealAndEmbeddedKmers(realKmers, kmer100Unique)

	#print "Embedded K-mers:", list(set(realKmers))
	print "Predicted K-mers : 25 ", list(kmer25Unique)
	print "Predicted K-mers : 50 ", list(kmer50Unique)
	print "Predicted K-mers : 100 ", list(kmer100Unique)

	numPredictedKmersFound25, predictedKmersFound25 = GetNumPredictedKmersFoundInReal(kmerDict25, realDict);
	print predictedKmersFound25
	numPredictedKmersFound50, predictedKmersFound50 = GetNumPredictedKmersFoundInReal(kmerDict50, realDict);
	numPredictedKmersFound100, predictedKmersFound100 = GetNumPredictedKmersFoundInReal(kmerDict100, realDict);
	#print "Total Predicted Kmers: ", str(len(kmerDict)), ", found : ", str(numPredictedKmersFound);
