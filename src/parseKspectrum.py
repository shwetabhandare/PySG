import sys
import re
import parseRealKmers

def read_file(dreme_file):
	with open (dreme_file, "r") as myfile:
		data=myfile.readlines()
	return data;

def findKmers(file_contents, maxKmers=None):
	kmerDict = dict();
	pattern2 = re.compile('[^-](\d+\.\d+)\,([ATGC]+)');
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

def GetUniqueRealKmers(realDict):
	realKmers = list();
	for seqid, value in realDict.iteritems():
		realKmers.append(value[0]);

	result = set (realKmers);
	return list(result);

def GetNumPredictedKmersFoundInReal(predictedDict, realDict):

	uniqueRealKmers = GetUniqueRealKmers(realDict);

	numPredictedKmersFound = 0;
	for predictedKmer, value in predictedDict.iteritems():
		#print "Looking for predicted kmer: ", predictedKmer, ", real kmers: ", uniqueRealKmers
		matching = [s for s in uniqueRealKmers if predictedKmer in s]
		if len(matching) > 0:
			numPredictedKmersFound = numPredictedKmersFound + 1;

	return numPredictedKmersFound;

if __name__ == "__main__":
	import sys
	kmerDict = FindKspectrumKmers(sys.argv[1]);
	#print kmerDict;
	realDict = parseRealKmers.GetRealKmerDict(sys.argv[2]);
	numPredictedKmersFound = GetNumPredictedKmersFoundInReal(kmerDict, realDict);
	print "Total Predicted Kmers: ", str(len(kmerDict)), ", found : ", str(numPredictedKmersFound);
