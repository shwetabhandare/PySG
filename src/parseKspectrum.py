import sys
import re

def read_file(dreme_file):
	with open (dreme_file, "r") as myfile:
		data=myfile.readlines()
	return data;

def findKmers(file_contents):
	kmerDict = dict();
	pattern2 = re.compile('(\d+\.\d+)\,([ATGC]+)')
	kmerCount = 0;
	for match2 in pattern2.finditer(file_contents):
		kmer_score = match2.group(1);
		kmer = match2.group(2);

		if kmerCount <= 20:
			#print "%s:%s" %(match2.group(1), match2.group(2))
			kmerDict[kmer] = kmer_score;
			kmerCount = kmerCount + 1;

	return kmerDict;

def FindKspectrumKmers(kspectrumFeatureFile):
	fileContents = str(read_file(kspectrumFeatureFile))
	return findKmers(fileContents);

if __name__ == "__main__":
	import sys
	kmerDict = FindKspectrumKmers(sys.argv[1]);
	for key, value in kmerDict.iteritems():
		print key, value
