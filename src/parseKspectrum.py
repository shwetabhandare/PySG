import sys
import re

def read_file(dreme_file):
	with open (dreme_file, "r") as myfile:
		data=myfile.readlines()
	return data;

def findKmers(file_contents):
	kmerDict = dict();
	pattern2 = re.compile('(\d+\.\d+)\,([ATGC]+)')
	for match2 in pattern2.finditer(file_contents):
		kmer_score = match2.group(1);
		kmer = match2.group(2);
		if float(kmer_score) >= 1.0:
			#print "%s:%s" %(match2.group(1), match2.group(2))
			kmerDict[kmer] = kmer_score;

	return kmerDict;

def FindDremeKmers(dremeResultFile):
	fileContents = str(read_file(dremeResultFile))
	return findKmers(fileContents);

if __name__ == "__main__":
	import sys
	kmerDict = FindDremeKmers(sys.argv[1]);
	#for key, value in kmerDict.iteritems():
		#print key, value
