import sys
import re

def read_file(dreme_file):
	with open (dreme_file, "r") as myfile:
		data=myfile.readlines()
	return data;

def findKmers(file_contents):
	kmerDict = dict();
	pattern2 = re.compile('#\s(?:BEST\s+|\s+)([ATGC]+)\s+([ATGC]+)\s+(\d+)\s+(\d+)')
	for match2 in pattern2.finditer(file_contents):
		kmer = match2.group(1);
		kmer_rc = match2.group(2);
		kmer_pos_count = match2.group(3);
		kmer_neg_count = match2.group(4);
		#print "%s:%s:%s:%s" %(match2.group(1), match2.group(2), match2.group(3), match2.group(4))
		kmerDict[kmer] = [kmer_rc, kmer_pos_count, kmer_neg_count];
	return kmerDict;

def getPredictedDremeMotifs(file_contents):
	motifs = []
	pattern = re.compile(r"MOTIF\s+([ACGTURYSWKMBDHVN]+)")
	for match in pattern.finditer(str([file_contents])):
		motifs.append(match.group(1))

	print motifs
	return motifs

def FindDremeKmers(dremeResultFile):
	fileContents = str(read_file(dremeResultFile))
	return findKmers(fileContents);

if __name__ == "__main__":
	import sys
	#kmerDict = FindDremeKmers(sys.argv[1]);
	fileContents = read_file(sys.argv[1])
	FindDremeMotifs(fileContents)
	#for key, value in kmerDict.iteritems():
	#	print key, value