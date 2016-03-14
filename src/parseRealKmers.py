import csv

def GetRealKmerDict(realKmersCsvFile):
	realKmersDict = dict();
	with open(realKmersCsvFile, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			seq_id = row[0]
			kmer = row[1]
			location = row[2]
			realKmersDict[seq_id] = [kmer, location];
	return realKmersDict;

if __name__ == "__main__":
	import sys
	realKmersDict = GetRealKmerDict(sys.argv[1])
	for key, value in realKmersDict.iteritems():
		print key;
		for seqid, location in zip(value[0], value[1]):
			print seqid, location;