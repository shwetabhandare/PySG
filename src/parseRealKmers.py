import csv

def GetRealKmerDict(realKmersCsvFile):
	realKmersDict = dict();
	with open(realKmersCsvFile, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			key = row[0]
			kmer = row[1]
			location = row[2]
			realKmersDict[key] = [kmer, location];
	return realKmersDict;

if __name__ == "__main__":
	import sys
	realKmersDict = GetRealKmerDict(sys.argv[1])
