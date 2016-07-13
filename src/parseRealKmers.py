import csv

def GetRealKmerDict(realKmersCsvFile):
	realKmersDict = dict();
	with open(realKmersCsvFile, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if len(row) < 3:
				print "Found invalid row:", row;
			else:
				seq_id = row[0]
				kmer = row[1]
				location = row[2]
				realKmersDict[seq_id] = [kmer, location];
	return realKmersDict;

if __name__ == "__main__":
	import sys
	realKmersDict = GetRealKmerDict(sys.argv[1])
	embeddedKmers = list();
	for key, value in realKmersDict.iteritems():
		#print key, value;
		kmer = value[0]
		embeddedKmers.append(kmer);
	print "Unique Kmers embedded: ", set(embeddedKmers)
