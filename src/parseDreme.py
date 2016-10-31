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

def getPSSMListFromDremeFile(dremeFile):
	pattern = re.compile('letter-probability matrix: alength=\s+(\d+)\s+w=\s+(\d+)');
	fileContents = str(read_file(dremeFile))
	fileOffset = 0;
	pssmLinesList = [];
	numPSSMs = 0;
	for match in pattern.finditer(fileContents[fileOffset:]):

		if numPSSMs > 5:
			break;

		lengthOfPSSM = match.group(2)
		fileOffset = match.end();
		numPSSMs = numPSSMs + 1;


		pssmRe = "(\d.\d+\s+\d.\d+\s+\d.\d+\s+\d.\d+)"
		pattern2 = re.compile(pssmRe);
		count = 0;
		pssmLines = ""
		for match2 in pattern2.finditer(fileContents[fileOffset:]):
			fileOffset = match2.end();
			pssmLine = match2.group(1)

			if count < int(lengthOfPSSM):
				count = count + 1;
				pssmLines = pssmLines + pssmLine + "\n";

			if count >= int(lengthOfPSSM):
				pssmLines = pssmLines[:-1];
				pssmLinesList.append(pssmLines);
				break;

	#print "Total PSSMs: ", numPSSMs;
	return pssmLinesList;

def getPredictedDremeMotifs(dremeFile):
	fileContents = read_file(dremeFile)
	motifs = []
	pattern = re.compile(r"MOTIF\s+([ACGTURYSWKMBDHVN]+)")
	for match in pattern.finditer(str([fileContents])):
		motifs.append(match.group(1))
	return motifs

def FindDremeKmers(dremeResultFile):
	fileContents = str(read_file(dremeResultFile))
	return findKmers(fileContents);

if __name__ == "__main__":
	import sys
	dremeFile = sys.argv[1]
	#kmerDict = FindDremeKmers(sys.argv[1]);
	
	#getPredictedDremeMotifs(dremeFile)
	getPSSMListFromDremeFile(dremeFile)
	#for key, value in kmerDict.iteritems():
	#	print key, value
