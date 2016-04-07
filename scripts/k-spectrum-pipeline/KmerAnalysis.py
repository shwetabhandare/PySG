import sys
import re
import fasta
import GenerateFlankingRegions


def getSeqDict(fasta_file):
	return GenerateFlankingRegions.fasta_read(fasta_file);

def createKmerStringRE(kmerList):
	kmerReToSearchFor = "("
	for kmer in kmerList:
		kmerReToSearchFor = kmerReToSearchFor + kmer
		kmerReToSearchFor = kmerReToSearchFor + '|'
	kmerReToSearchFor = kmerReToSearchFor[:-1]
	kmerReToSearchFor = kmerReToSearchFor + ")"
	#print "Kmer RE: ", kmerReToSearchFor
	return kmerReToSearchFor;

def read_file(kmerFeatureFile):
	with open (kmerFeatureFile, "r") as myfile:
		data=myfile.readlines()
	return data;

def findKmers(file_contents):
	kmerDict = dict();
	pattern2 = re.compile('(\d+\.\d+)\,([HuR_|TTP_|ATGC]+)')
	kmerCount = 0
	for match2 in pattern2.finditer(file_contents):
		kmer_score = match2.group(1);
		kmer = match2.group(2);
		if float(kmer_score) >= 1.0 and kmerCount <= 50:
			#print "%s:%s" %(match2.group(1), match2.group(2))
			kmerDict[kmer] = kmer_score;
			kmerCount = kmerCount + 1;

	return kmerDict;

def getNumKmersFoundInDict(KmerList, SeqDict):
	kmerREString = createKmerStringRE(KmerList)
	numFoundInDict = 0;

	for seqid, seq in SeqDict.iteritems():
		for m in re.finditer(kmerREString, seq):
			if int(m.start()) != int(m.end()):
				#print "Found: ", m.start(), m.end(), m.group(1)
				numFoundInDict = numFoundInDict + 1;

	return numFoundInDict;

def checkKmersInSeqDict(HuR_Kmers, TTP_Kmers, Common_Kmers, Common_Seq_Dict, HuROnly_Seq_Dict, TTP_Only_Seq_Dict):

	hurFoundInHuR = getNumKmersFoundInDict(HuR_Kmers, HuROnly_Seq_Dict);
	hurFoundInTTP = getNumKmersFoundInDict(HuR_Kmers, TTP_Only_Seq_Dict);
	hurFoundInCommon = getNumKmersFoundInDict(HuR_Kmers, Common_Seq_Dict);

	ttpFoundInHuR = getNumKmersFoundInDict(TTP_Kmers, HuROnly_Seq_Dict);
	ttpFoundInTTP = getNumKmersFoundInDict(TTP_Kmers, TTP_Only_Seq_Dict);
	ttpFoundInCommon = getNumKmersFoundInDict(TTP_Kmers, Common_Seq_Dict);

	commonFoundInHuR = getNumKmersFoundInDict(Common_Kmers, HuROnly_Seq_Dict);
	commonFoundInTTP = getNumKmersFoundInDict(Common_Kmers, TTP_Only_Seq_Dict);
	commonFoundInCommon = getNumKmersFoundInDict(Common_Kmers, Common_Seq_Dict);


	print "HuR Kmers: Only HuR: ", str(hurFoundInHuR), ", Only TTP: ",str(hurFoundInTTP), ", Common: ", str(hurFoundInCommon)
	print "TTP Kmers: Only HuR: ", str(ttpFoundInHuR), ", Only TTP: ",str(ttpFoundInTTP), ", Common: ", str(ttpFoundInCommon)
	print "Common Kmers: Only HuR: ", str(commonFoundInHuR), ", Only TTP: ",str(commonFoundInTTP), ", Common: ", str(commonFoundInCommon)

def sortKmers(kmerDict):
	HuR_Kmers = list();
	TTP_Kmers = list();
	Common_Kmers = list();

	for key, value in kmerDict.iteritems():
		if "HuR" in key:
			HuR_Kmers.append(key[4:])
		elif "TTP" in key:
			TTP_Kmers.append(key[4:])
		else:
			Common_Kmers.append(key)

	return HuR_Kmers, TTP_Kmers, Common_Kmers

if __name__ == "__main__":

	kmerFeatureFile = sys.argv[1]

	commonFaFile = sys.argv[2]
	Only_HuR_File = sys.argv[3]
	Only_TTP_File = sys.argv[4]

	kmerFeatureContents = str(read_file(kmerFeatureFile))
	kmerDict = findKmers(kmerFeatureContents)
	HuR_Kmers, TTP_Kmers, Common_Kmers = sortKmers(kmerDict)

	Common_Seq_Dict = getSeqDict(commonFaFile)
	HuROnly_Seq_Dict = getSeqDict(Only_HuR_File)
	TTP_Only_Seq_Dict = getSeqDict(Only_TTP_File)

	checkKmersInSeqDict(HuR_Kmers, TTP_Kmers, Common_Kmers, Common_Seq_Dict, HuROnly_Seq_Dict, TTP_Only_Seq_Dict)



