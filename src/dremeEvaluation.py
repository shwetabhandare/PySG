import sys
import parseDreme
import SeqGenUtils
import compareKmerCommon


def FindBestScoreSeqForMotif(seqDict, pssmList):
	for header, sequence in seqDict.iteritems():
		kmerREString = compareKmerCommon.getKmerFromPSSM(pssmList, sequence)
		print kmerREString


if __name__ == "__main__":
	import sys
	dremeFile = sys.argv[1]
	seqFile = sys.argv[2]
	seqDict  = SeqGenUtils.fasta_read(seqFile);
	pssmList = parseDreme.getPSSMListFromDremeFile(dremeFile)
	FindBestScoreSeqForMotif(seqDict, pssmList)
