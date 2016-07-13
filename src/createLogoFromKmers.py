import csv
import TAMO
from   TAMO		import MotifTools
from   TAMO.seq import Fasta
import parseKspectrum;

#
#A Content: 27.2844670516
#U Content: 29.8989489198
#G Content: 21.4580725792
#C Content: 21.3585114494
    

if __name__ == "__main__":
	import sys
	kmerFile = sys.argv[1]
	numKmers = int(sys.argv[2])
	logoFileName = sys.argv[3]

	kmerDict = parseKspectrum.FindKspectrumKmers(kmerFile, numKmers)
	kmerList = kmerDict.keys();
	print "Kmer List: ", kmerList
	m_msa = MotifTools.Motif(kmerList, {'A': 0.28, 'C': 0.21, 'G': 0.21, 'T': 0.3})
	print m_msa;
	m_msa._print_counts()
	filename = m_msa.giflogo(logoFileName);

