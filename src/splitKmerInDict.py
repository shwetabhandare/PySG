import sys
import parseKspectrum

kmerFile = sys.argv[1]
kmerDict = parseKspectrum.FindKspectrumKmers(kmerFile)
kmerList = kmerDict.keys();
kmerDictList = list();
posDict = dict();

for kmer in kmerList:
	print kmer
	for counter, c in enumerate(kmer):
		posDict[counter] = c
	kmerDictList.append(posDict)

for kDict in kmerDictList:
	print kDict

