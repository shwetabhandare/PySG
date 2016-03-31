import sys
import parseKspectrum

def zerolistmaker(n):
	listofzeros = [0] * n;
	return listofzeros;

def createKmerDictList(kmerList):
	kmerDictList = list();
	posDict = dict();

	maxKmerLen = 0;
	for kmer in kmerList:
		if maxKmerLen < len(kmer):
			maxKmerLen = len(kmer);
		posDict = dict();
		for counter, c in enumerate(kmer):
			posDict[counter] = c
		kmerDictList.append(posDict)
	return maxKmerLen, kmerDictList;

def getInitializedLists(numKeys):
	aList = zerolistmaker(numKeys);
	tList = zerolistmaker(numKeys);
	cList = zerolistmaker(numKeys);
	gList = zerolistmaker(numKeys);

	return aList, tList, gList, cList;

def getUpdatedCounts(aList, tList, gList, cList, kmerDictList):
	for kDict in kmerDictList:
		for key, value in kDict.iteritems():
			if value == 'A':
				aList[key] = aList[key] + 1;
			elif value == 'T':
				tList[key] = tList[key] + 1;
			elif value == 'C':
				cList[key] = cList[key] + 1;
			elif value == 'G':
				gList[key] = gList[key] + 1;


if __name__ == "__main__":
	import sys
	kmerFile = sys.argv[1]
	kmerDict = parseKspectrum.FindKspectrumKmers(kmerFile)
	kmerList = kmerDict.keys();
	maxKmerLen, kmerDictList = createKmerDictList(kmerList);
	aList, tList, gList, cList = getInitializedLists(maxKmerLen);
