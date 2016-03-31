import sys
import parseKspectrum
import TAMO
from   TAMO    import MotifTools
from   TAMO.seq import Fasta
import TAMO_Motif

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

def getUpdatedCounts(maxKmerLen, kmerDictList):
	aList, tList, gList, cList = getInitializedLists(maxKmerLen);

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

	print aList
	print cList
	print gList
	print tList

	return aList, tList, gList, cList;

def createNormalizedLists(aList, tList, gList, cList, numKmers):
	aList[:] = [round(float(x)/numKmers, 3) for x in aList]
	tList[:] = [round(float(x)/numKmers, 3) for x in tList]
	gList[:] = [round(float(x)/numKmers, 3) for x in gList]
	cList[:] = [round(float(x)/numKmers, 3) for x in cList]

	print aList
	print cList
	print gList
	print tList

	return aList, tList, gList, cList;

def createPwm(aList, tList, gList, cList):

	pwm = list()
	for i in range(1,len(aList)):
		vals = [float(aList[i]), float(cList[i]), float(gList[i]), float(tList[i])]
		pwm.append(vals)

	print pwm
	return pwm;

def GetKspectrumPWM(kmerFile):
	kmerDict = parseKspectrum.FindKspectrumKmers(kmerFile)
	kmerList = kmerDict.keys();
	maxKmerLen, kmerDictList = createKmerDictList(kmerList);
	aList, tList, gList, cList = getUpdatedCounts(maxKmerLen, kmerDictList);
	numKmers = len(kmerList)

	aList, tList, gList, cList = createNormalizedLists(aList, tList, gList, cList, len(kmerList));
	pwm = createPwm(aList, tList, gList, cList)
	return pwm;

if __name__ == "__main__":
	import sys
	kmerFile = sys.argv[1]
	pwm = GetKspectrumPWM(sys.argv[1])

