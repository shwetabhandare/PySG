import fnmatch
import os
import pyngram;
import Distribution_Utils;
import numpy as np;
from scipy import stats;


def GetFastaFiles(dirName, searchRe):
	matches = []
	for root, dirnames, filenames in os.walk(dirName):
	    for filename in fnmatch.filter(filenames, searchRe):
	        matches.append(os.path.join(root, filename))
	return matches

def GetNGramDistributionForFiles(fastaFiles, nLen):
	freqList = []
	probList = []

	for faFile in fastaFiles:
		freq, prob  = Distribution_Utils.GetNgramDistributionForFile(faFile, nLen);
		freqList.append(freq);
		probList.append(prob);

	return freqList, probList;

def CreateNGramList(faList):
	cummulativeNGramList = {};
	for seqDict in faList:
		for key, value in seqDict.iteritems():
			if cummulativeNGramList.has_key(key):
				cummulativeNGramList[key].append(value)
			else:
				cummulativeNGramList[key] = [];
				cummulativeNGramList[key] = [value]

	return cummulativeNGramList;

def ComputeMeanAndStdError(SignalNGramList):
	for key, value in SignalNGramList.iteritems():
		SignalNGramList[key] = [round(np.mean(value), 4), round(stats.sem(value), 4)]

	return SignalNGramList;


if __name__ == "__main__":
	import sys;
	dirName = sys.argv[1]
	nLen = int(sys.argv[2])

	SignalFiles = GetFastaFiles(dirName, "Signal*.fa")
	NoSignalFiles = GetFastaFiles(dirName, "Signal*.fa")

	SignalFreqList, SignalProbList = GetNGramDistributionForFiles(SignalFiles, nLen)
	NoSignalFreqList, NoSignalProbList = GetNGramDistributionForFiles(NoSignalFiles, nLen)

	SignalNGramList = CreateNGramList(SignalProbList);
	NoSignalNGramList = CreateNGramList(NoSignalProbList);

	SignalNGramList = ComputeMeanAndStdError(SignalNGramList);
	NoSignalNGramList = ComputeMeanAndStdError(NoSignalNGramList);

	print SignalNGramList, NoSignalNGramList


