import fnmatch
import os
import pyngram;
import Distribution_Utils;
import numpy as np;
from scipy import stats;
import generateGraphs;


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


def ComputeNgramDistribution(dirName, nLen, graphName):
	#three_utr_freq, three_utr_prob = Distribution_Utils.Compute3UtrNgramDistibution(nLen);

	three_utr_prob = {'AA': 0.0878, 'AC': 0.0478, 'GT': 0.0542, 'AG': 0.0674, 'CC': 0.0604, 'TT': 0.1023, 'CG': 0.0121, 'GG': 0.0567, 'GC': 0.0468, 'AT': 0.0693, 'GA': 0.0568, 'TG': 0.0781, 'CT': 0.073, 'CA': 0.0678, 'TC': 0.0583, 'TA': 0.0602}
	SignalFiles = GetFastaFiles(dirName, "Signal*.fa")
	NoSignalFiles = GetFastaFiles(dirName, "NoSignal*.fa")

	SignalFreqList, SignalProbList = GetNGramDistributionForFiles(SignalFiles, nLen)
	NoSignalFreqList, NoSignalProbList = GetNGramDistributionForFiles(NoSignalFiles, nLen)

	SignalNGramList = CreateNGramList(SignalProbList);
	NoSignalNGramList = CreateNGramList(NoSignalProbList);

	SignalNGramDict = ComputeMeanAndStdError(SignalNGramList);
	NoSignalNGramDict = ComputeMeanAndStdError(NoSignalNGramList);

	graphTitle = str(nLen) + "-gram Distribution for 3'UTR and Generated Fasta files"
	graphName = dirName + "/" + graphName;
	generateGraphs.PlotNGrams(SignalNGramDict, NoSignalNGramDict, three_utr_prob, graphTitle, graphName);


if __name__ == "__main__":
	import sys;
	dirName = sys.argv[1]
	nLen = int(sys.argv[2])
	graphName = sys.argv[3]

	ComputeNgramDistribution(dirName, nLen, graphName)

