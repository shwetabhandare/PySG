#import re, ahocorasick,random,time
from __future__ import division
import parseDreme
import parseRealKmers
import parseKspectrum
import fasta

import splitKmerInDict;
import compareKmerCommon;


#from itertools import imap

def ComparePWMKmers(realKmerDict, predictedKspectrumFile, posFile, negFile, pssmList):
	numFP = numTP = numFN = 0	


	pwm = splitKmerInDict.GetKspectrumPWM(predictedKspectrumFile);
	if len(pwm) == 0:
		numFP = compareKmerCommon.getTotalFNKmerNotFound(realKmerDict);		
	else:
		numTP, numFP, numFN = compareKmerCommon.GetTotalNumbers(realKmerDict, posFile, negFile, None, pwm, None)

	return numTP, numFP, numFN;

def GetDremePSSM(predictedDremeFile):
	pssmList = parseDreme.getPSSMListFromDremeFile(predictedDremeFile);
	return pssmList;

def ComparePSSMKmers(realKmerDict, predictedDremeFile, posFile, negFile, pssmList):
	numTP = numFP = numFN = 0;
	if len(pssmList) == 0:
		numFP = compareKmerCommon.getTotalFNKmerNotFound(realKmerDict);		
	else:
		numTP, numFP, numFN = compareKmerCommon.GetTotalNumbers(realKmerDict, posFile, negFile, pssmList, None, None)

	return numTP, numFP, numFN;

def CompareTextMotifKmers(realKmerDict, predictedKmers, posFile, negFile):
	numFP = numTP = numFN = 0	
	if len(predictedKmers) == 0:
		numFN = compareKmerCommon.getTotalFNKmerNotFound(realKmerDict);
	else:
		numTP, numFP, numFN = compareKmerCommon.GetTotalNumbers(realKmerDict,  posFile, negFile, None, None, predictedKmers)
	return numTP, numFP, numFN;

def CompareKspectrumPredictedKmers(realCsvFile, predictedKspectrumFile, posFile, negFile, numKmers):
	realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);
	predictedKmerDict	= parseKspectrum.FindKspectrumKmers(predictedKspectrumFile, numKmers);
	numTP = 0
	numFP = 0
	numFN = 0

	if len(predictedKmerDict) == 0:
		numFP = compareKmerCommon.getTotalFNKmerNotFound(realKmerDict);		
	else:
		numTP, numFP, numFN = compareKmerCommon.GetTotalNumbers(realKmerDict, posFile, negFile, None, None, None, predictedKmerDict)

	return numTP, numFP, numFN

def CompareKspectrumKmers(realCsvFile, predictedKspectrumFile, posFile, negFile):

	realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);
	numTP = 0
	numFP = 0
	numFN = 0

	pwm	= splitKmerInDict.GetKspectrumPWM(predictedKspectrumFile);
	numTP, numFP, numFN = ComparePWMKmers(realKmerDict, predictedKspectrumFile, posFile, negFile, pwm);
	
	return numTP, numFP, numFN;

def CompareDremeKmers(realCsvFile, predictedDremeFile, posFile, negFile, textMotif):

	realKmerDict = parseRealKmers.GetRealKmerDict(realCsvFile);

	numTP = 0
	numFP = 0
	numFN = 0

	if textMotif == 1:
		predictedKmers = parseDreme.getPredictedDremeMotifs(predictedDremeFile);
		numTP, numFP, numFN = CompareTextMotifKmers(realKmerDict, predictedKmers, posFile, negFile)
	else:
		pssmList = GetDremePSSM(predictedDremeFile);
		numTP, numFP, numFN = ComparePSSMKmers(realKmerDict, predictedDremeFile, posFile, negFile, pssmList)
	
	return numTP, numFP, numFN;

def GetAccuracy(numTP, numTN, total):
	if ((numTP + numTN) == 0):
		accuracy = 0;
	else:
		accuracy = (numTP + numTN) / total;

	return accuracy;

def GetSpecificity(numFP, totalNeg):
	if numFP == 0:
		specificity = 0;
	else:
		specificity = numFP / totalNeg;

	return specificity;

def GetSensitivityAndPPV(numTP, numFP, numFN):
	if ((numTP + numFN) == 0):
		sensitivity = 0;
	else:
		sensitivity = numTP / (numTP + numFN)
	
	if ((numTP + numFP) == 0):
		ppv = 0;
	else:
		ppv =  numTP / (numTP + numFP)

	return sensitivity, ppv;

if __name__ == "__main__":
	import sys
	realCsvFile = sys.argv[1]
	predictedDremeFile= sys.argv[2]
	posFile = sys.argv[3]
	negFile = sys.argv[4]
	numKmers = int(sys.argv[5])
	#textMotif = int(sys.argv[5])

	numTP, numFP, numFN = CompareKspectrumPredictedKmers(realCsvFile, predictedDremeFile, posFile, negFile, numKmers)
	print "Num TP: ", str(numTP), ", Num FP: ", str(numFP), ", Num FN: ", str(numFN)
	sensitivity, ppv = GetSensitivityAndPPV(numTP, numFP, numFN);
	print sensitivity, ppv
