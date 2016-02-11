import os
import numpy
import sys
import itertools
import csv

from PyML.classifiers import svm,multi,ridgeRegression,knn,composite,modelSelection
from PyML.classifiers.svm import SVM
from PyML.feature_selection import featsel
from PyML.containers import ker,labels
from PyML.containers import vectorDatasets
from PyML.containers.aggregate import Aggregate
from PyML.containers.kernelData import KernelData
from PyML.containers.sequenceData import SequenceData
from PyML.classifiers import platt
from PyML.preproc import preproc
from PyML.utils import fasta


from PyML.evaluators import resultsObjects
from PyML.evaluators import roc as roc_mod

from PyML import *
import demo_utils
import generate_model

import yaml
import string

# Read train data from the best model.
# Read model file for best parameters.
# Read test data - create spectrum data with k1, k2 from the model.
# 

def readResult(results_yml):
	with open(results_yml) as f:
		return yaml.load(f)

def createSpectrumData(testSeqFile, resultsConf):

	filename_list = string.split(testSeqFile, "_")
	list_len = len(filename_list)


	testPosLen = int(filename_list[list_len - 3])
	testNegLen = int(filename_list[list_len - 2])

	k1 = int(resultsConf["results"]["k1"])
	k2 = int(resultsConf["results"]["k2"])

	testSpectrumData = demo_utils.get_spectrum_data(testSeqFile, k1, k2, 
		testPosLen, testNegLen, True);

	return testSpectrumData;
	

def inDomainTest(trainModelFile, trainSpectrumData, testSpectrumData, resultsFile):
	new_svm = SVM()
	new_svm.load(trainModelFile, trainSpectrumData)
	results = new_svm.test(testSpectrumData)

	demo_utils.print_results(results, resultsFile)

def getFeatures(trainSpectrumData, testSpectrumData, resultsConf, testFeaturesFile):
	trainC = float(resultsConf["results"]["C"])
	s = svm.SVM(C=trainC)
	s.train(trainSpectrumData);
	r = s.test(testSpectrumData);
	features = demo_utils.get_features(s.model.w, trainSpectrumData);
	demo_utils.writeFeaturesToFile(features, testFeaturesFile)
	

if __name__ == '__main__':
	import sys
	resultsYml = sys.argv[1]
	testSeqFile = sys.argv[2]
	trainModelFile = sys.argv[3]
	trainSeqFile = sys.argv[4]
	resultsFile = sys.argv[5]
	testFeaturesFile = sys.argv[6]
	testSpectrumFile = sys.argv[7]

	resultsConf = readResult(resultsYml);
	print resultsConf;

	trainSpectrumData = createSpectrumData(trainSeqFile, resultsConf);
	testSpectrumData = createSpectrumData(testSeqFile, resultsConf);
	testSpectrumData.save(testSpectrumFile);
	inDomainTest(trainModelFile, trainSpectrumData, testSpectrumData, resultsFile)
	getFeatures(trainSpectrumData, testSpectrumData, resultsConf, testFeaturesFile)
