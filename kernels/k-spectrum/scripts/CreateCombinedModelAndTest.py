import os
import sys
from PyML.classifiers.svm import SVM
from PyML.containers import ker,labels
from PyML.containers.kernelData import KernelData
from PyML.containers.sequenceData import SequenceData
from PyML import *
import demo_utils


def CreateSparsetDataSet(trainDataFile, testDataFile):

	trainData = SparseDataSet(trainDataFile);
	testData = SparseDataSet(testDataFile);

	return trainData, testData

def BuildModelAndTest(trainData, testData, modelFile, buildModel):

	if (buildModel == 'true'):
		print "* *** BUILDING MODEL *****"
		Cs = [ 10**x for x in xrange( -10, 5 ) ]
		param = modelSelection.Param(svm.SVM(), 'C', Cs)

		m = modelSelection.ModelSelector(param, measure ='roc')
		m.train(trainData)
		results = m.test(testData); #combined HuR, TTP test data
		m.save(modelFile);
		print "**** MODEL OUTPUT **** "
		print m.classifier.C
		print m.log
	else:
		s = SVM();
		s.load(modelFile, trainData);
		results = s.test(testData)

	return results;

def QuickTest(trainData, testData):
	s = svm.SVM(); 
	s.train(trainData);
	results = s.test(testData); #combined HuR, TTP test data
	return results;

def PrintResults(results, resultsFile, featureFile):
	demo_utils.print_results(results, resultsFile);
	print results.getLog();
	print results.getSuccessRate

	#Find features
	demo_utils.find_features(trainData, featureFile);

if __name__ == '__main__':
	import sys

	trainDataFile = sys.argv[1];
	testDataFile  = sys.argv[2];
	featureFile =  sys.argv[3];
	resultsFile = sys.argv[4]
	modelFile = sys.argv[5]
	buildModel = sys.argv[6]


	trainData, testData = CreateSparsetDataSet(trainDataFile, testDataFile);
	results = BuildModelAndTest(trainData, testData, modelFile, buildModel)
	PrintResults(results, resultsFile, featureFile);
