import sys
import os
import numpy
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

from PyML.evaluators import resultsObjects
from PyML.utils import misc

from PyML import *

def read_kernel_matrix(kernelInFile):
	kdata = KernelData(kernelInFile);
	l = Labels(kernelInFile)
	kdata.attachLabels(l)
	#print kdata;
	return kdata;

def run_svm (trainDataFile, trainData, testData):
	s = SVM();
	#param = modelSelection.Param(s, 'C', [0.0001,0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000])
	#m = modelSelection.ModelSelector(param, measure='roc');
	#m.train(trainData);
	#m.save("trainModel.pyml");
	s.train(trainData);
	s.save("trainModel.pyml");
	new_svm = SVM();
	new_svm.load("trainModel.pyml", trainData);
	r = new_svm.test(testData);
	return r

def print_results(r, res_file):
	res_wr = csv.writer(open(res_file, 'wb'))
	res_wr.writerow(['Success Rate','Balanced Success Rate','Sensitivity','PPV','ROC'])
	res_wr.writerow([r.getSuccessRate(),r.balancedSuccessRate,r.getSensitivity(),r.ppv,r.getROC()])

	print 'Success Rate:  %5.3f.' % r.getSuccessRate()
	print 'Balanced Success Rate:  %5.3f.' % r.getSuccessRate()
	print 'Sensitivity:  %5.3f.' % r.getSensitivity()
	print 'PPV:  %5.3f.' % r.ppv
	print 'ROC:  %5.3f.' % r.getROC()

	roc_file = res_file + '.eps'
	print roc_file
	r.plotROC(roc_file, show=False);



trainDataFile = sys.argv[1]
testDataFile = sys.argv[2]

trainData = read_kernel_matrix(trainDataFile)
testData = read_kernel_matrix(testDataFile)

r = run_svm(trainDataFile, trainData, testDataFile);
#print_results(r, 'results.csv')
