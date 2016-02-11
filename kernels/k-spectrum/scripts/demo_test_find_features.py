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

from PyML import *
import demo_utils



## Program starts here.
numArgs = len(sys.argv);
print numArgs
if numArgs < 4:
	print "USAGE: python demo_test.py <trainFeatureFile> <testFeatureFile> <outputFeatureFile>"
	exit(1);

trainFeatureFile = sys.argv[1];
testFeatureFile = sys.argv[2];
outFeatureFile = sys.argv[3];

trainData = SparseDataSet(trainFeatureFile);
testData = SparseDataSet(testFeatureFile);

rfe = featsel.RFE()
rfe.select(trainData);

s = svm.SVM();
s.train(trainData);
s.test(testData);

