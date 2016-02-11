#This script creates k-spectrum data.
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




## Program starts here.
numArgs = len(sys.argv);
print numArgs
if numArgs < 4:
	print "USAGE: python createTestSpectrumData.py <test seq file> <lenTest> <k1> <k2> <testDataFile>"
	exit(1);

testSeqFile = sys.argv[1];
testLen = int(sys.argv[2]); #equal pos/neg len
k1 = int(sys.argv[3]);
k2 = int(sys.argv[4]);
testDatFile = sys.argv[5];

testData = generate_model.get_spectrum_data(testSeqFile, k1, k2, testLen, testLen, True);
testData.save(testDatFile);

if __name__ == '__main__':
	import sys
	combinedFile = sys.argv[1]
	print "Combined File: ", combinedFile;
	prefix = sys.argv[2]
	confFile = sys.argv[3]
	CreateConf(combinedFile, prefix, confFile)
