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
if numArgs < 11:
	print "USAGE: python demo_test.py <k1> <k2> <trainFile> <testFilee>  \
	      <trainPosLen> <trainNegLen> <testPosLen> <testNegLen> <trainFeatureFile> <testFeatureFile>";
	
	exit(1);
	
Cs = [ 10**x for x in xrange( -3, 4 ) ]
k1 = int(sys.argv[1])
k2 = int(sys.argv[2])
trainFile = sys.argv[3];
testFile = sys.argv[4];

trainPosLen = int(sys.argv[5])
trainNegLen = int(sys.argv[6])

testPosLen = int(sys.argv[7])
testNegLen = int(sys.argv[8])

trainFeatureData = sys.argv[9];
testFeatureData = sys.argv[10];
	
trainData = demo_utils.get_spectrum_data(trainFile, k1, k2, trainPosLen, trainNegLen, True);
trainData.save(featureData);
for C in Cs:

   print C;
   #print "Train for C : " + str(C);
   s = svm.SVM(C=C);
   s.train(trainData);
   testData = demo_utils.get_spectrum_data(testFile, k1, k2, testPosLen, testNegLen, True);
   results = s.test(testData);

exit(1);

for C in Cs:

   print "Train for C : " + C;
   s = svm.SVM(C=C);
   s.train(trainData);
   fileBaseName = os.path.basename(trainFile);
   modelFileName = os.path.splitext(os.path.basename(trainFile))[0] + "_model_" + C + ".pyml";
   s.save(modelFileName);

   testData = demo_utils.get_spectrum_data(testFile, k1, k2, testPosLen, testNegLen, True);
   #s2 = SVM();
   #s2.load(modelFileName, trainData);
   results = s.test(testData);
   print results;
