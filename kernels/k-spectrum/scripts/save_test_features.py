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
k1 = int(sys.argv[1])
k2 = int(sys.argv[2])
testFile = sys.argv[3];
testPosLen = int(sys.argv[4])
testNegLen = int(sys.argv[5])
featureData = sys.argv[6];

testData = demo_utils.get_spectrum_data(testFile, k1, k2, testPosLen, testNegLen, True);
testData.save(featureData);
