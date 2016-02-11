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

def find_features(w, total_data, featureVectorFile):
	print "Finding features:"
	count = 0
	kmer_weight = dict()
	for val in w:
		featureID = total_data.featureID[count]
		#print featureID, val
		count = count + 1;
		kmer_weight[featureID] = val

	items = [(v, k) for k, v in kmer_weight.items()]
	items.sort()
	items.reverse()

	file = open(featureVectorFile, "w");
	for item in items:
		#print item[0], item[1]
		writeStr = str(item[0]) + "," + str(item[1]) + "\n";
		file.write(writeStr);
		#file.write(item);
		#print item
	file.close();

## Program starts here.
numArgs = len(sys.argv);
print numArgs
if numArgs < 4:
	print "USAGE: python demo_test.py <train seq file> <test seq file> <lenTrain> <lenTest>  <outputFeatureFile>"
	exit(1);

#Add input train/test file, and posLen, negLen.
#call get_spectrum_file.
Cs = [ 10**x for x in xrange( -10, 9 ) ]
#Cs = [ 10**x for x in xrange( -2, 2 ) ]

trainSeqFile = sys.argv[1];
testSeqFile = sys.argv[2];
trainLen = int(sys.argv[3]); #equal pos/neg len
testLen = int(sys.argv[4]);
outFeatureFile = sys.argv[5];


bestC = None
bestAUC = -numpy.Inf
bestFP = None
bestTP = None
K1 = [7, 8, 9, 10, 11, 12, 13]
K2 = [7, 8, 9, 10, 11, 12, 13]

result_file = open("K-spectrum.txt", 'w');
for k1 in K1:
   for k2 in K2:
      for C in Cs:
         print "**** Train/Test with K1: " + str(k1) + ", k2: " + str(k2) +", C: " + str(C);
         trainData = generate_model.get_spectrum_data(trainSeqFile, k1, k2, trainLen, trainLen, True);
         folds = [];
         s = svm.SVM(C=C);
         s.train(trainData);

         #testData = SparseDataSet(testFeatureFile);
         testData = demo_utils.get_spectrum_data(testSeqFile, k1, k2, testLen, testLen, True);

         results = s.test(testData);
         labels = results.getGivenClass();
         dvals = results.getDecisionFunction();
         folds.append( (dvals, labels) )

         demo_utils.print_results(results);
         print "Results Log: ";
         results.getLog();
         fpc, tpc, area = roc_mod.roc_VA(folds, None); 
         print "Area: " + str(area);
         if area > bestAUC:
            bestAUC = area;
            bestFP = fpc;
            bestTP = tpc;
            bestC = C; 
            ofile = open(  'roc%s.txt' % (str(C)) , 'w'  );
            ofile.write( "area: " + str( area ) + "\n" )
            ofile.write( "bestFP: " + str( bestFP ) + "\n" )
            ofile.write( "bestTP: " + str( bestTP ) + "\n" )
            ofile.write( "bestC: " + str( bestC ) + "\n" )
            ofile.close();
            print ("area: %s, bestFP: %s, bestTP: %s, bestC: %s"%(str(area), str(bestFP), str(bestTP), str(bestC)));

         result_file.write("%s, %s, %s, %s, %s, %s, %s\n"%(str(k1), str(k2), str(results.getSuccessRate()), str(results.getROC()), str(bestAUC), str(C), str(bestC)));
result_file.close();
exit(1);

find_features(s.model.w, trainData, outFeatureFile);
