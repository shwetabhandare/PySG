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
from PyML.classifiers.composite import FeatureSelect


from PyML.evaluators import resultsObjects
from PyML.evaluators import roc as roc_mod

from PyML import *
import demo_utils
from demo_utils import *
import generate_model
import yaml

configFile = sys.argv[1]
confMap = readConfFile(configFile);


## Program starts here.
trainSeqFile = confMap["data"]["dataFile"]
trainPosLen = int(confMap["data"]["posLen"])
trainNegLen = int(confMap["data"]["negLen"])
numFolds = int(confMap["kspectrum"]["folds"]);
modelFile = confMap["output"]["modelFile"];
dataFile = confMap["output"]["trainDataFile"]
resultFile = confMap["output"]["resultsFile"]
paramsFile = confMap["output"]["paramsFile"]
featureFile = confMap["output"]["featureFile"]


Cs = [ 10**x for x in xrange( -10, 5 ) ]
# trainSeqFile = sys.argv[1];
# trainLen = int(sys.argv[2]); #equal pos/neg len
# modelFile = sys.argv[3];
# dataFile = sys.argv[4];
# bestModelFile = sys.argv[5];
# outputParamsFile = sys.argv[6]
bestCVResults = None
bestROC = -numpy.Inf
bestROC50 = -numpy.Inf
bestSuccessRate = None
bestW = None;
bestK1 = 0;
bestK2 = 0;
K1 = [8, 9, 10, 11, 12, 13]
K2 = [8, 9, 10, 11, 12, 13]


for k1 in K1:
   for k2 in [x for x in K2 if x >= k1]:
      for C in Cs:
         print "**** Train/Test with K1: " + str(k1) + ", k2: " + str(k2) +", C: " + str(C);
         trainData = generate_model.get_spectrum_data(trainSeqFile, k1, k2, trainPosLen, trainNegLen, normalize=True, repeatCount=True);
         s = svm.SVM(C=C); 
         s.train(trainData)
         results = s.stratifiedCV(trainData, numFolds=numFolds, seen=1);
         labels = results.getGivenClass();

         roc  = results.getROC();
         roc50 = results.getROCn();
         successRate = results.getSuccessRate();
         print "ROC: " + str(roc)  + ", ROC50: " + str(roc50) + ", Success Rate: " + str(successRate);
         if roc > bestROC: 
            # Delete the old data file.
            os.remove(dataFile) if os.path.exists(dataFile) else None
            # Save the training data for K1, K2 combination.
            trainData.save(dataFile);

            # delete the old model file
            os.remove(modelFile) if os.path.exists(modelFile) else None
            # Save the model file
            s.save(modelFile);

            # Updateww the best numbers


            bestROC50 = roc50;
            bestSuccessRate = successRate;
            bestROC = roc;
            bestC = C; 
            bestCVResults = results;
            bestK1 = k1;
            bestK2 = k2;

            # Save results, and parameters in the paramsFile.
            print_results(results, paramsFile, k1, k2, bestC);
            find_features(trainData, featureFile)

				# ofile = open(  outputParamsFile,  'w'  );
				# ofile.write( "roc50: " + str( roc50 ) + "\n" )
				# ofile.write( "roc: " + str( roc ) + "\n" )
				# ofile.write( "bestSuccessRate: " + str( bestSuccessRate ) + "\n" )
				# ofile.write( "K1: " + str( k1 ) + ", K2: " + str(k2) + "\n" )
				# ofile.write( "bestC: " + str( bestC ) + "\n" )
				# ofile.write( "Features: " + str( len(trainData.featureID) ) + "\n" )
				# ofile.close();



print ("roc: %s, roc50: %s, SuccessRate: %s, bestC: %s"%(str(roc), str(roc50), str(successRate), str(bestC)));
print_results(bestCVResults, resultFile, bestK1, bestK2, bestC);
