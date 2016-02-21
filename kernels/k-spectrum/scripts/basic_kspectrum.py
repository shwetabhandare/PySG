#Input is consolidated file that contains positive and negative data set.
#Run k-spectrum kernel - basic cross validation, and return top k-mers.

import sys
import generate_model
from demo_utils import *
from generate_model import *
import yaml


def BuildModel(trainData, modelFile, numFolds):
	Cs = [ 10**x for x in xrange( -10, 5 ) ]
	param = modelSelection.Param(svm.SVM(), 'C', Cs)
	m = modelSelection.ModelSelector(param, measure ='roc')
	m.train(trainData)
	print "Best C: " + str(m.classifier.C);

	results = m.stratifiedCV(trainData, numFolds=numFolds, seen=1);

	return results, m;


def CrossValidate(trainData, numFolds):
	s = svm.SVM();
	results = s.stratifiedCV(trainData, numFolds=numFolds, seen=1)

	return results;

def ReadConfFile(configFile):
	confMap = readConfFile(configFile);
	return confMap;


if __name__ == '__main__':
	import sys
	confMap = ReadConfFile(sys.argv[1]);

	dataFile = confMap["data"]["dataFile"]
	k1 = int(confMap["kspectrum"]["k1"])
	k2 = int(confMap["kspectrum"]["k2"])
	C = float(confMap["kspectrum"]["C"])
	posLen = int(confMap["data"]["posLen"])
	negLen = int(confMap["data"]["negLen"])
	numFolds = int(confMap["kspectrum"]["folds"]);
	featureFile = confMap["output"]["featureFile"];
	resultFile = confMap["output"]["resultsFile"]
	modelFile = confMap["output"]["modelFile"]

	buildModel = False;

	trainData = generate_model.get_spectrum_data(dataFile, k1, k2, posLen, negLen, normalize=True, repeatCount=True);
	if buildModel:
		results, m = BuildModel(trainData, modelFile, numFolds);
		print_results(results, resultFile, k1, k2)
		find_features(trainData, featureFile, m.classifier.C);
	else:
		results = CrossValidate(trainData, numFolds)
		print_results(results, resultFile, k1, k2)
		find_features(trainData, featureFile);
