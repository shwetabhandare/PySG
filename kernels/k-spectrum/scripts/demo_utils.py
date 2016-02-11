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
import yaml;

from PyML.evaluators import resultsObjects

from PyML import *

def readConfFile(configFile):
	f = open(configFile)
	confMap = {}
	confMap = yaml.load(f)
	f.close()
	return confMap;
	
def print_results(r, target=None, k1=1, k2=1, C=1):
	successRate = r.getSuccessRate()
	balancedSuccessRate = r.getBalancedSuccessRate()
	sensitivity = r.getSensitivity();
	print successRate;
	if target != None:
		with open(target, 'w') as h:
			data = dict(
		 		results = dict (
			  		SuccessRate = "{0:5.3f}".format(successRate),
			  		BalancedSuccessRate = "{0:5.3f}".format(balancedSuccessRate),
			  		Sensitivity = "{0:5.3f}".format(sensitivity),
			  		PPV = "{0:5.3f}".format(r.ppv),
			  		ROC = "{0:5.3f}".format(r.getROC()),
			  		k1 = k1,
			  		k2 = k2,
			  		C = C,
			))
			print "Data: ", data
			print yaml.dump(data)
			with open(target, 'w') as outfile:
		 		outfile.write( yaml.dump(data ))
	else:
		print 'Success Rate:  %5.3f.' % r.getSuccessRate()
		print 'Balanced Success Rate:  %5.3f.' % r.getBalancedSuccessRate()
		print 'Sensitivity:  %5.3f.' % r.getSensitivity()
		print 'PPV:  %5.3f.' % r.ppv
		print 'ROC:  %5.3f.' % r.getROC()

def run_svm (data, num_folds=10, C=1, seed=1):
	#s.train(total_data)
	r = SVM(C=C).stratifiedCV(data, numFolds=num_folds,  seen=seed);

	return r

def get_spectrum_data (fasta_file, k1, k2, pos_data_len, neg_data_len, normalize=True, prefix='None'):
	# sequence data file, K1, K2 where k1 and k2 are k-mer lengths.
	if prefix == 'None':
		total_data = sequenceData.spectrum_data(fasta_file, k1, k2);
	else:
		total_data = sequenceData.spectrum_data(fasta_file, k1, k2, prefix=prefix);


	#Generate labels
	L = ["+1"] * (pos_data_len);
	L[pos_data_len:neg_data_len] = ["-1"] *(neg_data_len);
	total_data.attachLabels(Labels(L));

	if normalize:
		total_data.normalize();

	print total_data
	return total_data;

def RFE_FeatureSelect(trainData):
	rfe = featsel.RFE()
	rfe.select(trainData)
	print trainData.featureID

def get_features(w, total_data):
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
	return items;

def writeFeaturesToFile(items, featureVectorFile):

	file = open(featureVectorFile, "w");
	for item in items:
		#print item[0], item[1]
		writeStr = str(item[0]) + "," + str(item[1]) + "\n";
		file.write(writeStr);
		#file.write(item);
		#print item
	file.close();

def find_features(total_data, featureVectorFile, C=10):
	#rfe = featsel.RFE()
	#rfe.select(total_data)
	#print total_data;
	s = svm.SVM()
	s.C = C;
	s.train(total_data)
	w = s.model.w

	items = get_features(w, total_data);
	writeFeaturesToFile(items, featureVectorFile);



