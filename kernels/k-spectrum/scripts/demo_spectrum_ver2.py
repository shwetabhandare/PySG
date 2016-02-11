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

three_utr_combined_file = '/projects/debra/workspace/Data/HuR/Kishore2011-HuR_CLIP-Seq/Combined_Human_Seq/VD_3UTR_combined.txt'
five_utr_combined_file = '/projects/debra/workspace/Data/HuR/Kishore2011-HuR_CLIP-Seq/Combined_Human_Seq/VD_5UTR_combined.txt'
unspliced_utr_combined_file = '/projects/debra/workspace/Data/HuR/Kishore2011-HuR_CLIP-Seq/Combined_Human_Seq/VD_Unspliced_combined.txt'

def fasta_read(file_name) :
	"""read the sequence from a file in fasta format"""
	return [record.sequence for record in fasta.fasta_itr(file_name)]

def get_seq_length_feature(fasta_file):
	seqLengths = [];
	sequences = fasta_read(fasta_file);
	for s in sequences:
		seqLengths.append(len(s));
	return seqLengths;


def get_motif_count(haystack, needle):
	i = 0
	l = list()
	count = 0;
	try:
		while True:
			i = haystack.index(needle, i)
			count = count + 1;
			l.append(i)
			#print i, haystack[i:i+len(needle)];
			i += len(needle)
	except ValueError:
		pass

	return count, l;

def getSequenceFeatures(fasta_file, rna):
	seqLengths = [];
	motifCount = [];
	motifDist = [];
	#nonamer = "TTATTTATT";
	#nonamer = "TATTTAT";
	if rna == 1:
		nonamer = "UAUUUAU";
	else:
		nonamer = "TATTTAT";
	sequences = fasta_read(fasta_file);
	for s in sequences:
		seqLengths.append(len(s));
		count, l = get_motif_count(s, nonamer);
		motifCount.append(count);
		if len(l) == 0:
			motifDist.append(-1);
		elif len(l) == 1:
			motifDist.append(0);
		else:
			#print count;
			d = [];
			distToAdd = 0;
			for x,y in zip(l,l[1:]):
				diff = y - x;
				#print diff;
				d.append(diff);
			distToAdd = min(d);
			#print distToAdd
			motifDist.append(distToAdd);

	return seqLengths, motifCount, motifDist;

## 000 - no features to add.
## 001 - Sequence length
## 010 - Number of nonamers
## 100 - Distance between nonamers

def addSequenceFeatures(fasta_file, total_data, seqFeatureMask, rna):
	if seqFeatureMask == 0:
		return total_data;
	else:
		seqLengths, motifPos, motifDist = getSequenceFeatures(fasta_file, rna);

	print "SeqLen array: " + str(len(seqLengths)) + " Motif Pos array: " + str(len(motifPos)) + " Motif Dist array: " + str(len(motifDist));

	if seqFeatureMask & 0x1:
		print "Adding Seq Len Feature";
		total_data.addFeature(2, seqLengths);
	if seqFeatureMask & 0x2:
		print "Adding Motif Count Feature";
		total_data.addFeature(3, motifPos);
	if seqFeatureMask & 0x4:
		print "Adding Motif Distance Feature";
		total_data.addFeature(4, motifDist);

	return total_data;


def score_features(total_data):
	fs = featsel.FeatureScore ('golub')
	f = featsel.Filter (fs, sigma = 2)
	f.train(total_data)
	#print "Scoring features:";
	#print total_data.numFeatures, total_data.featureID

def get_spectrum_data (fasta_file, k1, k2, pos_data_len, neg_data_len, normalize=False, prefix='None'):
	# sequence data file, K1, K2 where k1 and k2 are k-mer lengths.
	if prefix == 'None':
		total_data = sequenceData.spectrum_data(fasta_file, k1, k2, normalize=normalize);
	else:
		total_data = sequenceData.spectrum_data(fasta_file, k1, k2, normalize=normalize, prefix=prefix);

	print total_data


	#Generate labels
	L = ["+1"] * (pos_data_len);
	L[pos_data_len:neg_data_len] = ["-1"] *(neg_data_len);
	total_data.attachLabels(Labels(L));

	if normalize:
		total_data.normalize();
	return total_data;


def run_svm (data, num_folds=10, C=1, seed=1):
	#s.train(total_data)
	r = SVM(C=C).stratifiedCV(data, numFolds=num_folds,  seen=seed);

	return r

def get_incorrectly_predicted_ids(pid, l):
	ids = []
	for i,j in itertools.izip(pid, l):
		p = [];
		q = []
		for x,y in itertools.izip(i,j):
			if y == 0:
				p.append(x)
			else:
				q.append(x);
		ids.append(p);
	return ids;

def convertIdsToTxIds(errorIds):
	"""read the sequence from a file in fasta format"""
	txList = list()
	for record in fasta.fasta_itr(three_utr_combined_file):
		header = record.header
		geneId, transcriptId, geneName = header.split('|')
		txList.append(transcriptId);

	for list1 in itertools.izip(errorIds):
		print list1
		for list2 in itertools.izip(list1[0]):
			idx = list2[0]
			print txList[idx]
	return txList;


def print_results(r, res_file):
	res_wr = csv.writer(open(res_file, 'wb'))
	res_wr.writerow(['Success Rate','Balanced Success Rate','Sensitivity','PPV','ROC'])
	res_wr.writerow([r.getSuccessRate(),r.balancedSuccessRate,r.getSensitivity(),r.ppv,r.getROC()])

	print 'Success Rate:  %5.3f.' % r.getSuccessRate()
	print 'Balanced Success Rate:  %5.3f.' % r.getSuccessRate()
	print 'Sensitivity:  %5.3f.' % r.getSensitivity()
	print 'PPV:  %5.3f.' % r.ppv
	print 'ROC:  %5.3f.' % r.getROC()

	## Get the incorrectly predicted labels.
	pid = r.getPatternID()
	pl = r.getPredictedClass();
	gl = r.getGivenClass();
	print "Results Log: ";
	print r.getLog();
	print "Predicted class length: %d, Given Class Length: %d"%(len(pl), len(gl));
	l = compare_predicted_given_labels(r, pl,gl)
	ids = get_incorrectly_predicted_ids(pid, l);
	print "Incorrectly predicted labels: "
	print ids
	#convertIdsToTxIds(ids);

def compare_predicted_given_labels(r, pl, gl):
	l = []
	for i, j in itertools.izip(pl,gl):
		p = []
		for x,y in itertools.izip(i,j):
			# x, y are given, and actual class labels.
			if (x==y):
				p.append(1);
			else:
				p.append(0);
		l.append(p);
	return l


## Program starts here.
numArgs = len(sys.argv);
print numArgs
if numArgs < 8:
	print "USAGE: python demo_spectrum_ver2.py <filename> <k1> <k2> <posLen> <negLen> <seqFeatureMask> <rna|dna>";	 
	print "<k1>: Size of k-mer, <k2> can be zero.";
	print "<posLen>: The size of the positive data set";
	print "<negLen>: The size of the negative data set";
	print "<numFolds>: Number of folds for cross validation. If 0, specifiy train in arg 1, and test in last arg."
	print "<seqFeatureMask>: 0 - no seq features, 0x1 - seq length, 0x2 - number of motifs, 0x4 - distance between motifs. \nYou could get all features if you did 111";
	print "<rna|dna>: 1 - RNA, 0 - DNA";
	print "<featurevector>: If no file specified, FeatureVectors.txt will be used."
	print "<DataFile> : Specify the name of the file to save the spectrum data. Default: kspectrum.data";

	exit(1);

three_utr_combined_file = sys.argv[1];
k1 = int(sys.argv[2])
k2 = int(sys.argv[3])
posLen = int(sys.argv[4])
negLen = int(sys.argv[5])
numFolds = int(sys.argv[6]);
seqMask = int(sys.argv[7]);
rna = int(sys.argv[8]);
if numArgs > 10:
	featureVectorFile=sys.argv[9]
else:
	featureVectorFile = "FeatureVectors.txt"; 

if numArgs == 11:
	dataFile=sys.argv[10]
else:
	dataFile = "spectrum.data";

print "=== Three UTR, k-mer  length: 5 ===\n"
three_utr_data = get_spectrum_data(three_utr_combined_file, k1, k2, posLen, negLen, True);
three_utr_data.save(dataFile);
#three_utr_data = addSequenceFeatures(three_utr_combined_file, three_utr_data, seqMask, rna);
r1 = run_svm(three_utr_data, numFolds, 1, 1)
print_results(r1, 'three_utr_5_spectrum.csv');
demo_utils.find_features(three_utr_data, featureVectorFile);
score_features(three_utr_data);

#print "=== Three UTR, k-mer  length: 7 ===\n"
#three_utr_data = get_spectrum_data(three_utr_combined_file, 7, None, 275, 286, True);
#print three_utr_data
#r2 = run_svm(three_utr_data, 5, 1, 1)
#print_results(r2, "three_utr_7_spectrum.csv");
