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

from PyML.evaluators import resultsObjects
from PyML.classifiers.composite import FeatureSelect
from PyML.utils import fasta

from PyML import *
#import createFeatures

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

def getSequenceFeatures(fasta_file):
	seqLengths = [];
	motifCount = [];
	motifDist = [];
	#nonamer = "TTATTTATT";
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
			print count;
			d = [];
			distToAdd = 0;
			for x,y in zip(l,l[1:]):
				diff = y - x;
				print diff;
				d.append(diff);
			distToAdd = min(d);
			print distToAdd
			motifDist.append(distToAdd);

	print len(motifDist), len(motifCount);
	return seqLengths, motifCount, motifDist;

def get_spectrum_data (fasta_file, k1, k2, pos_data_len, neg_data_len, normalize=False, prefix='None'):
	# sequence data file, K1, K2 where k1 and k2 are k-mer lengths.
	seqLengths = [];
	if prefix == 'None':
		if k2 == 'None':
			total_data = sequenceData.spectrum_data(fasta_file, k1, normalize=normalize);
		else:
			total_data = sequenceData.spectrum_data(fasta_file, k1, k2, normalize=normalize);
	else:
		if k2 == 'None':
			total_data = sequenceData.spectrum_data(fasta_file, k1, normalize=normalize, prefix=prefix);
		else:
			total_data = sequenceData.spectrum_data(fasta_file, k1, k2, normalize=normalize, prefix=prefix);

	seqLengths, motifPos, motifDist = getSequenceFeatures(fasta_file);
	
	total_data.addFeature(2, seqLengths);
	total_data.addFeature(3, motifPos);

	#print motifPos
	#print motifDist
	total_data.addFeature(4, motifDist);
	#Generate labels
	L = ["+1"] * (pos_data_len);

	L[pos_data_len:neg_data_len] = ["-1"] *(neg_data_len);

	total_data.attachLabels(Labels(L));


	if normalize:
		total_data.normalize();
	return total_data;


def run_svm (data, num_folds=5, C=1, seed=1):
	r = SVM(C=C).stratifiedCV(data, numFolds=num_folds,  seen=seed);
	return r;


def find_features(total_data):
	rfe = featsel.RFE()
	rfe.select(total_data)
	s = svm.SVM();
	s.train(total_data)
	w = s.model.w

	print "Weight vector"
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

	for item in items:
		print item
	

def score_features(total_data):
	fs = featsel.FeatureScore ('golub')
	f = featsel.Filter (fs, sigma = 2)
	f.train(total_data)
	print total_data.numFeatures, total_data.featureID

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

def print_results(r, res_file):
	res_wr = csv.writer(open(res_file, 'wb'))
	res_wr.writerow(['Success Rate','Balanced Success Rate','Sensitivity','PPV','ROC', 'ROC50'])
	res_wr.writerow([r.getSuccessRate(),r.balancedSuccessRate,r.getSensitivity(),r.ppv,r.getROC(), r.roc50]) 

	print 'Success Rate:  %5.3f.' % r.getSuccessRate()
	print 'Balanced Success Rate:  %5.3f.' % r.getSuccessRate()
	print 'Sensitivity:  %5.3f.' % r.getSensitivity()
	print 'PPV:  %5.3f.' % r.ppv
	print 'ROC:  %5.3f.' % r.getROC()
	print 'ROC50:  %5.3f.' % r.roc50

	## Get the incorrectly predicted labels.
	#pid = r.getPatternID()
	#pl = r.getPredictedClass();
	#gl = r.getGivenClass();
	#print r.getLog();
	#l = compare_predicted_given_labels(pl,gl)
	#ids = get_incorrectly_predicted_ids(pid, l);
	#print "Incorrectly predicted labels: "
	#print ids

	#Let's plot the ROC 
	roc_file = res_file + '.pdf'
	roc_file = res_file + '.eps'
	print roc_file
	#r.plotROC(roc_file, show=False);

#three_utr_homo_ttp_file = '/projects/debra/workspace/Data/TTP/Homo/DendriticPaper/ttpHomoComplete3UTR.txt'

#three_utr_combined_file = sys.argv[1]
#five_utr_combined_file = sys.argv[2]
#unspliced_utr_combined_file = '/projects/debra/workspace/Data/HuR/Kishore2011-HuR_CLIP-Seq/Combined_Human_Seq/VD_Unspliced_combined.txt'

def compare_predicted_given_labels(pl, gl):
	l = []
	for i, j in itertools.izip(pl,gl):
		p = []
		for x,y in itertools.izip(i,j):
			if (x==y):
				p.append(1);
			else:
				p.append(0);
		l.append(p);
	return l

def run_mixed_spectrum_kernel(k1, k2):

	#print " == Mixed Spectrum Kernel : Three UTR, %d, %d == \n"%(k1,k2)
	#three_utr_mixed_data = get_spectrum_data(three_utr_combined_file, k1, k2, 275, 286, True)
	#print three_utr_mixed_data
	#r1 = run_svm(three_utr_mixed_data, 5, 1, 1)
	#print_results(r1, 'mixed_three_utr_spectrum.csv');
	#find_features(three_utr_mixed_data);

	#print " == Mixed Spectrum Kernel : Five UTR, 5, 7 == \n"
	#five_utr_mixed_data = get_spectrum_data(five_utr_combined_file, k1, k2, 273, 286, True)
	#print five_utr_mixed_data
	#r2 = run_svm(five_utr_mixed_data, 5, 1, 1)
	#print_results(r2, 'mixed_five_utr_5_7_spectrum.csv');
	#find_features(five_utr_mixed_data);

	print " == Mixed Spectrum Kernel : Unspliced UTR, 5, 7 == \n"
	unspliced_utr_mixed_data = get_spectrum_data(unspliced_utr_combined_file, k1, k2, 279, 288, True);
	print unspliced_utr_mixed_data

	r3 = run_svm(unspliced_utr_mixed_data, 5, 1, 1)
	print_results(r3, 'mixed_five_utr_5_7_spectrum.csv');
	find_features(unspliced_utr_mixed_data);

def run_three_utr_spectrum_kernel(three_utr_homo_ttp_file, k1, k2, posLen, negLen):
	## HuR binders - 3utr : pos: 275, neg: 286
	if k2 == 0:
		k2 = None
	three_utr_five_spectrum_data = get_spectrum_data(three_utr_homo_ttp_file, k1, k2, posLen, negLen, True)
	print three_utr_five_spectrum_data

	r1 = run_svm(three_utr_five_spectrum_data, 5, 1, 1)
	print " ===== RESULTS ======"
	print_results(r1, 'three_utr_5_spectrum.csv');
	print " ===== SCORING FEATURES ======"
	find_features(three_utr_five_spectrum_data);
	score_features(three_utr_five_spectrum_data);

	#print "=== Three UTR, k-mer  length: 7 ===\n"
	#three_utr_seven_spectrum_data = get_spectrum_data(three_utr_mouse_ttp_file, 7, None, 495, 958, True)
	#print three_utr_seven_spectrum_data
	#r2 = run_svm(three_utr_seven_spectrum_data, 5, 1, 1)
	#print " ===== RESULTS ======"
	#print_results(r2, "three_utr_7_spectrum.csv");
	#print " ===== SCORING FEATURES ======"
	#find_features(three_utr_seven_spectrum_data);


def run_five_utr_spectrum_kernel():
	print "=== Five UTR, k-mer  length: 5 ===\n"
	## HuR binders - 3utr : pos: 273, neg: 286
	#five_utr_five_data = get_spectrum_data(five_utr_combined_file, 5, None, 273, 286, True)

	## HuR binders (random neg genes) - 3utr : pos: 273, neg: 256
	#five_utr_five_data = get_spectrum_data(five_utr_combined_file, 5, None, 273, 256, True)

	## HuR binders (random neg genes) - 5utr : pos: 273, neg: 2733
	five_utr_five_data = get_spectrum_data(five_utr_combined_file, 5, None, 273, 2733, True)
	print five_utr_five_data

	r1 = run_svm(five_utr_five_data, 5, 1, 1)
	print_results(r1, "five_utr_5_spectrum.csv");
	find_features(five_utr_five_data)

	print "=== Five UTR, k-mer  length: 7 ===\n"
	## HuR binders - 3utr : pos: 273, neg: 286
	#five_utr_seven_data = get_spectrum_data(five_utr_combined_file, 7, None, 273, 286, True)

	## HuR binders (random neg genes) - 3utr : pos: 273, neg: 256
	#five_utr_seven_data = get_spectrum_data(five_utr_combined_file, 7, None, 273, 256, True)

	## HuR binders (random neg genes) - 5utr : pos: 273, neg: 2733
	five_utr_seven_data = get_spectrum_data(five_utr_combined_file, 7, None, 273, 2733, True)
	print five_utr_seven_data

	r2 = run_svm(five_utr_seven_data, 5, 1, 1)
	print_results(r2, "five_utr_7_spectrum.csv");
	find_features(five_utr_seven_data)

def run_unspliced_spectrum_kernel():
	print "=== Unspliced , k-mer  length: 5 ===\n"
	unspliced_five_data = get_spectrum_data(unspliced_utr_combined_file, 5, None, 279, 288, True)
	print unspliced_five_data
	r1 = run_svm(unspliced_five_data, 5, 1, 1)
	print_results(r1, 'unspliced_5_spectrum.csv');
	find_features(unspliced_five_data)

	print "=== Unspliced, k-mer  length: 7 ===\n"
	unspliced_seven_data = get_spectrum_data(unspliced_utr_combined_file, 7, None, 279, 288, True)
	print unspliced_seven_data
	r2 = run_svm(unspliced_seven_data, 5, 1, 1)
	print_results(r2, 'unspliced_7_spectrum.csv');
	find_features(unspliced_seven_data)

def run_three_five_utr_5_spectrum_kernel():
	three_utr_data = get_spectrum_data(three_utr_combined_file, 5, None, 275, 286, True, '3utr');
	five_utr_data = get_spectrum_data(five_utr_combined_file, 5, None, 273, 286, True, '5utr');
	three_utr_data = three_utr_data.addFeatures(five_utr_data);
	r = run_svm(three_utr_data, 5, 1, 1)
	print_results(r, "three_five_utr_5_spectrum.csv");

def run_three_five_utr_7_spectrum_kernel():
	three_utr_data = get_spectrum_data(three_utr_combined_file, 7, None, 275, 286, True, '3utr');
	five_utr_data = get_spectrum_data(five_utr_combined_file, 7, None, 273, 286, True, '5utr');
	three_utr_data.addFeatures(five_utr_data);
	r = run_svm(three_utr_data, 5, 1, 1)
	print_results(r, "three_five_utr_7_spectrum.csv");

three_utr_homo_ttp_file = sys.argv[1];
k1 = int(sys.argv[2])
k2 = int(sys.argv[3])
posLen = int(sys.argv[4])
negLen = int(sys.argv[5])

if k2 == 0:
	k2 = 'None'
run_three_utr_spectrum_kernel(three_utr_homo_ttp_file, k1, k2, posLen, negLen)
#run_five_utr_spectrum_kernel()
#run_unspliced_spectrum_kernel()
#run_mixed_spectrum_kernel(5,10);

#run_three_five_utr_5_spectrum_kernel()
#run_three_five_utr_7_spectrum_kernel()

