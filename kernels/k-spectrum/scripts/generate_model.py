from __future__ import division
import os
import numpy
import sys
import itertools
import csv
import re

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

from itertools import groupby
from collections import Counter
from PyML.utils import fasta
import operator

def longest_repitition(alist):
   try:
      return Counter(alist).most_common(1)[0][0]
   except IndexError:
      return None

def fasta_read(file_name) :
	"""read the sequence from a file in fasta format"""
	return [record.sequence for record in fasta.fasta_itr(file_name)]

def get_ATGC_counts(fasta_file):
   sequences = fasta_read(fasta_file);
   a_count = [];
   t_count = [];
   g_count = [];
   c_count = [];
   a_repeated_count = [];
   t_repeated_count = [];
   au_count = [];

   for s in sequences:
      seq_len = len(s);
      lower_s = s.lower();
      num_a = lower_s.count('a');
      num_t = lower_s.count('t');
      num_g = lower_s.count('g');
      num_c = lower_s.count('c');
      num_au_list = getAURich(lower_s);
      num_au = len(num_au_list)

      num_a_repeated_count, num_t_repeated_count  = getATRepeatedCount(s); 
      #print "ACOUNT: " + str(num_a_repeated_count) + ", TCOUNT: " + str(num_t_repeated_count);

      a_count.append(1.0 * (num_a/seq_len));
      t_count.append(1.0 * (num_t/seq_len));
      g_count.append(1.0 * (num_g/seq_len));
      c_count.append(1.0 * (num_c/seq_len));

      a_repeated_count.append(1.0 * (num_a_repeated_count/seq_len));
      t_repeated_count.append(1.0 * (num_t_repeated_count/seq_len));

      au_count.append(1.0  * (num_au/(seq_len - 1)));


   return a_count, t_count, g_count, c_count, a_repeated_count, t_repeated_count, au_count;

def get_repeated_AT_counts(l):
	sorted_list = sorted(l, key=operator.itemgetter(0,1))
	search = 'A';
	result = [element for element in sorted_list if element[0] == search]
	if len(result) == 0:
		ACount = 0;
	else:
		ACount = result[len(result) - 1][1];

	TCount = sorted_list[len(sorted_list) - 1][1];
	return ACount, TCount;

def getATRepeatedCount(seq):
	unsorted_list = [[k,len(list(g))] for k, g in groupby(seq)];
	l = sorted(unsorted_list, key = lambda x: int(x[1]));
	ACount, TCount = get_repeated_AT_counts(l); 
	return ACount, TCount;

# Look for count of AA, AT, TA, TT 
# Count number of times each is found, add them up.
# Follow the same normalization and pass the value to get_normalized_array.
# 
def getAURich(seq):
	p = re.compile( '(aa|at|ta|tt)')
	numAUElements = p.findall(seq);
	return numAUElements;

def get_normalized_array(a_repeated_count):
	mean_a_repeated = numpy.mean(a_repeated_count);

	var_a_repeated = numpy.var(a_repeated_count);

	sqrt_a_repeated = numpy.sqrt(var_a_repeated);

	print "Mean: " + str(mean_a_repeated)
	print "Varianece: " + str(var_a_repeated)
	print "SQRT: " + str(sqrt_a_repeated)

	print "COUNT: ", a_repeated_count
	normalized_a_repeated = (a_repeated_count - mean_a_repeated) / sqrt_a_repeated;
	print "Normalized: ", normalized_a_repeated


	return normalized_a_repeated;

# Create k-spectrum features.
def get_spectrum_data (fasta_file, k1, k2, pos_data_len, neg_data_len, normalize=False, prefix='None', repeatCount=False):
	# sequence data file, K1, K2 where k1 and k2 are k-mer lengths.
	if prefix == 'None':
		total_data = sequenceData.spectrum_data(fasta_file, k1, k2, normalize=normalize);
	else:
		total_data = sequenceData.spectrum_data(fasta_file, k1, k2, normalize=normalize, prefix=prefix);

	if repeatCount:
	   # add ATGC counts.
		a_count  = [];
		t_count  = [];
		g_count  = [];
		c_count  = [];
		au_count = [];

		a_count, t_count, c_count, g_count, a_repeated_count, t_repeated_count, au_count = get_ATGC_counts(fasta_file); 

		normalized_a_count = get_normalized_array(a_count);
		normalized_t_count = get_normalized_array(t_count);
		normalized_c_count = get_normalized_array(c_count);
		normalized_g_count = get_normalized_array(g_count);
		normalized_au_count = get_normalized_array(au_count);

		total_data.addFeature('aCount', normalized_a_count);
		total_data.addFeature('tCount', normalized_t_count); 
		total_data.addFeature('gCount', normalized_g_count); 
		total_data.addFeature('cCount', normalized_c_count); 
		total_data.addFeature('auCount', normalized_au_count); 


		normalized_a_repeated = get_normalized_array(a_repeated_count);
		normalized_t_repeated = get_normalized_array(t_repeated_count);

		total_data.addFeature('aRepeated', normalized_a_repeated); 
		total_data.addFeature('tRepeated', normalized_t_repeated); 


	#Generate labels
	L = ["+1"] * (pos_data_len);
	L[pos_data_len:neg_data_len] = ["-1"] *(neg_data_len);
	total_data.attachLabels(Labels(L));

	if normalize:
		total_data.normalize();

	print total_data
	return total_data;

