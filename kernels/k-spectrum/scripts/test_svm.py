import os
import numpy
import sys

from PyML.classifiers.svm import SVM
from PyML.containers import ker,labels
from PyML.containers.kernelData import KernelData
from PyML.containers.sequenceData import SequenceData
from PyML.containers import vectorDatasets
from PyML import *

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


## Main

seqFile = sys.argv[1];
posLen = int(sys.argv[2])
negLen = int(sys.argv[3])

#Iterate through the kmers.
seq_data = get_spectrum_data(seqFile, 5, 7, posLen, negLen, True);
#seq_data.save("seq.txt");

Cs = [ 10**x for x in xrange( -5, 5 ) ]
for C in Cs:
   print C;

for C in Cs:
   print C;
   s = svm.SVM(optimizer = 'liblinear', C=C);
   s.train(seq_data);

