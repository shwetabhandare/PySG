import os
import sys
import itertools
import csv
import string
from PyML.utils import fasta

def fasta_read(file_name) :
	"""read the sequence from a file in fasta format"""
	return [record.sequence for record in fasta.fasta_itr(file_name)]

def fasta_read_header(file_name) :
	"""read the sequence from a file in fasta format"""
	return [record.header for record in fasta.fasta_itr(file_name)]

def get_seq_length_feature(fasta_file):
	seqLengths = [];
	sequences = fasta_read(fasta_file);
	for s in sequences:
		seqLengths.append(len(s));
	return seqLengths;


def get_motif_positions(haystack, needle):
	i = 0
	l = list()
	try:
		while True:
			i = haystack.index(needle, i)
			l.append(i)
			#print i, haystack[i:i+len(needle)];
			i += len(needle)
	except ValueError:
		pass

	return l;

def getSequenceFeatures(fasta_file):
	seqLengths = [];
	motifPos = [];
	motifDist = [];
	nonamer = "TTATTTATT";
	sequences = fasta_read(fasta_file);
	for s in sequences:
		seqLengths.append(len(s));
		l = get_motif_positions(s, nonamer);
		motifPos.append(l);
		if len(l) > 1:
			d = [];
			for x,y in zip(l,l[1:]):
				dist = y - x;
				d.append(dist);
			motifDist.append(d);
		else:
			d = [];
			d.append(-99);
			motifDist.append(d);

	return seqLengths, motifPos, motifDist;


fasta_file = sys.argv[1];
seqLengths, motifPos, motifDist = getSequenceFeatures(fasta_file);
headers = fasta_read_header(fasta_file);
#for (l,d,h,len) in zip(motifPos, motifDist, headers, seqLengths):
for (h, l, mp, md) in zip(headers, seqLengths, motifPos, motifDist):
	fields = h.rsplit("|");
	fields[-1] = fields[-1].strip();
	#print fields[2], l, len(mp), md
	if len(mp) > 0:
		#print fields;
		#print fields[2], l, d
		#print h, l, d
		if len(fields) == 3:
			print "%s, %s, %d, %s"%(fields[2], l, len(mp), md)
		else:
			print "%s, %s, %d, %s"%(fields[0], l, len(mp), md)
		#print l, d
