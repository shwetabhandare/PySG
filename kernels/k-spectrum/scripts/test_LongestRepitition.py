import os
import sys
from itertools import groupby
from collections import Counter
from PyML.utils import fasta
import operator

def fasta_read(file_name) :
	"""read the sequence from a file in fasta format"""
	return [record.sequence for record in fasta.fasta_itr(file_name)]

def fasta_read_header(file_name) :
	"""read the sequence from a file in fasta format"""
	return [record.header for record in fasta.fasta_itr(file_name)]

def getATRepeatedCounts(l):
	sorted_list = sorted(l, key=operator.itemgetter(0,1))
	search = 'A';
	result = [element for element in sorted_list if element[0] == search]
	ACount = result[len(result) - 1];
	TCount = sorted_list[len(sorted_list) - 1];
	return ACount, TCount;

def longest_repitition(seq):
	unsorted_list = [[k,len(list(g))] for k, g in groupby(seq)];
	l = sorted(unsorted_list, key = lambda x: int(x[1]));
	ACount, TCount = getATRepeatedCounts(l); 
	return ACount, TCount;

def addLongestRepetition(fasta_file):
	longestCounts = [];
	sequences = fasta_read(fasta_file);
	for s in sequences: 
		ACount, TCount = longest_repitition(s); 
		l = [ACount, TCount]; 
		longestCounts.append(l); 
	return longestCounts;

fasta_file = sys.argv[1];
longestCounts = addLongestRepetition(fasta_file);
for l in longestCounts:
   print l;

