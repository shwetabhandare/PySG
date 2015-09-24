import sys
import os
import csv
import operator
import GenerateFlankingRegions
from GenerateFlankingRegions import *

#this script will write to two new .fa files, which need to be defined in the yml file.
#add output: TopKmerSeqFile.fa and NoTopKmerSeqFile.fa, for example


def SeparateSequences(SeqsWithKmers, SeqWithoutKmers, topKmersFile, seqFile):

	SeqWithTop10kmersFileWriter = open(SeqsWithKmers, 'w')
	SeqNoTop10kmersFileWriter = open(SeqWithoutKmers, 'w')
	TopKmersReader = open(topKmersFile, 'r')

	top10kmers = [] # will house the top ten lines from the
	headerList = [] # will house the headers containing top ten kmers

	for i in range(10):  #reads first 10 lines of topKmersfile
		line = TopKmersReader.next().strip() 
		line = line.split(",")  # since the file topKmersFile is a csv
		kmer = line[1]
		top10kmers.append(kmer)

	seqDict = fasta_read(seqFile) #create seq header keys with file from GenerateFlankingRegions.py


	for header, seq in seqDict.iteritems(): #for each header check:
		kmercount = 0
		faRows = ">" + header + "\n" + seq + "\n"    # to be written to the .fa file

		for kmer in top10kmers:             #each kmer
			kmerIndex = seq.find(kmer);
			if kmerIndex > 0:         #if kmer is in seq, then add header to list
				kmercount += 1;
			if kmercount > 0:  # i.e. top10 kmer is in sequence
				SeqWithTop10kmersFileWriter.write(faRows) #write to Top10Kmer File
			else:
				SeqNoTop10kmersFileWriter.write(faRows)

	SeqWithTop10kmersFileWriter.close()
	SeqNoTop10kmersFileWriter.close()
	TopKmersReader.close()



if __name__ == '__main__':
	import sys
	SeqsWithKmers = sys.argv[1]
	SeqWithoutKmers = sys.argv[2]
	topKmersFile = sys.argv[3]
	seqFile = sys.argv[4]
	SeparateSequences(SeqsWithKmers, SeqWithoutKmers, topKmersFile, seqFile);
