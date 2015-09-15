import sys
import os
import csv
import operator
import GenerateFlankingRegions
from GenerateFlankingRegions import *

#this script will write to two new .fa files, which need to be defined in the yml file.
#add output: TopKmerSeqFile.fa and NoTopKmerSeqFile.fa, for example

CreateFlankingRegions() # generate the top kmers file and the kmerflanking region file


confMap = ReadConfigFile();
seqFile = confMap['input']['dataFile']
kmerFile = confMap['input']['featureKmers']
topKmersFile = confMap['output']['topKmersFile']

seqDict = fasta_read(seqFile) #create seq header keys with file from GenerateFlankingRegions.py
TopKmerSeqFile = confMap['output']['TopKmerSeqFile']  # create an output filename in yml to write the sequences with top 10 kmers. will be a .fa file
NoTopKmerSeqFile = confMap['output']['FilteredSeqFile'] # create output filename in yml where the sequences that do not contain top 10 kmers will be written.  will be .fa file


SeqWithTop10kmersFileWriter = open(TopKmerSeqFile, 'w')
SeqNoTop10kmersFileWriter = open(NoTopKmerSeqFile, 'w')
TopKmersReader = open(topKmersFile, 'r')



top10kmers = [] # will house the top ten lines from the
headerList = [] # will house the headers containing top ten kmers

for i in range(10):  #reads first 10 lines of topKmersfile
    line = TopKmerReaders.next().strip()
    line = line.split(",")  # since the file topKmersFile is a csv
    kmer = line[1]
    top10kmers.append(kmer)



for header, seq in seqDict.iteritems(): #for each header check:
    kmercount = 0
    faRows = ">" + header + "\n" + seq + "\n"    # to be written to the .fa file
        
    for kmer in top10kmers:             #each kmer
        if seq.count(kmer) > 0:         #if kmer is in seq, then add header to list
            kmercount += kmercount
    if kmercount > 0:  # i.e. top10 kmer is in sequence
        SeqWithTop10kmersFileWriter.write(faRows) #write to Top10Kmer File
    
    else:
        SeqNoTop10kmersFileWriter.write(faRows)

SeqWithTop10kmersFileWriter.close()
SeqNoTop10kmersFileWriter.close()
TopKmersReader.close()

