#!/usr/bin/env python

# A simple Python n-gram calculator. 
#
# Given an arbitrary string, and the value of n
# as the size of the n-gram (int), this code 
# snip will show you the results, sorted from
# most to least frequently occuring n-gram.
#
# The 'sort by value' operation for the dict 
# follows the PEP 265 recommendation.
#
# Quick start:
#
# from pyngram import calc_ngram
# calc_ngram('bubble bobble, bubble bobble, bubble bobble', 3)
#
# Enjoy!
#
from __future__ import division

__version__ = '1.0'
__author__ = 'Jay Liew' # @jaysern from @websenselabs
__license__ = 'MIT'

from operator import itemgetter
import SeqGenUtils;
import sys

def getNgramListForSeq(inputstring, nlen):
    empty_list = [];
    if nlen < 1:
        print "Uh, n-grams have to be of size 1 or greater. Makes no sense to have a 0 length n-gram."
        return empty_list;

    if len(inputstring) < 1:
        print "umm yeah, ... the inputstring has to be longer than 1 char ;)"
        return empty_list;

    # now, fish out the n-grams from the input string
    seq_len = len(inputstring)
    ngram_list = [inputstring[x:x+nlen] for x in xrange(seq_len)]

    return ngram_list;

def getNGramFreqForCombinedList(nGramCombinedList, nlen):
    ngram_freq = {}

    for nGramList in nGramCombinedList:
        for n in nGramList: # collect the distinct n-grams and count
            if len(n) == nlen:
                if n in ngram_freq:
                    ngram_freq[n] += 1 
                else:
                    ngram_freq[n] = 1 # human counting numbers start at 1
            else:
                # these are the remainder strings not long enough to fit n
                pass
    return sorted(ngram_freq.iteritems(), key=itemgetter(1), reverse=True)


def ComputeNgramFrequencyAndProbability(seqFile, nLen):

    seqDict = SeqGenUtils.fasta_read(seqFile);
    nGramCombinedList = list();
    total_seq_length = 0;

    for header, sequence in seqDict.iteritems():
        ngram_list = getNgramListForSeq(sequence, nLen)
        if len(ngram_list) > 0:
            nGramCombinedList.append(ngram_list)
            total_seq_length = total_seq_length + len(sequence)
        else: 
            print "Found empty sequence for ", header;

    ngram_freq = getNGramFreqForCombinedList(nGramCombinedList, nLen);
    ngram_freq = dict(ngram_freq)
    ngram_prob = {}

    for ngram, frequency in ngram_freq.iteritems():
        ngram_prob[ngram] = round(frequency/total_seq_length, 4);

    # print str(nLen) + "-" + "gram frequences: \n", ngram_freq;
    #print str(nLen) + "-" + "gram probabilities: \n", ngram_prob;
    #print "Total Nucleotides: ", str(total_seq_length) 

    return ngram_freq, ngram_prob;

if __name__ == '__main__':
    import sys

    seqFile = sys.argv[1]
    nLen = int(sys.argv[2])
    ComputeNgramFrequencyAndProbability(seqFile, nLen)


