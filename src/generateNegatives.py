import Conf;
import SeqGen;
from Conf import *;
from SeqGen import *;
import generateYaml;
from generateYaml import *;
from SeqGenUtils import *;

import sys;

directory = sys.argv[1]
CreateConfFiles(directory);
for confFile in findFiles(directory, '*.yml'):

	conf = Conf(confFile);

	seqGen = SeqGen(conf);
	filename = os.path.splitext(os.path.basename(confFile))[0]

	negFastaFile = directory + "/" + filename + "_neg.fa";

	seqGen.SetNegFileName(negFastaFile)

	seqGen.GenerateRandomSequences("negative", 0);

	#seqGen.embedMotifInSequence();

	# print "Positive Set: "
	# for seq in seqGen.GetPositiveSet():
	# 	print seq;

	# print "Negative Set: "
	# for seq in seqGen.GetNegativeSet():
	# 	print seq;

	seqGen.writeNegativeFile();

'''
LenMatchNegatives.py

Requires Python 2.7 and fasta.py in the same directory

Step 1: Create dictionary of negative set of sequences.

Step 2: Create a set of sequences for positive set.

Creates length-matched negative sequences for each sequence in a list of positive transcripts

Takes two command line arguments:
Positive File Name (to match to, in fasta format)     Negative File Name (to draw sequences from, in fasta format)

Outputs file called NON + Positive File Name
'''
