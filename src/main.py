import Conf;
import SeqGen;
from Conf import *;
from SeqGen import *;
import generateYaml;
from generateYaml import *;
import sys;

directory = sys.argv[1]
CreateConfFiles(directory);

import os, fnmatch

def findFiles (path, filter):
   for root, dirs, files in os.walk(path):
      for file in fnmatch.filter(files, filter):
         yield os.path.join(root, file)

for confFile in findFiles(directory, '*.yml'):

	conf = Conf(confFile);

	seqGen = SeqGen(conf);
	filename = os.path.splitext(os.path.basename(confFile))[0]
	posFastaFile = directory + "/" + filename + "_pos.fa";
	seqGen.SetPosFileName(posFastaFile);

	negFastaFile = directory + "/" + filename + "_neg.fa";

	seqGen.SetNegFileName(negFastaFile)

	seqGen.GenerateRandomSequences("negative");
	seqGen.GenerateRandomSequences("positive");
	seqGen.embedMotifInSequence();

	# print "Positive Set: "
	# for seq in seqGen.GetPositiveSet():
	# 	print seq;

	# print "Negative Set: "
	# for seq in seqGen.GetNegativeSet():
	# 	print seq;

	seqGen.writePositiveFile();
	seqGen.writeNegativeFile();
