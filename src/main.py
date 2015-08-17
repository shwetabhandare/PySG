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
	posFastaFile = directory + "/" + filename + "_pos.fa";
	seqGen.SetPosFileName(posFastaFile);

	negFastaFile = directory + "/" + filename + "_neg.fa";

	seqGen.SetNegFileName(negFastaFile)

	seqGen.GenerateRandomSequences("negative", 0);
	seqGen.GenerateRandomSequences("positive", 1);
	#seqGen.embedMotifInSequence();

	# print "Positive Set: "
	# for seq in seqGen.GetPositiveSet():
	# 	print seq;

	# print "Negative Set: "
	# for seq in seqGen.GetNegativeSet():
	# 	print seq;

	seqGen.writePositiveFile();
	seqGen.writeNegativeFile();
