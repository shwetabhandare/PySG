from generateYaml import *;
from SeqGenUtils import *;
from NoSignal import *
from Kmer import *
import sys;
import glob, os
import subprocess

directory = sys.argv[1]
#CreateConfFiles(directory, 'shuffle');
CreateConfFiles(directory, 'dirichlet');


for confFile in findFiles(directory, 'Kmer.yaml'):

	print "Generating sequences for : ", confFile;
	CreateNoSignalFastaFile(confFile);
	CreateFastaWithSignal(confFile)

os.chdir(directory)
for signalFile in glob.glob("Signal*.fa"):
	dremeDir = signalFile + "_DremeOut"
	noSignalFile = "No" + signalFile;
	print "Calling DREME for ", signalFile, ", ", noSignalFile;
	subprocess.call(["dreme", "-p", signalFile, "-n", noSignalFile, "-oc", dremeDir])

#	seqGen = SeqGen(conf);
#	filename = os.path.splitext(os.path.basename(confFile))[0]
#	posFastaFile = directory + "/" + filename + "_pos.fa";
#	seqGen.SetPosFileName(posFastaFile);
#
#	negFastaFile = directory + "/" + filename + "_neg.fa";
#
#	seqGen.SetNegFileName(negFastaFile)
#
#	seqGen.GenerateRandomSequences("negative", 0);
#	seqGen.GenerateRandomSequences("positive", 1);
#	#seqGen.embedMotifInSequence();

	# print "Positive Set: "
	# for seq in seqGen.GetPositiveSet():
	# 	print seq;

	# print "Negative Set: "
	# for seq in seqGen.GetNegativeSet():
	# 	print seq;

#	seqGen.writePositiveFile();
#	seqGen.writeNegativeFile();
