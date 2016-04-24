
from generateYaml import *;
from SeqGenUtils import *;
from NoSignal import *
from Kmer import *


import compareKmers;

import RunComputationalMethods;

import sys;
import glob;




def GenerateFastaFiles(directory):
	for confFile in findFiles(directory, '*.yml'):
		print "Generating sequences for : ", confFile;
		CreateNoSignalFastaFile(confFile);
		CreateFastaWithSignal(confFile)

def RunComputationalTools(directory):
	os.chdir(directory)
	for signalFile in glob.glob("Signal*.fa"):
		noSignalFile = "No" + signalFile;

		dremeResultDir, realKmersCsvFile = RunComputationalMethods.RunDremeAndGetResults(signalFile, noSignalFile);
		kspectrumResultDir, realKmersCsvFile = RunComputationalMethods.RunKspectrumAndGetResults(signalFile, noSignalFile);
		RunComputationalMethods.CopyResults(signalFile, noSignalFile, realKmersCsvFile, dremeResultDir, kspectrumResultDir)

if __name__ == "__main__":
	import sys

	directory = sys.argv[1]
	randomType = sys.argv[2];
	#CreateConfFiles(directory, randomType);
	#GenerateFastaFiles(directory);
	RunComputationalTools(directory);





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
