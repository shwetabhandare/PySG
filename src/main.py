from __future__ import division

from generateYaml import *;
from SeqGenUtils import *;
from NoSignal import *
from Kmer import *


import compareKmers;


import sys;
import glob, os
import subprocess

def RunKspectrumKernel(signalFile, noSignalFile):
	baseNameOfSignalFile = os.path.splitext(os.path.basename(signalFile))[0]
	outDir = baseNameOfSignalFile + "_kspectrum"
	makeFile = "/projects/bhandare/workspace/PySG/scripts/k-spectrum-pipeline/Makefile.BasicKspectrum"
	posFileArg = "pos_file=" + signalFile;
	negFileArg = "neg_file=" + noSignalFile;
	prefixArg = "conf_prefix=" + baseNameOfSignalFile 
	resultDirArg = "result_dir=" + outDir
	subprocess.call(["make", "-f", makeFile, posFileArg, negFileArg, prefixArg, resultDirArg])



def RunDreme(signalFile, noSignalFile):
	dremeDir = signalFile + "_DremeOut"
	print "Calling DREME for ", signalFile, ", ", noSignalFile;
	subprocess.call(["dreme", "-p", signalFile, "-n", noSignalFile, "-oc", dremeDir])

	predictedDremeFile = dremeDir + "/" + "dreme.txt";
	realKmersCsvFile = os.path.splitext(signalFile)[0] + ".kmers"

	return predictedDremeFile, realKmersCsvFile;

def ComputeDremeResults(predictedDremeFile, realKmersCsvFile, signalFile, noSignalFile):

	numTP, numFP, numFN = compareKmers.CompareDremeKmers(realKmersCsvFile, predictedDremeFile,
		signalFile, noSignalFile, None)

	sensitivity = numTP / (numTP + numFN)
	ppv = numTP / (numTP + numFP)

	return sensitivity, ppv;

def WriteResults(sensitivity, ppv, resultFile):
	target = open(resultFile, 'w')
	resultStr = "Sensitivity:" + str(sensitivity) + ",PPV:" + str(ppv);
	target.write(resultStr)
	target.close();


def GenerateFastaFiles(directory):
	for confFile in findFiles(directory, 'Kmer90.yaml'):
		print "Generating sequences for : ", confFile;
		CreateNoSignalFastaFile(confFile);
		CreateFastaWithSignal(confFile)

def RunComputationalTools(directory):
	os.chdir(directory)
	for signalFile in glob.glob("Signal*90.fa"):
		noSignalFile = "No" + signalFile;
		resultFile = signalFile + ".results"

		#predictedDremeFile, realKmersCsvFile = RunDreme(signalFile, noSignalFile);
		#sensitivity, ppv = ComputeDremeResults(predictedDremeFile, realKmersCsvFile, signalFile, noSignalFile);
		#WriteResults(sensitivity, ppv, resultFile)

		RunKspectrumKernel(signalFile, noSignalFile)


if __name__ == "__main__":
	import sys
	directory = sys.argv[1]
	CreateConfFiles(directory, 'dirichlet');
	GenerateFastaFiles(directory);
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
