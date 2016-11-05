
from generateYaml import *;
from SeqGenUtils import *;
from NoSignal import *
from Kmer import *
import datetime;
import uuid



import compareKmers;

import RunComputationalMethods;

import sys;
import glob;


def GetCurrentDateTimeStr():
	import time
	dateTime= datetime.datetime.now();
	timestr = dateTime.strftime("%Y%m%d-%H%M%S-%f")
	timestr = timestr + "-" + str(dateTime.microsecond)
	return timestr;


def GenerateFastaFiles(directory):
	for confFile in findFiles(directory, '*.yml'):
		print "Generating sequences for : ", confFile;
		# Create positive set of sequences with signal embedded. 
		PosSeqDict = CreateFastaWithSignal(confFile)
		# Shuffle the sequences to generate negative set.
		NegSeqDict = ShuffleToCreateNoSignalSequences(PosSeqDict, confFile);

def RunComputationalTools(directory, ntBased=True):
	os.chdir(directory)
	for signalFile in findFiles(directory, "Signal*.fa"):
		dirName = os.path.dirname(signalFile);
		signalFileName = os.path.basename(signalFile);
		noSignalFile = dirName + "/No" + signalFileName;
		
		print "Signal File: ", signalFileName
		dremeResultDir, realKmersCsvFile = RunComputationalMethods.RunDremeAndGetResults(signalFile, noSignalFile, ntBased);
		kspectrumResultDir, realKmersCsvFile = RunComputationalMethods.RunKspectrumAndGetResults(signalFile, noSignalFile, ntBased);
		#RunComputationalMethods.CopyResults(signalFile, noSignalFile, realKmersCsvFile, dremeResultDir, kspectrumResultDir)


if __name__ == "__main__":
	import sys

	confFile = sys.argv[1];
	if len(sys.argv) == 1
		ntBased = True
	else:
		ntBased = False
	uuid_to_append = str(uuid.uuid4())
	generator = YamlFastaGenerator(confFile, uuid_to_append);
	targetDir = generator.GetTargetDir();
	generator.CreateConfFiles();
	GenerateFastaFiles(targetDir);
	RunComputationalTools(targetDir, ntBased);
	#ParseResultsAndGenerateGraphs(targetDir)

