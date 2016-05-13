
from generateYaml import *;
from SeqGenUtils import *;
from NoSignal import *
from Kmer import *


import compareKmers;

import RunComputationalMethods;

import sys;
import glob;


def GetCurrentDateTimeStr():
	import time
	timestr = time.strftime("%Y%m%d-%H%M%S")
	return timestr;


def GenerateFastaFiles(directory):
	for confFile in findFiles(directory, '*.yml'):
		print "Generating sequences for : ", confFile;
		CreateNoSignalFastaFile(confFile);
		CreateFastaWithSignal(confFile)

def RunComputationalTools(directory):
	os.chdir(directory)
	for signalFile in findFiles(directory, "Signal*.fa"):
		dirName = os.path.dirname(signalFile);
		signalFileName = os.path.basename(signalFile);
		noSignalFile = dirName + "/No" + signalFileName;
		dremeResultDir, realKmersCsvFile = RunComputationalMethods.RunDremeAndGetResults(signalFile, noSignalFile);
		kspectrumResultDir, realKmersCsvFile = RunComputationalMethods.RunKspectrumAndGetResults(signalFile, noSignalFile);
		#RunComputationalMethods.CopyResults(signalFile, noSignalFile, realKmersCsvFile, dremeResultDir, kspectrumResultDir)

if __name__ == "__main__":
	import sys

	confFile = sys.argv[1];
	generator = YamlFastaGenerator(confFile);
	targetDir = generator.GetTargetDir();
	timestr = GetCurrentDateTimeStr();

	newTargetDir = targetDir + "/" + timestr;
	generator.CreateConfFiles(timestr);
	GenerateFastaFiles(newTargetDir);
	RunComputationalTools(newTargetDir);

