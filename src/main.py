
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
		CreateNoSignalFastaFile(confFile);
		CreateFastaWithSignal(confFile)

def RunComputationalTools(directory):
	os.chdir(directory)
	for signalFile in findFiles(directory, "Signal*.fa"):
		dirName = os.path.dirname(signalFile);
		signalFileName = os.path.basename(signalFile);
		noSignalFile = dirName + "/No" + signalFileName;
		
		print "Signal File: ", signalFileName
		dremeResultDir, realKmersCsvFile = RunComputationalMethods.RunDremeAndGetResults(signalFile, noSignalFile);
		kspectrumResultDir, realKmersCsvFile = RunComputationalMethods.RunKspectrumAndGetResults(signalFile, noSignalFile);
		#RunComputationalMethods.CopyResults(signalFile, noSignalFile, realKmersCsvFile, dremeResultDir, kspectrumResultDir)

if __name__ == "__main__":
	import sys

	confFile = sys.argv[1];
	uuid_to_append = str(uuid.uuid4())
	generator = YamlFastaGenerator(confFile, uuid_to_append);
	targetDir = generator.GetTargetDir();
	generator.CreateConfFiles();
	GenerateFastaFiles(targetDir);
	RunComputationalTools(targetDir);

