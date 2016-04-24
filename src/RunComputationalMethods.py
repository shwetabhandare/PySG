from __future__ import division

import glob, os
import subprocess
import shutil
import compareKmers

def RunKspectrumKernel(signalFile, noSignalFile):
	baseNameOfSignalFile = os.path.splitext(os.path.basename(signalFile))[0]
	outDir = baseNameOfSignalFile + "_kspectrum"
	makeFile = "/projects/bhandare/workspace/PySG/scripts/k-spectrum-pipeline/Makefile.BasicKspectrum"
	posFileArg = "pos_file=" + signalFile;
	negFileArg = "neg_file=" + noSignalFile;
	prefixArg = "conf_prefix=" + baseNameOfSignalFile 
	resultDirArg = "result_dir=" + outDir
	subprocess.call(["make", "-f", makeFile, posFileArg, negFileArg, prefixArg, resultDirArg])

	predictedKspectrumFile = outDir + "/" + baseNameOfSignalFile + "_Features.dat"
	realKmersCsvFile = os.path.splitext(signalFile)[0] + ".kmers"

	return predictedKspectrumFile, realKmersCsvFile, outDir;


def RunDreme(signalFile, noSignalFile):
	dremeDir = signalFile + "_DremeOut"
	print "Calling DREME for ", signalFile, ", ", noSignalFile;
	subprocess.call(["dreme", "-p", signalFile, "-n", noSignalFile, "-oc", dremeDir])

	predictedDremeFile = dremeDir + "/" + "dreme.txt";
	realKmersCsvFile = os.path.splitext(signalFile)[0] + ".kmers"

	return predictedDremeFile, realKmersCsvFile, dremeDir;

def GetSensitivityAndPPV(numTP, numFP, numFN):
	if ((numTP + numFN) == 0):
		sensitivity = 0;
	else:
		sensitivity = numTP / (numTP + numFN)
	
	if ((numTP + numFP) == 0):
		ppv = 0;
	else:
		ppv =  numTP / (numTP + numFP)

	return sensitivity, ppv;

def ComputeDremeResults(predictedDremeFile, realKmersCsvFile, signalFile, noSignalFile):

	numTP, numFP, numFN = compareKmers.CompareDremeKmers(realKmersCsvFile, predictedDremeFile,
		signalFile, noSignalFile, None)

	sensitivity, ppv = GetSensitivityAndPPV(numTP, numFP, numFN);

	return sensitivity, ppv;

def ComputeKspectrumResults(predictedKspectrumFile, realKmersCsvFile, signalFile, noSignalFile):
	#numTP, numFP, numFN = compareKmers.CompareKspectrumKmers(realKmersCsvFile, predictedKspectrumFile,
		#signalFile, noSignalFile)
	numTP, numFP, numFN = compareKmers.CompareKspectrumPredictedKmers(realKmersCsvFile, predictedKspectrumFile, 
	signalFile, noSignalFile)

	sensitivity, ppv = GetSensitivityAndPPV(numTP, numFP, numFN);

	return sensitivity, ppv;

def WriteResults(sensitivity, ppv, signalFile, resultFile):
	target = open(resultFile, 'w')
	resultStr = "Sensitivity:" + str(sensitivity) + ",PPV:" + str(ppv);
	target.write(resultStr)
	target.close();


def RunDremeAndGetResults(signalFile, noSignalFile):

	predictedDremeFile, realKmersCsvFile, dremeResultDir = RunDreme(signalFile, noSignalFile);
	sensitivity, ppv = ComputeDremeResults(predictedDremeFile, realKmersCsvFile, signalFile, noSignalFile);
	print "DREME: ", str(sensitivity), str(ppv);
	dremeResultFile = dremeResultDir + "/" + signalFile + "_dreme.results"
	WriteResults(sensitivity, ppv, signalFile, dremeResultFile)

	return dremeResultDir, realKmersCsvFile;

def RunKspectrumAndGetResults(signalFile, noSignalFile):

	predictedKspectrumFile, realKmersCsvFile, kspectrumResultDir = RunKspectrumKernel(signalFile, noSignalFile)
	sensitivity, ppv = ComputeKspectrumResults(predictedKspectrumFile, realKmersCsvFile, signalFile, noSignalFile);
	kSpectrumResultFile = kspectrumResultDir + "/" + signalFile + "_kspectrum.results"
	WriteResults(sensitivity, ppv, signalFile, kSpectrumResultFile)
	print "K_SPECTRUM: ", str(sensitivity), str(ppv);

	return kspectrumResultDir, realKmersCsvFile;

def CopyResults(signalFile, noSignalFile, realKmersCsvFile, dremeResultDir, kspectrumResultDir):
	# Create an enclosing directory for the signal/nosignal file
	baseNameOfSignalFile = os.path.splitext(os.path.basename(signalFile))[0]
	baseNameOfNoSignalFile = os.path.splitext(os.path.basename(noSignalFile))[0]

	enclosingDir = baseNameOfSignalFile + "_" + baseNameOfNoSignalFile;
	if not os.path.exists(enclosingDir):
		os.makedirs(enclosingDir)

		shutil.move(signalFile, enclosingDir)
		shutil.move(noSignalFile, enclosingDir)
		shutil.move(kspectrumResultDir, enclosingDir)
		shutil.move(dremeResultDir, enclosingDir)
		shutil.move(realKmersCsvFile, enclosingDir)
