from __future__ import division

import glob, os
import subprocess
import shutil
import compareKmers
import dremeSequenceBasedComparison

def RunKspectrumKernel(signalFile, noSignalFile):

	signalDirName, signalFileName, outDir = GetResultDirName(signalFile, "_kspectrum");
	noSignalDirName, noSignalFileName, noSignalOut = GetResultDirName(noSignalFile, "_kspectrum");

	print "PRED: ", signalDirName, ", Signal: ", signalFileName, ", OUT: ", outDir;

	baseNameOfSignalFile = os.path.splitext(os.path.basename(signalFileName))[0]

	makeFile = "/projects/bhandare/workspace/PySG/scripts/k-spectrum-pipeline/Makefile.BasicKspectrum"
	posFileArg = "pos_file=" +  signalFileName;
	negFileArg = "neg_file="  + noSignalFileName;
	prefixArg = "conf_prefix=" + baseNameOfSignalFile 
	resultDirArg = "result_dir=" + outDir
	if not os.path.exists(outDir):
		os.makedirs(outDir);
	subprocess.call(["make", "-f", makeFile, posFileArg, negFileArg, prefixArg, resultDirArg])


	predictedKspectrumFile = outDir + "/" + baseNameOfSignalFile + "_Features.dat"
	realKmersCsvFile = signalDirName + "/" + os.path.splitext(signalFileName)[0] + ".kmers"

	return predictedKspectrumFile, realKmersCsvFile, outDir;

def GetResultDirName(signalFile, suffixStr):
	signalDirName = os.path.dirname(signalFile);
	signalFileName =  os.path.basename(signalFile);

	resultDir =  signalDirName + "/" + signalFileName + suffixStr;

	return signalDirName, signalFileName, resultDir;

def RunDreme(signalFile, noSignalFile):
	signalDirName, signalFileName, dremeDir =  GetResultDirName(signalFile, "_DremeOut");
	print "Calling DREME for ", signalFile, ", ", noSignalFile, ", ", dremeDir;
	
	subprocess.call(["dreme", "-p", signalFile, "-n", noSignalFile, "-oc", dremeDir])

	predictedDremeFile = dremeDir + "/" + "dreme.txt";
	realKmersCsvFile = signalDirName + "/" + os.path.splitext(signalFileName)[0] + ".kmers"

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
		signalFile, noSignalFile, 0)

	sensitivity, ppv = GetSensitivityAndPPV(numTP, numFP, numFN);

	return sensitivity, ppv;

def ComputeKspectrumResults(predictedKspectrumFile, realKmersCsvFile, signalFile, noSignalFile, numKmers):
	#numTP, numFP, numFN = compareKmers.CompareKspectrumKmers(realKmersCsvFile, predictedKspectrumFile,
		#signalFile, noSignalFile)
	numTP, numFP, numFN = compareKmers.CompareKspectrumPredictedKmers(realKmersCsvFile, predictedKspectrumFile, 
	signalFile, noSignalFile, numKmers)

	sensitivity, ppv = GetSensitivityAndPPV(numTP, numFP, numFN);

	return sensitivity, ppv;

def WriteResults(sensitivity, ppv, signalFile, resultFile):
	target = open(resultFile, 'w')
	resultStr = "Sensitivity:" + str(sensitivity) + ",PPV:" + str(ppv);
	target.write(resultStr)
	target.close();


def RunDremeAndGetResults(signalFile, noSignalFile, nucleotideBased=True):

	predictedDremeFile, realKmersCsvFile, dremeResultDir = RunDreme(signalFile, noSignalFile);
	if nucleotideBased:
		sensitivity, ppv = ComputeDremeResults(predictedDremeFile, realKmersCsvFile, signalFile, noSignalFile);
		print "NT BASED: DREME: ", str(sensitivity), str(ppv);
	else:
		sensitivity, ppv = computeSequenceBasedDREMEResults(predictedDremeFile, realKmersCsvFile, signalFile, noSignalFile);
		print "SEQ BASED: DREME: ", str(sensitivity), str(ppv);

	dremeResultFile = dremeResultDir + "/" + os.path.basename(signalFile) + "_dreme.results"
	WriteResults(sensitivity, ppv, signalFile, dremeResultFile)

	return dremeResultDir, realKmersCsvFile;

def ComputeKspectrumResultForNumKmers(predictedKspectrumFile, realKmersCsvFile, kspectrumResultDir, signalFile, noSignalFile, numKmers, ntBased):
	if ntBased:
		sensitivity, ppv = ComputeKspectrumResults(predictedKspectrumFile, realKmersCsvFile, signalFile, noSignalFile, numKmers);
	else:
		sensitivity, ppv = kspectrumSequenceBasedComparison.computeSequenceBasedKspectrumResults(predictedKspectrumFile, realKmersCsvFile, 
		signalFile, noSignalFile, numKmers):

	kSpectrumResultFile = kspectrumResultDir + "/" + os.path.basename(signalFile) + "_kspectrum_" + str(numKmers) + ".results"
	WriteResults(sensitivity, ppv, signalFile, kSpectrumResultFile)
	print "K_SPECTRUM: ", str(sensitivity), str(ppv);
	return kspectrumResultDir, realKmersCsvFile;

def RunKspectrumAndGetResults(signalFile, noSignalFile, ntBased=True):

	predictedKspectrumFile, realKmersCsvFile, kspectrumResultDir = RunKspectrumKernel(signalFile, noSignalFile)

	numKmers = 25;
	kspectrumResultDir, realKmersCsvFile = ComputeKspectrumResultForNumKmers(predictedKspectrumFile, realKmersCsvFile, kspectrumResultDir, signalFile, noSignalFile, numKmers, ntBased);

	numKmers = 50;
	kspectrumResultDir, realKmersCsvFile = ComputeKspectrumResultForNumKmers(predictedKspectrumFile, realKmersCsvFile, kspectrumResultDir, signalFile, noSignalFile, numKmers, ntBased);

	numKmers = 100;
	kspectrumResultDir, realKmersCsvFile = ComputeKspectrumResultForNumKmers(predictedKspectrumFile, realKmersCsvFile, kspectrumResultDir, signalFile, noSignalFile, numKmers, ntBased);


	return kspectrumResultDir, realKmersCsvFile;

def CopyResults(signalFile, noSignalFile, realKmersCsvFile, dremeResultDir, kspectrumResultDir):
	# Create an enclosing directory for the signal/nosignal file
	signalDirName, signalFileName, outDir = GetResultDirName(signalFile, "");
	noSignalDirName, noSignalFileName, outDir = GetResultDirName(noSignalFile, "");

	enclosingDir = signalFileName + "_" + signalFileName;
	if not os.path.exists(signalDirName + "/" + enclosingDir):
		os.makedirs(enclosingDir)
		yamlFile = signalFileName[7:-3] + ".yml"
		shutil.move(signalDirName + "/" + yamlFile, enclosingDir)
		shutil.move(signalDirName + "/" + signalFileName, enclosingDir)
		shutil.move(signalDirName + "/" + noSignalFileName, enclosingDir)
		shutil.move(kspectrumResultDir, enclosingDir)
		shutil.move(dremeResultDir, enclosingDir)
		shutil.move(realKmersCsvFile, enclosingDir)
