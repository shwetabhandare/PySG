import sys
import os
import csv
import glob

import numpy as np;
import math;

writeTitleDreme = False;
writeTitleKspectrum = False;

def getColumn(filename, column):
	results = csv.reader(open(filename), delimiter=",")
	return [result[column] for result in results]

def appendTitle(fileName):
	with open(fileName, 'wb') as csvfile:	
		spamwriter = csv.writer(csvfile, delimiter=',');
		spamwriter.writerow(["label", "value"]);

def appendResult(fileName, label, value):
	with open(fileName, 'a') as csvfile:	
		spamwriter = csv.writer(csvfile, delimiter=',');
		spamwriter.writerow([label, value])

def writeTitle(filepath, label, sensitivityResultFile, ppvResultFile):
	global writeTitleDreme;
	global writeTitleKspectrum;


	if writeTitleDreme == False:
		#print "Appending title sensitivity: ", sensitivityResultFile;
		appendTitle(sensitivityResultFile);
		writeTitleDreme = True;

	if writeTitleKspectrum == False:
		#print "Appending title ppv: ", ppvResultFile;
		appendTitle(ppvResultFile);
		writeTitleKspectrum = True;


def appendResultsToGraphFile(filepath, label, sensitivityResultFile, ppvResultFile):
	writeTitle(filepath, label, sensitivityResultFile, ppvResultFile);

	if filepath.find("dreme") != -1:
		suffix = "dreme";
	elif filepath.find("kspectrum") != -1:
		suffix = "kspectrum";
	else:
		suffix = "";

	with open(filepath, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			sensitivity = row[0].split(":")[1]
			ppv = row[1].split(":")[1]
			label = os.path.splitext(label)[0] + "_" + suffix;

			appendResult(sensitivityResultFile, label, sensitivity);
			appendResult(ppvResultFile, label, ppv);

def createTitle(label):
	tokens= label.split("_");
	seqLen = tokens[0]
	numSeq = tokens[1]

	title = "";
	title = "Number of Sequences: " + tokens[0] + ", Sequence Length: " + tokens[1];
	if len(tokens) == 5:
	 	alpha = tokens[2]
	 	signalPercent = tokens[3]
	 	signal = tokens[4]
	 	title = title + ", Signal Percent: " + signalPercent + ", Signal: " + signal
	else:
	 	signalPercent = tokens[2]
	 	signal = tokens[3]	
	 	title = title + ", Signal Percent: " + signalPercent + " Signal: " + signal	
	
	return title;

def getToolArray(labels, value):
	dremeDict = {}
	kspectrumDict = {}

	for idx, label in enumerate(labels):
		if label.find("dreme") != -1:
			newLabel = label.split("_")[1]
			dremeDict[newLabel] = value[idx]
		elif label.find("kspectrum") != -1:
			newLabel = label.split("_")[1]

			kspectrumDict[newLabel] = value[idx]

	print "DREME values: ", dremeDict.values()
	print "kspectrum values: ", kspectrumDict.values()
	return dremeDict, kspectrumDict;

def writeAndSavePlot(fileName, title, uniqueLabels, dremeValues, kspectrumValues):
	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	matplotlib.rcParams.update({'font.size': 8})
	
	fig = plt.figure()

	if fileName.find("sensitivity_graph") != -1:
		title = "Effect of sequence length on Sensitivity"
		plt.ylabel('Sensitivity');	
	else:
		title = "Effect of sequence length on PPV"
		plt.ylabel('PPV');	

	plt.title(title);
	line_up, = plt.plot(dremeValues, label="Dreme")
	line_down, = plt.plot(kspectrumValues, label = "k-spectrum")
	plt.legend([line_up, line_down], ['Dreme', 'k-spectrum'])

	plt.xticks(range(len(dremeValues)), list(uniqueLabels), rotation=90)
	figureName = os.path.splitext(fileName)[0] + ".png"
	plt.savefig(figureName);
	plt.close(fig)

def graphResults(fileName):

	labels = getColumn(fileName, 0)
	value = getColumn(fileName, 1)

	labels = labels[1:]
	value = value[1:]

	dremeDict, kspectrumDict = getToolArray(labels, value);
	uniqueLabels = createUniqueLabels(labels)

	dremeKeys = sorted(dremeDict)
	print "Type: ", type(dremeDict)
	print dremeDict
	title = createTitle(labels[0]);

	writeAndSavePlot(fileName, title, list(dremeKeys), list(dremeDict.values()), list(kspectrumDict.values()))


def splitall(path):
	allparts = []
	while 1:
		parts = os.path.split(path)
		if parts[0] == path:  # sentinel for absolute paths
			allparts.insert(0, parts[0])
			break
		elif parts[1] == path: # sentinel for relative paths
			allparts.insert(0, parts[1])
			break
		else:
			path = parts[0]
			allparts.insert(0, parts[1])
	return allparts

def sortResultFile(resultFile):
	f = open(resultFile,'r')
	lines = f.readlines()[1:]
	f.close()
	lines = sorted(lines)

	f = open(resultFile, "w");
	f.write("label, value\n");
	for line in lines:
		f.write(line)
	f.close();

def createUniqueLabels(labels):
	uniqueLabels = set();

	for label in labels:
		newLabel = label[:label.rfind("_")]
		uniqueLabels.add(newLabel)
	newLabels = [s.split("_")[1] for s in uniqueLabels]
		
	print newLabels;
	return newLabels

def getFileNameAndPrefixToGraph(filepath, subdir):

	fileNameToGraph =  os.path.splitext(os.path.basename(filepath))[0]


	fileNameToGraph = fileNameToGraph[7:]; # remove Signal_
	fileNameToGraph = os.path.splitext(fileNameToGraph)[0]

	allparts = splitall(subdir)
	graphPrefix = allparts[-2]	 #Get second to last directory.

	return fileNameToGraph, graphPrefix;

def parseDirectory(resultDir):

	global writeTitleDreme;
	global writeTitleKspectrum;

	#print resultDir;

	sensitivityResultFile = resultDir + "/sensitivity_graph.csv"
	ppvResultFile = resultDir + "/ppv_graph.csv"

	writeTitleKspectrum = writeTitleDreme = False;

	for subdir, dirs, files in os.walk(resultDir):
	
		for file in files:
			filepath = subdir + os.sep + file
			if filepath.endswith(".results"):
	 			fileNameToGraph, graphPrefix = getFileNameAndPrefixToGraph(filepath, subdir)
	 			appendResultsToGraphFile(filepath, fileNameToGraph, sensitivityResultFile, ppvResultFile)

	sortResultFile(sensitivityResultFile)
	sortResultFile(ppvResultFile)

	searchStr = resultDir + "/*graph.csv"
	for resultFile in glob.glob(searchStr):
		print "File to graph: ", resultFile;
		graphResults(resultFile)



if __name__ == "__main__":
	resultDir = sys.argv[1]
	
	parseDirectory(resultDir)
	#graphResults(dremeResultsFile);
	#graphResults(kspectrumResultsFile)
