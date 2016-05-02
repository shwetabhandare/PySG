import sys
import os
import csv

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
		spamwriter.writerow(["label", "sensitivity", "ppv"]);

def appendResult(fileName, label, sensitivity, ppv):
	with open(fileName, 'a') as csvfile:	
		spamwriter = csv.writer(csvfile, delimiter=',');
		spamwriter.writerow([label, sensitivity, ppv])

def getResultFile(filepath, label, dremeResultsFile, kspectrumResultsFile):
	global writeTitleDreme;
	global writeTitleKspectrum;

	if label.find("dreme") != -1:
		fileName = dremeResultsFile;
		if writeTitleDreme == False:
			appendTitle(fileName);
			writeTitleDreme = True;
	elif label.find("kspectrum") != -1:
		fileName = kspectrumResultsFile;
		if writeTitleKspectrum == False:
			appendTitle(fileName);
			writeTitleKspectrum = True;
	else:
		fileName = "resultsToGraph.csv"

	return fileName;


def appendResultsToGraphFile(filepath, label, dremeResultsFile, kspectrumResultsFile):
	fileName = getResultFile(filepath, label, dremeResultsFile, kspectrumResultsFile);

	with open(filepath, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			sensitivity = row[0].split(":")[1]
			ppv = row[1].split(":")[1]
			label = os.path.splitext(label)[0];
			appendResult(fileName, label, sensitivity, ppv)

def graphResults(fileName):

	#data1=np.genfromtxt('data1.csv', skip_header=1) 
	#plt.plot(data1)


	labels = getColumn(fileName, 0)
	sensitivity = getColumn(fileName, 1)
	ppv = getColumn(fileName, 2)

	labels = labels[1:]
	sensitivity = sensitivity[1:]
	ppv = ppv[1:]

	print labels
	print sensitivity
	print ppv
	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	matplotlib.rcParams.update({'font.size': 18})
	
	#plt.ion()
	fig = plt.figure()
	#fig = plt.figure(figsize=(20,30))
	plt.plot(sensitivity)
	plt.plot(ppv)
	plt.xticks(range(len(sensitivity)), labels, rotation=90)
	figureName = os.path.splitext(fileName)[0] + ".png"
	plt.savefig(figureName);
	plt.close(fig)

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

def getFileNameAndPrefixToGraph(filepath, subdir):
	fileNameToGraph =  os.path.splitext(os.path.basename(filepath))[0]
	fileNameToGraph = fileNameToGraph[7:]; # remove Signal_

	allparts = splitall(subdir)
	graphPrefix = allparts[-2]	

	return fileNameToGraph, graphPrefix;

def parseDirectory(resultDir):

	for subdir, dirs, files in os.walk(resultDir):
		graphFlag = False;
		# if subdir != resultDir:
		# 	print "Sub dir: ", subdir;

		for file in files:
			#print os.path.join(subdir, file)
			filepath = subdir + os.sep + file
			if filepath.endswith(".results"):
				graphFlag = True;
				print "Result file: ", filepath;
				fileNameToGraph, graphPrefix = getFileNameAndPrefixToGraph(filepath, subdir)
				dremeResultsFile = graphPrefix + "_dreme_graph.csv";
				kspectrumResultsFile  = graphPrefix + "_kspectrum_graph.csv"
				appendResultsToGraphFile(filepath, fileNameToGraph, dremeResultsFile, kspectrumResultsFile)

		if graphFlag:
			if os.path.isfile(dremeResultsFile):  
				graphResults(dremeResultsFile)
			if os.path.isfile(kspectrumResultsFile): 
				graphResults(kspectrumResultsFile)
			graphFlag = False;

if __name__ == "__main__":
	resultDir = sys.argv[1]
	parseDirectory(resultDir)
	#graphResults(dremeResultsFile);
	#graphResults(kspectrumResultsFile)
