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
			#print "Label ", label;
			appendResult(fileName, label, sensitivity, ppv)

def createTitle(label):
	print "Label: ", label;
	x = label.split("_");
	seqLen = x[0]
	numSeq = x[1]

	title = "";

	title = "Number of Sequences: ", x[0], ", Sequence Length: ", x[1];
	if len(x) == 5:
		alpha = x[2]
		signalPercent = x[3]
		signal = x[4]
		title = title + ", Signal Percent: ", signalPercent, ", Signal: ", signal
	else:
		signalPercent = x[2]
		signal = x[3]	
		title = title + ", Signal Percent: ", signalPercent, ", Signal: ", signal	

	
	return title;

def graphResults(fileName):

	#data1=np.genfromtxt('data1.csv', skip_header=1) 
	#plt.plot(data1)


	labels = getColumn(fileName, 0)
	sensitivity = getColumn(fileName, 1)
	ppv = getColumn(fileName, 2)

	labels = labels[1:]
	title = createTitle(labels[0]);
	#print "before: ", labels;
	labels[:] = [l.split("_")[-1] for l in labels]
	#print "after: ", labels;

	sensitivity = sensitivity[1:]
	ppv = ppv[1:]

	#print "Title: ", title;
	#print labels
	#print sensitivity
	#print ppv
	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	#matplotlib.rcParams.update({'font.size': 8})
	
	#plt.ion()
	fig = plt.figure()
	#yerr = 0.1 + 0.2*np.sqrt(sensitivity)
	#yerr = 0.1 + 0.2*np.sqrt(sensitivity)

	plt.ylabel('Sensitivity/PPV');
	plt.title(title);
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
	graphPrefix = allparts[-2]	 #Get second to last directory.

	return fileNameToGraph, graphPrefix;

def parseDirectory(resultDir):

	global writeTitleDreme;
	global writeTitleKspectrum;

	for subdir, dirs, files in os.walk(resultDir):
		
		# if subdir != resultDir:
		# 	print "Sub dir: ", subdir;

		#print "Files: ", files;
		if subdir == resultDir or subdir.find("Signal") == -1:
			print "Sub-Dir: ", subdir;
			writeTitleDreme=False;
			writeTitleKspectrum=False;

			# if graphFlag:
			# 	print "Graph Flag: ", graphFlag;
			# 	if os.path.isfile(dremeResultsFile):  
			# 		graphResults(dremeResultsFile)
			# 	if os.path.isfile(kspectrumResultsFile): 
			# 		graphResults(kspectrumResultsFile)
			# 	graphFlag = False;


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


	for resultFile in glob.glob("*graph*.csv"):
		print "File to graph: ", resultFile;
		graphResults(resultFile)



if __name__ == "__main__":
	resultDir = sys.argv[1]
	parseDirectory(resultDir)
	#graphResults(dremeResultsFile);
	#graphResults(kspectrumResultsFile)
