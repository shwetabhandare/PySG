import sys
import os
import csv
import matplotlib
from matplotlib import pyplot as plt
import numpy as np;
import math;

writeTitleDreme = False;
writeTitleKspectrum = False;

dremeResultsFile = "dremeResultsToGraph.csv"
kspectrumResultsFile="kspectrumResultsToGraph.csv"


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

def getResultFile(filepath, label):
	global writeTitleDreme;
	global writeTitleKspectrum;
	global dremeResultsFile;
	global kspectrumResultsFile;

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


def getResultsFromFile(filepath, label):


	fileName = getResultFile(filepath, label);

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
	
	matplotlib.rcParams.update({'font.size': 12})

	plt.ion()
	#plt.figure()
	plt.figure(figsize=(20,30))

	plt.plot(sensitivity)
	plt.plot(ppv)
	plt.xticks(range(len(sensitivity)), labels, rotation=90)
	figureName = os.path.splitext(fileName)[0] + ".png"
	plt.savefig(figureName);


def parseDirectory(resultDir):

	for subdir, dirs, files in os.walk(resultDir):
	    for file in files:
	        #print os.path.join(subdir, file)
	        filepath = subdir + os.sep + file
	        if filepath.endswith(".results"):
	            #print "RESULTS FILE: ", filepath
	            fileNameWithoutPath = os.path.basename(filepath);
	            fileNameToGraph =  os.path.splitext(fileNameWithoutPath)[0]
	            fileNameToGraph = fileNameToGraph[7:]; # remove Signal_
	            getResultsFromFile(filepath, fileNameToGraph)


if __name__ == "__main__":
	global dremeResultsFile;
	global kspectrumResultsFile;
	resultDir = sys.argv[1]
	parseDirectory(resultDir)
	graphResults(dremeResultsFile);
	graphResults(kspectrumResultsFile)
