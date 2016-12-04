import sys
import os
import csv
import glob

import numpy as np;
from scipy import stats
import math;
import collections;
from decimal import *
import datetime;
import generateGraphs;

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
		justFileName = os.path.basename(filepath)
		#split of the file name returns ['Signal_100_256_100_HuR_ParClip1', 'pwm', 'fa_kspectrum_25', 'results']
		suffix = justFileName.split(".")[2][3:];
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

def getToolArray(resultDict, index):
	dremeDict = {}
	kspectrumDict_25 = {}
	kspectrumDict_50 = {}
	kspectrumDict_100 = {}

	for key, value in resultDict.iteritems(): 
		if key.find("dreme") != -1:
			label = key.split("_")[index]
			dremeDict[int(label)] = value;
			#print label, value;

		elif key.find("kspectrum") != -1:
			label = key.split("_")[index]
			if key.rfind("kspectrum_10") != -1:
				kspectrumDict_25[int(label)] = value;
			elif key.rfind("kspectrum_15") != -1:
				kspectrumDict_50[int(label)] = value;
				#print kspectrumDict_50
			elif key.rfind("kspectrum_25") != -1:
				kspectrumDict_100[int(label)] = value;

	#print "100 dict: ", kspectrumDict_100
	#print "50 dict: ", kspectrumDict_50
	#print "25 dict: ", kspectrumDict_25

			
		
			#print label, value;

	#print "DREME values: ", dremeDict.values()
	#print "kspectrum values: ", kspectrumDict.values()
	dremeOd = collections.OrderedDict(sorted(dremeDict.items()))
	kspectrum25Od = collections.OrderedDict(sorted(kspectrumDict_25.items()))
	kspectrum50Od = collections.OrderedDict(sorted(kspectrumDict_50.items()))
	kspectrum100Od = collections.OrderedDict(sorted(kspectrumDict_100.items()))
	#print "DREME OD: ", dremeOd.values()
	#print "kspectrum OD: ", kspectrumOd.values()
	kspectrumOd_list = [kspectrum25Od, kspectrum50Od, kspectrum100Od]
	return dremeOd, kspectrumOd_list;

def writeAndSavePlot(fileName, title, xAxisTitle, index, labels, dremeMeanValues, kspectrumMeanValues, dremeErrorValues, kspectrumErrorValues):
	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	matplotlib.rcParams.update({'font.size': 8})
	
	fig = plt.figure()

	if fileName.find("sensitivity_graph") != -1:
		title = title + " on Sensitivity";
		plt.ylabel('Sensitivity');	
	else:
		title = title + " on PPV";
		plt.ylabel('PPV');	

	plt.xlabel(xAxisTitle);
	plt.title(title);

	#print "PLOTTING DREME: ", labels, dremeMeanValues, dremeErrorValues
	#print "PLOTTING kspectrum: ", labels, kspectrumMeanValues, kspectrumErrorValues
	
	#plt.xticks(range(len(labels)), list(labels), rotation=90)
	eb1 = plt.errorbar(labels, dremeMeanValues, dremeErrorValues, fmt='', color='b')
	eb2 = plt.errorbar(labels, kspectrumMeanValues[0], kspectrumErrorValues[0], fmt='', color='g')
	eb3 = plt.errorbar(labels, kspectrumMeanValues[1], kspectrumErrorValues[1], fmt='', color='r')
	eb4 = plt.errorbar(labels, kspectrumMeanValues[2], kspectrumErrorValues[2], fmt='', color='y')

	plt.legend([eb1, eb2, eb3, eb4, eb5], ['DREME', 'k-spectrum-25', 'k-spectrum-50', 'k-spectrum-100'])
	plt.savefig(fileName);
	plt.close(fig)

def ComputeMeanAndStdError(resultDict):
	meanDict = {}
	for key, value in resultDict.iteritems():
		results = [float(i) for i in resultDict[key]]
		meanValue = np.mean(results);
		errorValue = stats.sem(results);
		#print "MEAN, ERROR: ", str(meanValue), str(errorValue)
		meanDict[key] = [round(meanValue, 4), round(errorValue, 4)]
	#print meanDict
	return meanDict;	

def GetStdDeviationDict(resultDict):
	stdDict = {}
	for key, value in resultDict.iteritems():
		results = [float(i) for i in resultDict[key]]
		stdValue = np.std(results);
		stdDict[key] = stdValue;
	return stdDict;

def graphResults(fileName, resultDict, title, xAxisTitle, index):


	print resultDict
	stdDict = GetStdDeviationDict(resultDict)
	meanDict = ComputeMeanAndStdError(resultDict)
	dremeDict, kspectrumDictList = getToolArray(meanDict, index);

	#print dremeDict

	labels = dremeDict.keys()

	dremeMeanValues = [x[0] for x in dremeDict.values()]
	dremeErrorValues = [x[1] for x in dremeDict.values()]

	print "DREME  VALUES: ", dremeMeanValues, dremeErrorValues;
	kspectrum25MeanValues = [x[0] for x in kspectrumDictList[0].values()]
	kspectrum25ErrorValues = [x[1] for x in kspectrumDictList[0].values()]	

	kspectrum50MeanValues = [x[0] for x in kspectrumDictList[1].values()]
	kspectrum50ErrorValues = [x[1] for x in kspectrumDictList[1].values()]	

	kspectrum100MeanValues = [x[0] for x in kspectrumDictList[1].values()]
	kspectrum100ErrorValues = [x[1] for x in kspectrumDictList[1].values()]	


	kspectrumMeanValues = [kspectrum25MeanValues, kspectrum50MeanValues, kspectrum100ErrorValues]
	kspectrumErrorValues = [kspectrum25ErrorValues, kspectrum50ErrorValues, kspectrum100ErrorValues]
	print "kspectrum  VALUES: ", kspectrumMeanValues, kspectrumErrorValues;

	writeAndSavePlot(fileName, title, xAxisTitle, index, labels, dremeMeanValues, kspectrumMeanValues, dremeErrorValues, kspectrumErrorValues)



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
		
	#print newLabels;
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
	 			fileNameToGraph, graphPrefix = getFileNameAndPrefixToGraph(filepath, subdir);
	 			appendResultsToGraphFile(filepath, fileNameToGraph, sensitivityResultFile, ppvResultFile)

	sortResultFile(sensitivityResultFile)
	sortResultFile(ppvResultFile)

def ReadFileAndCreateDict(resultFile):
	resultDict = {};
	with open(resultFile, mode='r') as infile:
		next(infile)
		reader = csv.reader(infile)
		for rows in reader:
			label = rows[0];
			value = rows[1];
			#print label, value;	
			if label in resultDict:
				resultDict[label].append(value)
			else:
				resultDict[label] = [value];
	return resultDict;


def GetMeanAndStdErrorValues(resultDict, index):
	stdDict = GetStdDeviationDict(resultDict)
	meanDict = ComputeMeanAndStdError(resultDict)
	dremeDict, kspectrumDict = getToolArray(meanDict, index);

	dremeMeanValues = [x[0] for x in dremeDict.values()]
	dremeErrorValues = [x[1] for x in dremeDict.values()]

	#print "DREME  VALUES: ", dremeMeanValues, dremeErrorValues;
	kspectrum25MeanValues = [x[0] for x in kspectrumDict[0].values()]
	kspectrum25ErrorValues = [x[1] for x in kspectrumDict[0].values()]		

	kspectrum50MeanValues = [x[0] for x in kspectrumDict[1].values()]
	kspectrum50ErrorValues = [x[1] for x in kspectrumDict[1].values()]		

	kspectrum100MeanValues = [x[0] for x in kspectrumDict[2].values()]
	kspectrum100ErrorValues = [x[1] for x in kspectrumDict[2].values()]		

	labels = dremeDict.keys();

	kspectrumMeanValues = [kspectrum25MeanValues, kspectrum50MeanValues, kspectrum100MeanValues]
	kspectrumErrorValues = [kspectrum25ErrorValues, kspectrum50ErrorValues, kspectrum100ErrorValues]

	#print kspectrumMeanValues, kspectrumErrorValues
	return dremeMeanValues, dremeErrorValues, kspectrumMeanValues, kspectrumErrorValues, labels;

def GetDictForResultFile(searchStr):

	for resultFile in glob.glob(searchStr):
		#print "File to graph: ", resultFile;
		result_dict = ReadFileAndCreateDict(resultFile)
		#print "Result Dict: ", result_dict

	return result_dict;

def createDataToGraph(resultDir):
	searchStr = resultDir + "/ppv_graph.csv"
	ppv_dict = GetDictForResultFile(searchStr)

	searchStr = resultDir + "/sensitivity_graph.csv"
	sensitivity_dict = GetDictForResultFile(searchStr)

	#print sensitivity_dict, ppv_dict

	return sensitivity_dict, ppv_dict;

def PlotSingleGraph(sensitivity_dict, ppv_dict, title, xAxisTitle, index, graphNamePrefix):
	graphName = graphNamePrefix + ".png";
	generateGraphs.PlotSensitivityAndPPVGraphs(sensitivity_dict, ppv_dict, title, xAxisTitle, index, graphName);


def PlotSeparateGraphs(sensitivity_dict, ppv_dict, title, xAxisTitle, index, graphNamePrefix):
	graphName = graphNamePrefix + "_" + "sensitivity_graph.png"
	graphResults(graphName,sensitivity_dict, title, xAxisTitle, index);

	graphName = graphNamePrefix + "_" + "ppv_graph.png"
	graphResults(graphName, ppv_dict, title, xAxisTitle, index)


def GraphResults(resultDir, title, xAxisTitle, index):
	graphNamePrefix = os.path.basename(os.getcwd())
	parseDirectory(resultDir)
	sensitivity_dict, ppv_dict = createDataToGraph(resultDir)

	PlotSingleGraph(sensitivity_dict, ppv_dict, title, xAxisTitle, index, graphNamePrefix)

	#PlotSeparateGraphs(sensitivity_dict, ppv_dict, title, xAxisTitle, index, graphNamePrefix)




if __name__ == "__main__":
	resultDir = sys.argv[1]
	title = sys.argv[2]
	xAxisTitle = sys.argv[3]
	index = int(sys.argv[4])
	
	GraphResults(resultDir, title, xAxisTitle, index);

	#print resultDir;

	# graphNamePrefix = os.path.basename(os.getcwd())
	# parseDirectory(resultDir)
	# sensitivity_dict, ppv_dict = createDataToGraph(resultDir)

	# #dateStr = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
	# dateStr = datetime.datetime.now().strftime('%Y-%m-%d')
	# graphName = dateStr + "_" + graphNamePrefix + "_" + "sensitivity_graph.png"
	# graphResults(graphName,sensitivity_dict, title, index);

	# graphName = dateStr + "_" + graphNamePrefix + "_" + "ppv_graph.png"
	# graphResults(graphName, ppv_dict, title, index)
