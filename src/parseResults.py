import sys
import os
import csv
import matplotlib
from matplotlib import pyplot as plt
import numpy as np;
import math;

writeTitle = False;


def simple():
    x_data = np.linspace(0., 100., 1000)
    y_data = ["Signal_100_500_90_100_Cgd2_3490.pwm.fa_dreme", "Signal_100_500_90_100_Gal4.pwm.fa_dreme"]
    x_data = [0.54, 0.808]
    for counter, x in enumerate(x_data):
        y = x_data[counter]
        matplotlib.pyplot.scatter(x, y)


    axes = matplotlib.pyplot.gca()
    axes.set_xlabel('x')
    axes.set_ylabel('y')
    matplotlib.pyplot.show()

def getColumn(filename, column):
	results = csv.reader(open(filename), delimiter=",")
	return [result[column] for result in results]

def appendTitle():
	with open('resultsToGraph.csv', 'wb') as csvfile:	
		spamwriter = csv.writer(csvfile, delimiter=',');
		spamwriter.writerow(["label", "sensitivity", "ppv"]);

def appendResult(label, sensitivity, ppv):
	with open('resultsToGraph.csv', 'a') as csvfile:	
		spamwriter = csv.writer(csvfile, delimiter=',');
		spamwriter.writerow([label, sensitivity, ppv])

def getResultsFromFile(filepath, label):
	global writeTitle;

	if writeTitle == False:
		appendTitle();
		writeTitle = True;


	with open(filepath, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			sensitivity = row[0].split(":")[1]
			ppv = row[1].split(":")[1]
			print label, sensitivity, ppv
			appendResult(label, sensitivity, ppv)

def graphResults():

	#data1=np.genfromtxt('data1.csv', skip_header=1) 
	#plt.plot(data1)


	label = getColumn("resultsToGraph.csv", 0)
	sensitivity = getColumn("resultsToGraph.csv", 1)
	label = label[1:]
	sensitivity = sensitivity[1:]
	print label, sensitivity
	ppv = getColumn("resultsToGraph.csv", 2)
	ppv = ppv[1:]

	axes = matplotlib.pyplot.gca()
	axes.set_xlabel('x')
	axes.set_ylabel('y')

	#plt.figure("sensitivity/ppv")
	#plt.xlabel("filename")
	#plt.ylabel("sensitivity/ppv")
	plt.plot(label, sensitivity);
	plt.show();


def parseDirectory(resultDir):

	for subdir, dirs, files in os.walk(resultDir):
	    for file in files:
	        #print os.path.join(subdir, file)
	        filepath = subdir + os.sep + file

	        if filepath.endswith(".results"):
	            #print (filepath)
	            fileNameWithoutPath = os.path.basename(filepath);
	            fileNameToGraph =  os.path.splitext(fileNameWithoutPath)[0]
	            print fileNameToGraph
	            getResultsFromFile(filepath, fileNameToGraph)


if __name__ == "__main__":
	resultDir = sys.argv[1]
	parseDirectory(resultDir)
	graphResults();
	#simple();