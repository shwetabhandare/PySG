import sys
import os
import csv

def getResultsFromFile(filepath):
	with open(filepath, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			sensitivity = row[0].split(":")[1]
			ppv = row[1].split(":")[1]

			print sensitivity
			print ppv


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
	            getResultsFromFile(filepath)


if __name__ == "__main__":
	resultDir = sys.argv[1]
	parseDirectory(resultDir)