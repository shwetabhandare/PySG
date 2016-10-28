import sys
import os
import generateGraphs
import SeqGenUtils;
import gcContent;
import numpy as np;
from scipy import stats
import csv;

def writeDictToFile(gcContentMeanStd, resultFileName):
	with open(resultFileName, 'wb') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in gcContentMeanStd.items():
			writer.writerow([key, value])	

def parseSubDirectories(resultDir, level=1):
	gcContentMap = dict();
	resultFileName = resultDir + resultDir[:-1] + "_GC_Content.out"
	print "Result FileName: ", resultFileName
	for signalFile in SeqGenUtils.findFiles(resultDir, "Signal*.fa"):
		print "Signal File: ", signalFile;
		expt_name = os.path.dirname(signalFile).split("/")[2]
		gcContentValue, atContentValue = gcContent.getNucleotideComposition(signalFile)
		print gcContentValue, atContentValue;
		if expt_name in gcContentMap.keys():
			gcContentMap[str(expt_name)].append(gcContentValue);
		else:
			gcContentMap[str(expt_name)] = [gcContentValue];

	gcContentMeanStd = dict();
	for key, value in gcContentMap.iteritems():
		meanValue = np.mean(value);
		varianceValue = np.var(value);
		gcContentMeanStd[key] = [meanValue, varianceValue]

	writeDictToFile(gcContentMeanStd, resultFileName);

if __name__ == "__main__":
	resultDir = sys.argv[1]
	parseSubDirectories(resultDir)
