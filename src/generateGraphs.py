import sys;
import os, glob;
import parseResults;
import shutil;
import datetime;
import numpy as np;
import shuffle_utils;
import Distribution_Utils;

threeUtrValues = list();
threeUtrErrors = list();
threeUtrNGrams = dict();

def isDirResultDir(currentDir):
	validSubStrings = ["Percent", "Length", "Alpha", "Dirichlet", "Shuffle"]
	if any(x in currentDir for x in validSubStrings):
		return True

	return False; 

def createTitleFromDirName(dirName):
	#example: HuR_Test_SignalPercent_Dirichlet
	fileNameTokens = os.path.basename(dirName).split("_");
	if fileNameTokens[2] == "SeqLength":
		index = 1;
		xAxisTitle = "Sequence Length";
	elif fileNameTokens[2] == "SignalPercent":
		index = 2;
		xAxisTitle = "Signal Percent";
	else:
		xAxisTitle = fileNameTokens[2];
		index = 3;

	title = fileNameTokens[0] + ": Effect of " + xAxisTitle;
	return xAxisTitle, title, index;

def CopyGraphToTargetDir():
	graphDir = "/projects/bhandare/workspace/PySG/src/tmp_results";
	for file in glob.glob("*.png"):
		#print "Graph File: ", file;
		dateStr = datetime.datetime.now().strftime('%Y-%m-%d')
		targetDir = graphDir + "/" + dateStr;
		print "Copying ", file, " to ", targetDir;
		if not os.path.exists(targetDir):
			os.makedirs(targetDir)
		shutil.copy2(file, targetDir)

def ParseResultsAndGenerateGraph(path):
	#print "Path to Graph: ", path;
	xAxisTitle, title, index = createTitleFromDirName(path)
	savedDir = os.getcwd();
	mycwd = os.chdir(path)
	parseResults.GraphResults(".", title, xAxisTitle, index);
	CopyGraphToTargetDir();
	os.chdir(savedDir);

def GetSubDirCount(path, map = {}):
	count = 0
	for f in os.listdir(path):
		child = os.path.join(path, f)
		if os.path.isdir(child):
			count = count + 1 # unless include self
		else:
			continue;
	map[path] = count
	return count


def distribution_plot(ax, title, meanValues, errorValues, labels):
	global threeUtrValues, threeUtrErrors;

	#ax.set_ylabel(yLabel)
	ax.set_title(title);
	xValues = np.arange(len(meanValues))

	eb1 = ax.errorbar(xValues, meanValues, yerr=errorValues, fmt='', color='b')
	eb2 = ax.errorbar(xValues, threeUtrValues, yerr=threeUtrErrors, xerr=None, fmt='', color='g')

	ax.set_xticks(xValues)
	ax.set_xticklabels(labels, fontdict=None, minor=False, size='small');

	return eb1, eb2;

def sensitivity_ppv_plot(ax, title, dremeMeanValues, dremeErrorValues, kspectrumMeanValues, kspectrumErrorValues, labels, xLabel, yLabel):

	ax.set_ylabel(yLabel)
	#ax.set_xlabel(xLabel)
	ax.set_title(title);
	xValues = np.arange(len(dremeMeanValues))

	eb1 = ax.errorbar(labels, dremeMeanValues, dremeErrorValues, fmt='', color='b')
	#print kspectrumMeanValues[0], kspectrumErrorValues[0]
	eb2 = ax.errorbar(labels, kspectrumMeanValues[0], kspectrumErrorValues[0], fmt='go', color='g')
	#print kspectrumMeanValues[1], kspectrumErrorValues[1]
	eb3 = ax.errorbar(labels, kspectrumMeanValues[1], kspectrumErrorValues[1], fmt='k^', color='r')
	eb4 = ax.errorbar(labels, kspectrumMeanValues[2], kspectrumErrorValues[2], fmt='r--', color='y')
	ax.set_xticks(labels)
	ax.set_xticklabels(labels, fontdict=None, minor=False, size='small');

	return eb1, eb2, eb3, eb4;

def ngram_plot(title, posMeanValues, posErrorValues, negMeanValues, negErrorValues, threeUtrValues, labels, xLabel, yLabel, graphFileName):

	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	matplotlib.rcParams.update({'font.size': 8})

	threeUtrErrorValues = np.zeros(len(posErrorValues))	
	fig = plt.figure()

	plt.title(title);

	xValues = np.arange(len(posMeanValues))
	print "X-Values: ", xValues, len(xValues);

	fig, (ax) = plt.subplots(nrows=1, ncols=1)
	eb1 = ax.errorbar(xValues, posMeanValues, yerr=posErrorValues, fmt='', color='b', label="Positive File")
	eb2 = ax.errorbar(xValues, negMeanValues, yerr=negErrorValues, fmt='o', color='r', label="Negative File")
	eb3  = ax.errorbar(xValues, threeUtrValues, yerr=threeUtrErrorValues, fmt='', color='g', label="3'UTR")

	ax.legend()
	ax.set_ylabel(yLabel)
	ax.set_xlabel(xLabel)
	fig.suptitle("Comparison of N-gram distribution: 3'UTR vs Generated Files")


	plt.xticks(xValues, labels, rotation='vertical')

	plt.savefig(graphFileName);
	plt.close(fig)

	return eb1, eb2;

def PlotNGrams(posDict, negDict, threeUtrDict, graphTitle, graphFileName):
	
	import matplotlib
	matplotlib.use('Agg')

	from matplotlib import pyplot as plt	

	three_utr_values = threeUtrDict.values();
	labels = threeUtrDict.keys();

	posMeanValues = [x[0] for x in posDict.values()]
	posErrorValues = [x[1] for x in posDict.values()]	

	negMeanValues = [x[0] for x in negDict.values()]
	negErrorValues = [x[1] for x in negDict.values()]	

	ngram_plot(graphTitle, posMeanValues, posErrorValues, negMeanValues, negErrorValues, three_utr_values, labels, "N-Grams", "Frequency", graphFileName)	




def PlotSensitivityAndPPVGraphs(sensitivityDict, ppvDict, graphTitle, xAxisTitle, index, graphFileName):
	
	sensDremeMeanValues, sensDremeErrorValues, sensKspectrumMeanValues, sensKspectrumErrorValues, labels = parseResults.GetMeanAndStdErrorValues(sensitivityDict, index);
	ppvDremeMeanValues, ppvDremeErrorValues, ppvKspectrumMeanValues, ppvKspectrumErrorValues, labels = parseResults.GetMeanAndStdErrorValues(ppvDict, index);

	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	matplotlib.rcParams.update({'font.size': 8})

	#Two subplots, the axes array is 1-d
	fig, (ax1, ax2) = plt.subplots(nrows=2)
	eb1, eb2, eb3, eb4 = sensitivity_ppv_plot(ax1, graphTitle + " on Sensitivity", sensDremeMeanValues, sensDremeErrorValues, sensKspectrumMeanValues, sensKspectrumErrorValues, labels, xAxisTitle, "Sensitivity");
	eb1, eb2, eb3, eb4= sensitivity_ppv_plot(ax2, graphTitle + " on PPV", ppvDremeMeanValues, ppvDremeErrorValues, ppvKspectrumMeanValues, ppvKspectrumErrorValues, labels, xAxisTitle, "PPV");	

	ax2.set_xlabel(xAxisTitle)

	plt.figlegend((eb1, eb2, eb3, eb4), ("DREME", "k-spectrum-25", "k-spectrum-50", "k-spectrum-100"), loc = 'lower right');
	plt.savefig(graphFileName);
	plt.close(fig)

def PlotGraphs(posDistList, negDiNucDist, graphFileName):
	global threeUtrValues, threeUtrErrors;

	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	matplotlib.rcParams.update({'font.size': 8})

	posMeanValues = [x[0] for x in posDistList.values()]
	posErrorValues = [x[1] for x in posDistList.values()]	

	negMeanValues = [x[0] for x in negDiNucDist.values()]
	negErrorValues = [x[1] for x in negDiNucDist.values()]	

	#Two subplots, the axes array is 1-d
	fig, (ax1, ax2) = plt.subplots(nrows=2)
	fig.suptitle("Di-nucleotide distribution: 3'UTR vs Generated Files", fontsize=8)
	labels = posDistList.keys();
	eb1, eb2 = distribution_plot(ax1, "Positive Set", posMeanValues, posErrorValues, labels);
	eb1, eb2 = distribution_plot(ax2, "Negative Set", negMeanValues, negErrorValues, labels);

	plt.figlegend((eb1, eb2), ("Generated", "3'UTR"), loc = 'lower right');
	plt.savefig(graphFileName);
	plt.close(fig)



def GenerateDiNucleotideGraphs(resultDir, currDir, posDiNucDist, negDiNucDist):

	targetDir = resultDir + "/" + "dinucleotide_distribution"
	if not os.path.exists(targetDir):
		os.makedirs(targetDir);

	graphFileName = targetDir + "/" + currDir + ".png"
	PlotGraphs(posDiNucDist, negDiNucDist, graphFileName)
	


def computeDiNucleotideDistribution(resultDir, level=1):
	posDistList = list()
	negDistList = list();
	global threeUtrValues, threeUtrErrors;

	threeUtrValues, threeUtrErrors = Distribution_Utils.Compute3UtrDistibution();

	num_sep = resultDir.count(os.path.sep)
	for root, dirs, files in os.walk(resultDir):
		num_sep_this = root.count(os.path.sep)
		if num_sep + level > num_sep_this:
			for dir in dirs:
				if isDirResultDir(dir):
					posDiNucDist, negDiNucDist = Distribution_Utils.GetDiNucleotideDistribution(root + "/" + dir)
					GenerateDiNucleotideGraphs(root, dir, posDiNucDist, negDiNucDist);
						

def parseSubDirectories(resultDir, level=1):
	num_sep = resultDir.count(os.path.sep)
	for root, dirs, files in os.walk(resultDir):
		num_sep_this = root.count(os.path.sep)
		if num_sep + level > num_sep_this:
			for dir in dirs:
				if isDirResultDir(dir):
					numSubDirs = GetSubDirCount(root + "/" + dir)
					print "Found Result Directory: ", dir, ", num sub-dir: ", str(numSubDirs)
					if numSubDirs == 5:
						ParseResultsAndGenerateGraph(root + "/" + dir)
					else:
						print "Look at directory: ", dir;
						continue;

if __name__ == "__main__":
	resultDir = sys.argv[1]
	parseSubDirectories(resultDir)
	#computeDiNucleotideDistribution(resultDir)
