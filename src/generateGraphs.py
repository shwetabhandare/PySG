import sys;
import os, glob;
import parseResults;
import shutil;
import datetime;
import SeqGenUtils;
import numpy as np;
import shuffle_utils;
from scipy import stats;

threeUtrValues = list();
threeUtrErrors = list();

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
	print "Title: ", title;
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

def GetDistribListForDirectory(directory, filesToFind):
	distribList = list();
	for seqFile in SeqGenUtils.findFiles(directory, filesToFind):
		seqs, gc_list, fg_lengths = shuffle_utils.get_seqs(seqFile)
		dinuc_distrib = shuffle_utils.compute_dinuc_distrib(seqs, True)			
		distribList.append(dinuc_distrib)		
	return distribList;

def ComputeIndividualDistributions(posDistList):

	aa_values = [x['AA'] for x in posDistList]
	ac_values = [x['AC'] for x in posDistList]
	gt_values = [x['GT'] for x in posDistList]
	ag_values = [x['AG'] for x in posDistList]
	cc_values = [x['CC'] for x in posDistList]
	ca_values = [x['CA'] for x in posDistList]
	cg_values = [x['CG'] for x in posDistList]
	tt_values = [x['TT'] for x in posDistList]
	gg_values = [x['GG'] for x in posDistList]
	gc_values = [x['GC'] for x in posDistList]
	at_values = [x['AT'] for x in posDistList]
	ga_values = [x['GA'] for x in posDistList]
	gt_values = [x['GT'] for x in posDistList]
	ta_values = [x['TA'] for x in posDistList]
	tc_values = [x['TC'] for x in posDistList]
	ct_values = [x['CT'] for x in posDistList]

	diNucDist = dict();

	diNucDist['AA'] = [round(np.mean(aa_values), 4), round(stats.sem(aa_values), 4)]
	diNucDist['AC'] = [round(np.mean(ac_values), 4),  round(stats.sem(ac_values), 4)]
	diNucDist['GT'] = [round(np.mean(gt_values), 4),  round(stats.sem(gt_values), 4)]
	diNucDist['AG'] = [round(np.mean(ag_values), 4), round(stats.sem(ag_values), 4)]
	diNucDist['CC'] = [round(np.mean(cc_values), 4),  round(stats.sem(cc_values), 4)]
	diNucDist['CA'] = [round(np.mean(ca_values), 4),  round(stats.sem(ca_values), 4)]
	diNucDist['CG'] = [round(np.mean(cg_values), 4),  round(stats.sem(cg_values), 4)]
	diNucDist['TT'] = [round(np.mean(tt_values), 4),  round(stats.sem(tt_values), 4)]
	diNucDist['GG'] = [round(np.mean(gg_values), 4),  round(stats.sem(gg_values), 4)]
	diNucDist['GC'] = [round(np.mean(gc_values), 4),  round(stats.sem(gc_values), 4)]
	diNucDist['AT'] = [round(np.mean(at_values), 4),  round(stats.sem(at_values), 4)]
	diNucDist['GA'] = [round(np.mean(ga_values), 4),  round(stats.sem(ga_values), 4)]
	diNucDist['GT'] = [round(np.mean(gt_values), 4),  round(stats.sem(gt_values), 4)]
	diNucDist['TA'] = [round(np.mean(ta_values), 4),  round(stats.sem(ta_values), 4)]
	diNucDist['TC'] = [round(np.mean(tc_values), 4),  round(stats.sem(tc_values), 4)]
	diNucDist['CT'] = [round(np.mean(ct_values), 4),  round(stats.sem(ct_values), 4)]


	return diNucDist;

def Compute3UtrDistibution():
	global threeUtrValues, threeUtrErrors;

	three_utr_file_dir = "/projects/bhandare/workspace/scripts/NegFileCreator/";
	three_utr_file = "3UTR_transcripts_Human.txt"

	threeUtrDist = GetDistributionForFiles(three_utr_file_dir, three_utr_file);
	threeUtrValues = [x[0] for x in threeUtrDist.values()]
	threeUtrErrors = [x[1] for x in threeUtrDist.values()]	

	return threeUtrValues, threeUtrErrors;

def distribution_plot(ax, title, meanValues, errorValues, labels):
	global threeUtrValues, threeUtrErrors;

	ax.set_ylabel(yLabel)
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
	eb2 = ax.errorbar(labels, kspectrumMeanValues, kspectrumErrorValues, fmt='', color='g')

	ax.set_xticks(labels)
	ax.set_xticklabels(labels, fontdict=None, minor=False, size='small');

	return eb1, eb2;



def PlotSensitivityAndPPVGraphs(sensitivityDict, ppvDict, graphTitle, xAxisTitle, index, graphFileName):
	
	sensDremeMeanValues, sensDremeErrorValues, sensKspectrumMeanValues, sensKspectrumErrorValues, labels = parseResults.GetMeanAndStdErrorValues(sensitivityDict, index);
	ppvDremeMeanValues, ppvDremeErrorValues, ppvKspectrumMeanValues, ppvKspectrumErrorValues, labels = parseResults.GetMeanAndStdErrorValues(ppvDict, index);

	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	
	matplotlib.rcParams.update({'font.size': 8})

	#Two subplots, the axes array is 1-d
	fig, (ax1, ax2) = plt.subplots(nrows=2)
	eb1, eb2 = sensitivity_ppv_plot(ax1, graphTitle + " on Sensitivity", sensDremeMeanValues, sensDremeErrorValues, sensKspectrumMeanValues, sensKspectrumErrorValues, labels, xAxisTitle, "Sensitivity");
	eb1, eb2 = sensitivity_ppv_plot(ax2, graphTitle + " on PPV", ppvDremeMeanValues, ppvDremeErrorValues, ppvKspectrumMeanValues, ppvKspectrumErrorValues, labels, xAxisTitle, "PPV");	

	ax2.set_xlabel(xAxisTitle)

	plt.figlegend((eb1, eb2), ("DREME", "k-spectrum"), loc = 'lower right');
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


def GetDistributionForFiles(resultDir, searchRe):

	posDistList = GetDistribListForDirectory(resultDir, searchRe)
	posDiNucDist = ComputeIndividualDistributions(posDistList)
	return posDiNucDist;

def GetDiNucleotideDistribution(resultDir):

	posDiNucDist = GetDistributionForFiles(resultDir, "Signal*.fa")
	negDiNucDist = GetDistributionForFiles(resultDir, "NoSignal*.fa")

	return posDiNucDist, negDiNucDist

def GenerateDiNucleotideGraphs(resultDir, currDir, posDiNucDist, negDiNucDist):

	targetDir = resultDir + "/" + "dinucleotide_distribution"
	if not os.path.exists(targetDir):
		os.makedirs(targetDir);

	# posPngFile = targetDir + "/" + currDir + "_pos.png"
	# negPngFile = targetDir + "/" + currDir + "_neg.png"
	graphFileName = targetDir + "/" + currDir + ".png"
	PlotGraphs(posDiNucDist, negDiNucDist, graphFileName)
	

def computeDiNucleotideDistribution(resultDir, level=1):
	posDistList = list()
	negDistList = list();

	Compute3UtrDistibution();

	num_sep = resultDir.count(os.path.sep)
	for root, dirs, files in os.walk(resultDir):
		num_sep_this = root.count(os.path.sep)
		if num_sep + level > num_sep_this:
			for dir in dirs:
				if isDirResultDir(dir) == False:
					posDiNucDist, negDiNucDist = GetDiNucleotideDistribution(root + "/" + dir)
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
