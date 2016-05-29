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

def PlotGraphs(diNucDist, fileName):
	global threeUtrValues, threeUtrErrors;

	import matplotlib
	matplotlib.use('Agg')
	
	from matplotlib import pyplot as plt	

	meanValues = [x[0] for x in diNucDist.values()]
	errorValues = [x[1] for x in diNucDist.values()]	



	fig = plt.figure()
	plt.title("Di-Nucleotide Distributions");
	plt.ylabel("Distribution");
	plt.xlabel("Di-nucleotides")
	labels = diNucDist.keys();
	xValues = np.arange(len(meanValues))

	#print labels;
	#print meanValues
	#print errorValues;


	plt.xticks(range(len(labels)), list(labels));	

	eb1 = plt.errorbar(xValues, meanValues, yerr=errorValues, xerr=None, fmt='', color='b')
	eb2 = plt.errorbar(xValues, threeUtrValues, yerr=threeUtrErrors, xerr=None, fmt='', color='g')

	plt.legend([eb1, eb2], ['Generated Files', "3'UTR"])

	plt.savefig(fileName);
	#plt.show();
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

	posPngFile = targetDir + "/" + currDir + "_pos.png"
	negPngFile = targetDir + "/" + currDir + "_neg.png"

	PlotGraphs(posDiNucDist, posPngFile)
	PlotGraphs(negDiNucDist, negPngFile)
	

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
	computeDiNucleotideDistribution(resultDir)
