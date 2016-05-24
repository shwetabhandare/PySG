import sys;
import os, glob;
import parseResults;
import shutil;
import datetime;


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
	print "calling walklevel"
	parseSubDirectories(resultDir)
