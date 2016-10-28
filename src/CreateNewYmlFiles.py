import glob;
import os;
import re;
import subprocess;

def CreateNewYmlFiles(directory, copyPrefix, createPrefix, extension):
	os.chdir(directory)
	searchStr = copyPrefix + "*." + extension;
	print searchStr
	reStr = copyPrefix + "_" + "(.*)\." + extension;
	print reStr;

	for ymlFile in glob.glob(searchStr):
		m = re.search(reStr, ymlFile);
		newFile = createPrefix + "_" + m.group(1) + "." + extension;
		print "copying file : ", ymlFile, " to ", newFile;
		subprocess.call(["cp", ymlFile, newFile])

if __name__ == "__main__":
	import sys
	copyPrefix = sys.argv[1]
	createPrefix = sys.argv[2]
	extension = sys.argv[3]
	directory = sys.argv[4];

	#directory = "/projects/bhandare/workspace/PySG/src/resources";
	CreateNewYmlFiles(directory, copyPrefix, createPrefix, extension)

	
