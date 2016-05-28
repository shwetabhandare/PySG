import re;
import numpy as np;
import os;

def read_file(dreme_file):
	pfmList = []
	for line in open(dreme_file): # no need to use readlines if you don't want to store them
		# use a list comprehension to build your array on the fly
		pwmLine = [];
		tokens = line.split();
		for token in tokens:
			token = token.strip();
			pwmLine.append(float(token));

		pfmList.append(pwmLine);

	return np.transpose(pfmList).tolist();

def writePwmFile(pfmList, fileName):
	count = 0;
	target = open(fileName, 'w')
	firstLine = "Gene:\t" + os.path.basename(fileName);
	target.write(firstLine)
	target.write("\n");

	print type(pfmList)
	for pfmLine in pfmList:
		print pfmLine;
		pfmLineAsStrings = ["%.4f" % number for number in pfmLine]

		if count == 0:

			lineToAppend = "A:" + "\t".join(pfmLineAsStrings);
			count = count + 1;
		elif count == 1:
			lineToAppend = "C:" + "\t".join(pfmLineAsStrings);
			count = count + 1;
		elif count == 2:
			lineToAppend = "G:" + "\t".join(pfmLineAsStrings);
			count = count + 1;
		elif count == 3:
			lineToAppend = "T:" + "\t".join(pfmLineAsStrings);
			count = count + 1;
		target.write(lineToAppend)
		target.write("\n")

	target.close()

def createPFMFile(pfmFile, outFile):
	pfmLines = read_file(pfmFile)
	print type(pfmLines)
	writePwmFile(pfmLines, outFile)

if __name__ == "__main__":
	import sys
	pfmFile = sys.argv[1]
	outFile = sys.argv[2]
	createPFMFile(pfmFile, outFile)
