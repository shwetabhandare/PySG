import sys
import os
import generateGraphs

def writePredictedMotifsToFile(dremeResultFile, dremePredictedFile, firstTime, range1, range2):
	if firstTime == True:
		fo = open(dremeResultFile, "w")
		firstTime = False;
	else:
		fo = open(dremeResultFile, "a")

	N = 100;
	from itertools import islice	
	with open(dremePredictedFile, "r") as text_file:
		head = list(islice(text_file, range1, range2));
		fo.write(''.join(head))

	fo.write("----------------------------------------------------------------")
	# Close opend file
	fo.close()
	return firstTime;

def parseResultDirectory(resultDir):

	firstTimeDreme = True;
	firstTimeKspectrum = True;

	#print resultDir;

	dremeResultFile = resultDir + "/predictedDremeMotifs.txt"
	kspectrumResultFile = resultDir + "/predictedKspectrumMotifs.txt"

	for subdir, dirs, files in os.walk(resultDir):
	
		for file in files:
			filepath = subdir + os.sep + file
			if filepath.endswith("dreme.txt"):
				firstTimeDreme = writePredictedMotifsToFile(dremeResultFile, filepath, firstTimeDreme, 10, 100)

	 		if filepath.endswith(".dat"):
				firstTimeKspectrum = writePredictedMotifsToFile(kspectrumResultFile, filepath, firstTimeKspectrum, 0, 100);




def printPredictedMotifs(resultDir, level=1):
	num_sep = resultDir.count(os.path.sep)
	for root, dirs, files in os.walk(resultDir):
		num_sep_this = root.count(os.path.sep)
		if num_sep + level > num_sep_this:
			for dir in dirs:
				if generateGraphs.isDirResultDir(dir):
					numSubDirs = generateGraphs.GetSubDirCount(root + "/" + dir)
					print "Found Result Directory: ", dir, ", num sub-dir: ", str(numSubDirs)
					if numSubDirs == 5:
						parseResultDirectory(root + "/" + dir)
					else:
						print "Look at directory: ", dir;
						continue;

if __name__ == "__main__":
	resultDir = sys.argv[1]
	printPredictedMotifs(resultDir)
	#computeDiNucleotideDistribution(resultDir)