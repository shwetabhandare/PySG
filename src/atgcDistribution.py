import re
import string
import sys
 
def printATGCDistribution(fastaFile):
	file = open(fastaFile, "r")
	gcCount = 0
	atCount = 0
	totalBaseCount = 0
	for line in file:
		line = line.strip("\n")
		if not line.startswith(">"):
			gcCount += len(re.findall("[GC]", line))
			atCount += len(re.findall("[AT]", line))
			totalBaseCount += len(re.findall("[GCTA]", line))

	gcFraction = float(gcCount) / totalBaseCount
	atFraction = float(atCount) / totalBaseCount

	print "AT Content: " + str(atFraction * 100);
	print "GC Content: " + str(gcFraction * 100);


if __name__ == '__main__':
	import sys
	printATGCDistribution(sys.argv[1])
