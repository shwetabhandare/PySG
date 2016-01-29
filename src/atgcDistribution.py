import re
import string
import sys
 
def printATGCDistribution(fastaFile):
	file = open(fastaFile, "r")
	gcCount = 0
	atCount = 0
	aCount = tCount = gCount = cCount = 0
	totalBaseCount = 0
	for line in file:
		line = line.strip("\n")
		if not line.startswith(">"):
			line = line.upper();
			gcCount += len(re.findall("[GC]", line))
			atCount += len(re.findall("[AT]", line))
			aCount += len(re.findall("[A]", line))
			tCount += len(re.findall("[T]", line))
			gCount += len(re.findall("[G]", line))
			cCount += len(re.findall("[C]", line))
			totalBaseCount += len(re.findall("[GCTA]", line))

	gcFraction = float(gcCount) / totalBaseCount
	atFraction = float(atCount) / totalBaseCount

	aFraction = float(aCount) / totalBaseCount
	tFraction = float(tCount) / totalBaseCount
	gFraction = float(gCount) / totalBaseCount
	cFraction = float(cCount) / totalBaseCount

	print "AT Content: " + str(atFraction * 100);
	print "GC Content: " + str(gcFraction * 100);

	print "A Content: " + str(aFraction * 100);
	print "T Content: " + str(tFraction * 100);
	print "G Content: " + str(gFraction * 100);
	print "C Content: " + str(cFraction * 100);

	total = aFraction + tFraction + gFraction + cFraction;
	print "Total : " + str(total);


if __name__ == '__main__':
	import sys
	printATGCDistribution(sys.argv[1])
