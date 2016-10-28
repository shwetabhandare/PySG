import re
import string
import sys
 

def getNucleotideComposition(file):
	seqFile = open(file, "r")
	gcCount = 0
	atCount = 0
	aCount = gCount = cCount = uCount = 0;
	totalBaseCount = 0
	for line in seqFile:
		line = line.strip("\n")
		if not line.startswith(">"):
			gcCount += len(re.findall("[GC]", line))
			atCount += len(re.findall("[AT]", line))
			aCount += len(re.findall("[A]", line))
			uCount += len(re.findall("[T]", line))
			cCount += len(re.findall("[C]", line))
			gCount += len(re.findall("[G]", line))
			totalBaseCount += len(re.findall("[GCTA]", line))

	#print "Total Base count: ", totalBaseCount;
	gcFraction = float(gcCount) / totalBaseCount
	atFraction = float(atCount) / totalBaseCount

	aFraction = float(aCount) / totalBaseCount
	uFraction = float(uCount) / totalBaseCount
	gFraction = float(gCount) / totalBaseCount
	cFraction = float(cCount) / totalBaseCount

	print "AT Content: " + str(atFraction * 100);
	print "GC Content: " + str(gcFraction * 100);

	#print "A Content: " + str(aFraction * 100);
	#print "U Content: " + str(uFraction * 100);
	#print "G Content: " + str(gFraction * 100);
	#print "C Content: " + str(cFraction * 100);

	return gcFraction * 100, atFraction * 100;
	
def main():
	gcContent, atContent = getNucleotideComposition(file);
	#print gcContent, atContent;


if __name__ == '__main__':
	main()
