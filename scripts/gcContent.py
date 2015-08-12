import re
import string
import sys
 
def main():
	file = open(sys.argv[1], "r")
	gcCount = 0
	totalBaseCount = 0
	for line in file:
		line = line.strip("\n")
		if not line.startswith(">"):
			gcCount += len(re.findall("[GC]", line))
			totalBaseCount += len(re.findall("[GCTA]", line))
			gcFraction = float(gcCount) / totalBaseCount
			print( gcFraction * 100 )


if __name__ == '__main__':
	main()
