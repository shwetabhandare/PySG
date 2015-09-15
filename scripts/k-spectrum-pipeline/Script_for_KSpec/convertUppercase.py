import os
import sys


inputFile = open(sys.argv[1], "r")
content = inputFile.read()
with open(sys.argv[1], "wb") as outputFile:
	outputFile.write(content.upper())
