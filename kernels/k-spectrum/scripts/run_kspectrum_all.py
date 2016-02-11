import os;
import sys;
import subprocess;

seqFile = sys.argv[1];
isRNA = sys.argv[2];

kmerList = [[5,7],[9,16]];
fileFields = seqFile.split("_");
posLen = fileFields[-3];
negLen = fileFields[-2];

#print posLen, negLen;
for kmerPair in kmerList:
	for count in range(0,2):
		print kmerPair[0], kmerPair[1], posLen, negLen, 10, count, isRNA
		k1=int(kmerPair[0]);
		k2=int(kmerPair[1]);
		#python demo_spectrum_ver2.py seqFile kmerPair[0] kmerPair[1] posLen negLen 10 count isRNA
		#subprocess.call("python demo_spectrum_ver2.py seqFile k1 k2 posLen negLen 10 count isRNA", shell=True);

