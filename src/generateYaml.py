import yaml
import os
from os import path

#alpha = [0.1, 1, 10, 100, 1000]
#numSeq = [1000, 2000, 3000, 4000, 5000]
#seqLen = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
alpha = [100, 1000]
numSeq = [100, 200]
seqLen = [500, 600]
seqWithSignalPercent = [90]
pwmFileDirectory = "/projects/bhandare/workspace/PySG/data/pwm"
utrDist = dict(
	A = 0.27,
	T = 0.22,
	G = 0.21,
	C = 0.30,
)

# generate A, T, G, C using s = np.random.dirichlet((0.1,0.1,0.1,0.1),1)
# expect sparse distribution, there is no apriori information.
# tp = np.random.dirichlet((0.1,0.1,0.1,0.1),1)[0]
# tp
#array([  4.61012266e-01,   5.37356121e-01,   1.59085119e-03,
#         4.07617397e-05])
# hur_tp = np.random.dirichlet(tp,1)[0]
# hur_tp
#array([ 0.06761885,  0.93238115,  0.        ,  0.        ])
# hur_tp = np.random.dirichlet(tp*100,1)[0]
# hur_tp
#array([  4.13762328e-01,   5.85305180e-01,   9.32492390e-04,
#         1.97868086e-59])
# ttp_tp = np.random.dirichlet(tp*100,1)[0]
# ttp_tp
#array([  5.14829673e-01,   4.85166046e-01,   4.28061017e-06,
#         1.59080129e-58])
# Feed hur_tp, ttp_tp to the sequence generator for A, T, G, C percents.
# when gamma is large, the ATGC percentages will be similar, and when its smaller they will differ more.
# here gammma = 100.



def getNoSignalDictWithAlpha(location, numberSeq, seqLength, alphaValue, fileId):
	outFileName = "NoSignal_" + fileId
	data = dict(
		  seqLen = seqLength,
			 numSeq = numberSeq,
		  outNoSignalFastaFile = location + "/" + outFileName + ".fa",
		  seqBackGround = utrDist,
		  alpha = alphaValue,
	)
	return data;

def getNoSignalDict(location, numberSeq, seqLength, inpFastaFile, fileId):
	outFileName = "NoSignal_" + fileId
	data = dict(
		  seqLen = seqLength,
		  numSeq = numberSeq,
		  fastaFile = inpFastaFile,
		  outNoSignalFastaFile = location + "/" + outFileName + ".fa",
	)

	return data;

def getMotifDict(location, numberSeq, seqLength, signalSeq, pwmFileName, fileId):
	signalLocation = 10;
	outFileName = "Signal_" + fileId;
	pwmFileToAdd = pwmFileDirectory + "/" + pwmFileName;
	data = dict(
			outSignalFile = location + "/" + outFileName + ".fa",
			seqBackGround = utrDist,
			seqWithSignal = signalSeq,
			locationFromStart = signalLocation,
			pwmFile = pwmFileToAdd,
		)
	yamlFileName = fileId + ".yml"
	return yamlFileName, data;

def writeYamlFile(location, yamlFileName, noSignalDict, motifDict):
	from yaml.representer import Representer
	import collections;
	import sys
	yaml.add_representer(collections.defaultdict, Representer.represent_dict)
	data = dict(
		sequence = dict(
			nosignal = noSignalDict,
			signal = motifDict,
		)
	)
	noalias_dumper = yaml.dumper.SafeDumper
	noalias_dumper.ignore_aliases = lambda self, data: True

	with open(location + "/" + yamlFileName, 'w') as outfile:
		 outfile.write(yaml.dump(data, default_flow_style=False, Dumper=noalias_dumper))

def getPwmFiles(pwmDir):
	pwmFiles = [];
	for f in os.listdir(pwmDir):
		if f.lower().endswith(('.pwm')):
			pwmFiles.append(f)
	return pwmFiles;

def GenerateNoSignalWithDirichlet(location):
	pwmFiles = getPwmFiles(pwmFileDirectory);	
	for numSeqIdx, i in enumerate(numSeq):
		print str(numSeqIdx), str(i)
		for numSeqLenIdx, j in enumerate(seqLen):
			print str(numSeqLenIdx),  str(j)
			for a in alpha:
				print "ALPHA: ", str(a)
				for signalPercent in seqWithSignalPercent:
					print "signal percent: ", str(signalPercent)
					for pwmFile in pwmFiles:
						print "PWM: ", pwmFile
						fileId = str(i) + "_" + str(j) + "_" + str(signalPercent) + "_" + str(a) + "_" + pwmFile;
						print "FILE ID:", fileId;
						yamlFileName, motifDict = getMotifDict(location, i, j, signalPercent, pwmFile, fileId);
						print "YAML FILE:", yamlFileName;
						noSignalDict = getNoSignalDictWithAlpha(location, i, j, a, fileId);
						writeYamlFile(location, yamlFileName, noSignalDict, motifDict);

def GenerateNoSignalWithFasta(location):
	inpFastaFile = "/projects/bhandare/workspace/scripts/NegFileCreator/3UTR_transcripts_Human.txt"
	pwmFiles = getPwmFiles(pwmFileDirectory);	
	for numSeqIdx, i in enumerate(numSeq):
		for numSeqLenIdx, j in enumerate(seqLen):
			for signalPercent in seqWithSignalPercent:
				for pwmFile in pwmFiles:
					fileId = str(i) + "_" + str(j) + "_" + str(signalPercent) + pwmFile;
					yamlFileName, motifDict = getMotifDict(location, i, j, signalPercent, pwmFile, fileId);
					noSignalDict = getNoSignalDict(location, i, j, inpFastaFile, fileId);
					writeYamlFile(location, yamlFileName, noSignalDict, motifDict);

## Program stats here.
def CreateConfFiles(location, type):
	if type == 'shuffle':
		GenerateNoSignalWithFasta(location);
	elif type == 'dirichlet':
		GenerateNoSignalWithDirichlet(location);

if __name__ == "__main__":
	import sys
	location = sys.argv[1]
	type = sys.argv[2]
	CreateConfFiles(location, type)
