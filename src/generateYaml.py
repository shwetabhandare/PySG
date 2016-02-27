import yaml
import os
from os import path

alpha = [0.1, 1, 10, 100, 1000]
#numSeq = [1000, 2000, 3000, 4000, 5000]
#seqLen = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
numSeq = [100, 200, 300, 400, 500]
seqLen = [200, 300, 400, 500, 600]
seqWithSignalPercent = [0, 50, 75, 90, 100]
pwmFileDirectory = "/projects/bhandare/workspace/PySG/data/pwm"
utrDist = dict(
	seqBackGround = dict(
		A = 0.27,
		T = 0.22,
		G = 0.21,
		C = 0.30,
	)
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


def getNoSignalDictWithAlpha(location, numberSeq, seqLength, alphaValue):
	outFileName = "NoSignal_" + str(numberSeq) + "_" + str(seqLength) + "_" + str(alpha)
	data = dict(
		  seqLen = seqLength,
		  numSeq = numberSeq,
		  outNoSignalFastaFile = location + "/" + outFileName + ".fa",
		  seqBackGround = utrDist,
		  alpha = alphaValue,
	)
	return data;

def getNoSignalDict(location, numberSeq, seqLength, inpFastaFile):
	outFileName = "NoSignal_" + str(numberSeq) + "_" + str(seqLength)
	data = dict(
		  seqLen = seqLength,
		  numSeq = numberSeq,
		  fastaFile = inpFastaFile,
		  outNoSignalFastaFile = location + "/" + outFileName + ".fa",
	)

	return data;

def getMotifDict(location, numberSeq, seqLength, signalSeq, pwmFileName):
	signalLocation = 10;
	outFileName = "Motif_" + str(numberSeq) + "_" + str(seqLength) + "_" + str(signalSeq) + "_" + str(signalLocation) + "_" + pwmFileName;
	pwmFileToAdd = pwmFileDirectory + "/" + pwmFileName;
	data = dict(
			outSignalFile = location + "/" + outFileName + ".fa",
			seqBackGround = utrDist,
			seqWithSignal = signalSeq,
			locationFromStart = signalLocation,
			pwmFile = pwmFileToAdd,
		)
	yamlFileName = str(numberSeq) + "_" + str(seqLength) + "_" + str(signalSeq) + "_" + str(signalLocation) + "_" + pwmFileName + ".yml"
	return yamlFileName, data;

def writeYamlFile(location, yamlFileName, noSignalDict, motifDict):
	data = dict(
		sequence = dict(
				nosignal = noSignalDict,
				signal = motifDict,
		)
	)

	with open(location + "/" + yamlFileName, 'w') as outfile:
		 outfile.write( yaml.dump(data, default_flow_style=True) )

def getPwmFiles(pwmDir):
	pwmFiles = [];
	for f in os.listdir(pwmDir):
		if f.lower().endswith(('.pwm')):
			pwmFiles.append(f)
	print pwmFiles;
	return pwmFiles;

def GenerateNoSignalWithDirichlet(location):
	pwmFiles = getPwmFiles(pwmFileDirectory);	
	for numSeqIdx, i in enumerate(numSeq):
		for numSeqLenIdx, j in enumerate(seqLen):
			for a in alpha:
				noSignalDict = getNoSignalDictWithAlpha(location, i, j, a);
				for signalPercent in seqWithSignalPercent:
					print signalPercent
					for pwmFile in pwmFiles:
						yamlFileName, motifDict = getMotifDict(location, i, j, signalPercent, pwmFile);
						writeYamlFile(location, yamlFileName, noSignalDict, motifDict);

def GenerateNoSignalWithFasta(location):
	inpFastaFile = "/projects/bhandare/workspace/scripts/NegFileCreator/3UTR_transcripts_Human.txt"
	pwmFiles = getPwmFiles(pwmFileDirectory);	
	for numSeqIdx, i in enumerate(numSeq):
		for numSeqLenIdx, j in enumerate(seqLen):
			noSignalDict = getNoSignalDict(location, i, j, inpFastaFile);
			for signalPercent in seqWithSignalPercent:
				print signalPercent
				for pwmFile in pwmFiles:
					yamlFileName, motifDict = getMotifDict(location, i, j, signalPercent, pwmFile);
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
