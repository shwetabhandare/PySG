import yaml
import os
from os import path
import SeqGenUtils
import shutil;

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

class YamlFastaGenerator():
	alpha = list()
	numSeq = list()
	seqLen = list()
	pwmFiles = list()
	seqLen = []
	seqWithSignalPercent = []
	noSignalType = ""
	signalType = ""
	seqWithSignalPercent = list()
	pwmFileDirectory = ""
	utrDist = dict(
		A = 0.27,
		T = 0.22,
		G = 0.21,
		C = 0.30,
	)
	confMap = {}
	targetDir = ""

	def __init__(self, confFile):
		self.confMap = SeqGenUtils.GetConf(confFile)
		self.SetupTargetDir()
		#setupVariables(self.confMap)

	def GetAlpha(self):
		return self.alpha;
	def GetNumSeqList(self):
		return self.numSeq;
	def GetNumSeqLenList(self):
		return self.seqLen;
	def GetPwmFilesList(self):
		return self.pwmFiles;
	def GetSignalPercentList(self):
		return self.seqWithSignal;
	def GetNoSignalType(self):
		return self.noSignalType;
	def GetSignalType(self):
		return self.signalType;
	def GetUtrDistDict(self):
		return self.utrDist;
	def GetPwmDir(self):
		return pwmFileDirectory;
	def GetTargetDir(self):
		print "Target DIR: ", self.targetDir;
		return self.targetDir;

	def SetupTargetDir(self):
		if self.confMap.get("targetDir"):
			self.targetDir = self.confMap["targetDir"]
		else:
			self.targetDir = os.os.getcwd() + "/" + "tmp";

		if not os.path.exists(self.targetDir):
			print "Creating directory: ", targetDir
			os.makedirs(self.targetDir)
		else:
			self.deleteFiles()

		print "Target DIR: ", self.targetDir;

	def deleteFiles(self):
		for the_file in os.listdir(self.targetDir):
			file_path = os.path.join(self.targetDir, the_file)
			try:
				if os.path.isfile(file_path):
					os.unlink(file_path)
				elif os.path.isdir(file_path): 
					shutil.rmtree(file_path)
			except Exception as e:
				print(e)


	# def getNoSignalDictWithAlpha(location, numberSeq, seqLength, alphaValue, fileId):
	# 	outFileName = "NoSignal_" + fileId
	# 	data = dict(
	# 		  seqLen = seqLength,
	# 			 numSeq = numberSeq,
	# 		  outNoSignalFastaFile = location + "/" + outFileName + ".fa",
	# 		  seqBackGround = utrDist,
	# 		  alpha = alphaValue,
	# 	)
	# 	return data;

	# def getNoSignalDict(location, numberSeq, seqLength, inpFastaFile, fileId):
	# 	outFileName = "NoSignal_" + fileId
	# 	data = dict(
	# 		  seqLen = seqLength,
	# 		  numSeq = numberSeq,
	# 		  fastaFile = inpFastaFile,
	# 		  outNoSignalFastaFile = location + "/" + outFileName + ".fa",
	# 	)

	# 	return data;

	# def getMotifDict(location, numberSeq, seqLength, signalSeq, pwmFileName, fileId):
	# 	signalLocation = 10;
	# 	outFileName = "Signal_" + fileId;
	# 	pwmFileToAdd = pwmFileDirectory + "/" + pwmFileName;
	# 	data = dict(
	# 			outSignalFile = location + "/" + outFileName + ".fa",
	# 			seqBackGround = utrDist,
	# 			seqWithSignal = signalSeq,
	# 			locationFromStart = signalLocation,
	# 			pwmFile = pwmFileToAdd,
	# 		)
	# 	yamlFileName = fileId + ".yml"
	# 	return yamlFileName, data;

	# def writeYamlFile(location, yamlFileName, noSignalDict, motifDict):
	# 	from yaml.representer import Representer
	# 	import collections;
	# 	import sys
	# 	yaml.add_representer(collections.defaultdict, Representer.represent_dict)
	# 	data = dict(
	# 		sequence = dict(
	# 			nosignal = noSignalDict,
	# 			signal = motifDict,
	# 		)
	# 	)
	# 	noalias_dumper = yaml.dumper.SafeDumper
	# 	noalias_dumper.ignore_aliases = lambda self, data: True

	# 	with open(location + "/" + yamlFileName, 'w') as outfile:
	# 		 outfile.write(yaml.dump(data, default_flow_style=False, Dumper=noalias_dumper))

	# def getPwmFiles(pwmDir):
	# 	global pwmFiles;
	# 	for f in os.listdir(pwmDir):
	# 		if f.lower().endswith(('.pwm')):
	# 			pwmFiles.append(f)
	# 	print "PWM FILE: ", pwmFiles

	# 	return pwmFiles;

	# def GenerateNoSignalWithDirichlet(location):
	# 	pwmFiles = getPwmFiles(pwmFileDirectory);	
	# 	for numSeqIdx, i in enumerate(numSeq):
	# 		print str(numSeqIdx), str(i)
	# 		for numSeqLenIdx, j in enumerate(seqLen):
	# 			print str(numSeqLenIdx),  str(j)
	# 			for a in alpha:
	# 				print "ALPHA: ", str(a)
	# 				for signalPercent in seqWithSignalPercent:
	# 					print "signal percent: ", str(signalPercent)
	# 					for pwmFile in pwmFiles:
	# 						print "PWM: ", pwmFile
	# 						fileId = str(i) + "_" + str(j) + "_" + str(signalPercent) + "_" + str(a) + "_" + pwmFile;
	# 						print "FILE ID:", fileId;
	# 						yamlFileName, motifDict = getMotifDict(location, i, j, signalPercent, pwmFile, fileId);
	# 						print "YAML FILE:", yamlFileName;
	# 						noSignalDict = getNoSignalDictWithAlpha(location, i, j, a, fileId);
	# 						writeYamlFile(location, yamlFileName, noSignalDict, motifDict);

	# def GenerateNoSignalWithFasta(location):
	# 	inpFastaFile = "/projects/bhandare/workspace/scripts/NegFileCreator/3UTR_transcripts_Human.txt"
	# 	pwmFiles = getPwmFiles(pwmFileDirectory);	
	# 	for numSeqIdx, i in enumerate(numSeq):
	# 		for numSeqLenIdx, j in enumerate(seqLen):
	# 			for signalPercent in seqWithSignalPercent:
	# 				for pwmFile in pwmFiles:
	# 					fileId = str(i) + "_" + str(j) + "_" + str(signalPercent) + pwmFile;
	# 					yamlFileName, motifDict = getMotifDict(location, i, j, signalPercent, pwmFile, fileId);
	# 					noSignalDict = getNoSignalDict(location, i, j, inpFastaFile, fileId);
	# 					writeYamlFile(location, yamlFileName, noSignalDict, motifDict);

	# ## Program stats here.
	# def CreateConfFiles(targetDir, noSignalType):
	# 	print "NO SIGNAL TYPE:", noSignalType
	# 	if noSignalType == 'shuffle':
	# 		print "Shuffle type";
	# 		#GenerateNoSignalWithFasta(targetDir);
	# 	elif noSignalType == 'dirichlet':
	# 		print "dirichlet"
	# 		GenerateNoSignalWithDirichlet(targetDir);
	# 	else:
	# 		print "Invalid No Signal type", noSignalType;

	# def setupNoSignalVariables(confMap):
	# 	global noSignalType, alpha, noSignalFastaFile;
	# 	if confMap["nosignal"].get("type"):
	# 		noSignalType = confMap["nosignal"]["type"]
	# 		if noSignalType == "dirichlet":
	# 			alpha = confMap["nosignal"]["alpha"]
	# 		elif noSignalType == "shuffle":
	# 			noSignalFastaFile = confMap["nosignal"]["fastaFile"]
	# 		else:
	# 			print "Did not specify nosignal type"
	# 			return;
	# 	return noSignalType;

	# def getSignalType(confMap):
	# 	global signalType;

	# 	if confMap["signal"].get("type"):
	# 		signalType = confMap["signal"]["type"]
	# 		if signalType == 'pwmFiles':
	# 			pwmFiles = confMap["signal"]["pwmFiles"]
	# 		elif signalType == "pwmDir":
	# 			pwmDir =  confMap["signal"]["pwmDir"]
	# 			pwmFiles = getPwmFiles(pwmDir);
	# 		elif signalType == "textMotif":
	# 			textMotifs = confMap["signal"]["textMotif"]
	# 		elif signalType == "kmers":
	# 			kmers = confMap["signal"]["kmers"]
	# 		else:
	# 			print "Invalid signal Type:", signalType;
	# 	return signalType;

	# def setupSignalVariables(confMap):
	# 	global signalType, signalLocation, seqWithSignalPercent, pwmFiles;

	# 	signalType = getSignalType(confMap)

	# 	if confMap["signal"].get("location"):
	# 		signalLocation = confMap["signal"]["location"]
	# 	else:
	# 		signalLocation = [0]
		
	# 	if confMap["signal"].get("signalPercent"):
	# 		seqWithSignalPercent = confMap["signal"]["signalPercent"]		
	# 	else:
	# 		seqWithSignalPercent = [100]


	# 	return signalType;


	# def setupVariables(confMap):
	# 	numSeq = confMap["numSeq"]
	# 	seqLen = confMap["seqLen"]
	# 	signalType = setupSignalVariables(confMap)
	# 	noSignalType = setupNoSignalVariables(confMap)
	# 	return signalType, noSignalType;





	# def main():
	# 	import sys
	# 	confFile = sys.argv[1]
	# 	confMap = SeqGenUtils.GetConf(confFile)
	# 	targetDir = getTargetDir(confMap)
	# 	signalType, noSignalType = setupVariables(confMap)
	# 	CreateConfFiles(targetDir, noSignalType)

if __name__ == "__main__":
	main();