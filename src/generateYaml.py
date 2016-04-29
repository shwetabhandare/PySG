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
	noSignalFastaFile = ""

	def __init__(self, confFile):
		self.confMap = SeqGenUtils.GetConf(confFile)
		self.SetupTargetDir()
		self.setupVariables()

	def GetAlpha(self):
		return self.alpha;
	def GetNumSeqList(self):
		return self.numSeq;
	def GetNumSeqLenList(self):
		return self.seqLen;
	def GetPwmFilesList(self):
		return self.pwmFiles;
	def GetSignalPercentList(self):
		return self.seqWithSignalPercent;
	def GetNoSignalType(self):
		return self.noSignalType;
	def GetSignalType(self):
		return self.signalType;
	def GetUtrDistDict(self):
		return self.utrDist;
	def GetPwmDir(self):
		return self.pwmFileDirectory;
	def GetTargetDir(self):
		return self.targetDir;
	def GetNoSignalFastaFile(self):
		return self.noSignalFastaFile;

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

	def getSignalType(self):

		if self.confMap["signal"].get("type"):
			self.signalType = self.confMap["signal"]["type"]

			if self.signalType == 'pwmFiles':
				self.pwmFiles = self.confMap["signal"]["pwmFiles"]
			elif self.signalType == "pwmDir":
				self.pwmFileDirectory =  self.confMap["signal"]["pwmDir"]
				self.pwmFiles = self.getPwmFiles();
			elif self.signalType == "textMotif":
				self.textMotifs = self.confMap["signal"]["textMotif"]
			elif self.signalType == "kmers":
				self.kmers = confMap["signal"]["kmers"]
			else:
				print "Invalid signal Type:", self.signalType;

	def setupSignalVariables(self):
		self.getSignalType();

		if self.confMap["signal"].get("location"):
			self.signalLocation = self.confMap["signal"]["location"]
		else:
			self.signalLocation = [0]
		
		if self.confMap["signal"].get("signalPercent"):
			self.seqWithSignalPercent = self.confMap["signal"]["signalPercent"]		
		else:
			self.seqWithSignalPercent = [100]

	def getPwmFiles(self):
		for f in os.listdir(self.pwmFileDirectory):
			if f.lower().endswith(('.pwm')):
				self.pwmFiles.append(f)

	def setupVariables(self):
		self.numSeq = self.confMap["numSeq"]
		self.seqLen = self.confMap["seqLen"]
		self.setupSignalVariables()
		self.setupNoSignalVariables()

	def setupNoSignalVariables(self):
		if self.confMap["nosignal"].get("type"):
			self.noSignalType = self.confMap["nosignal"]["type"]
			if self.noSignalType == "dirichlet":
				self.alpha = self.confMap["nosignal"]["alpha"]
			elif self.noSignalType == "shuffle":
				self.noSignalFastaFile = self.confMap["nosignal"]["fastaFile"]
			else:
				print "Did not specify nosignal type"
				return;

	## Program stats here.
	def CreateConfFiles(self):
		if self.noSignalType == 'shuffle':
			print "Shuffle type";
			#GenerateNoSignalWithFasta(targetDir);
		elif self.noSignalType == 'dirichlet':
			self.GenerateNoSignalWithDirichlet();
		else:
			print "Invalid No Signal type", self.noSignalType;


	def GenerateNoSignalWithDirichlet(self):
		location = self.targetDir;
		pwmFiles = self.pwmFiles;
		numSeq = self.GetNumSeqList();
		seqLen = self.GetNumSeqLenList();
		alpha = self.GetAlpha();
		seqWithSignalPercent = self.GetSignalPercentList();

		for numSeqIdx, i in enumerate(numSeq):
			print str(numSeqIdx), str(i)
			for numSeqLenIdx, j in enumerate(seqLen):
				print str(numSeqLenIdx),  str(j)
				for a in alpha:
					#print "ALPHA: ", str(a)
					for signalPercent in seqWithSignalPercent:
						#print "signal percent: ", str(signalPercent)
						for pwmFile in pwmFiles:
							#print "PWM: ", pwmFile
							if os.path.isabs(pwmFile) == True:
								self.pwmFileDirectory = os.path.dirname(pwmFile);
								pwmFile = os.path.basename(pwmFile);

							fileId = str(i) + "_" + str(j) + "_" + str(a) + "_" + str(signalPercent) + "_" + pwmFile;
							#print "FILE ID:", fileId;
							yamlFileName, motifDict = self.getMotifDict(location, i, j, signalPercent, pwmFile, fileId);
							#print "YAML FILE:", yamlFileName;
							noSignalDict = self.getNoSignalDictWithAlpha(location, i, j, a, fileId);
							self.writeYamlFile(location, yamlFileName, noSignalDict, motifDict);

	def getMotifDict(self, location, numberSeq, seqLength, signalSeq, pwmFileName, fileId):
		signalLocation = 10;
		outFileName = "Signal_" + fileId;
		#print "PWM FILENAME: ", pwmFileName;

		pwmFileToAdd = self.pwmFileDirectory + "/" + pwmFileName;	
		data = dict(
				outSignalFile = location + "/" + outFileName + ".fa",
				seqBackGround = self.utrDist,
				seqWithSignal = signalSeq,
				locationFromStart = signalLocation,
				pwmFile = pwmFileToAdd,
			)
		yamlFileName = fileId + ".yml"
		return yamlFileName, data;

	def getNoSignalDictWithAlpha(self, location, numberSeq, seqLength, alphaValue, fileId):
		outFileName = "NoSignal_" + fileId
		data = dict(
			seqLen = seqLength,
			numSeq = numberSeq,
			outNoSignalFastaFile = location + "/" + outFileName + ".fa",
			seqBackGround = self.utrDist,
			alpha = alphaValue,
		)
		return data;

	def writeYamlFile(self, location, yamlFileName, noSignalDict, motifDict):
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


	# def getNoSignalDict(location, numberSeq, seqLength, inpFastaFile, fileId):
	# 	outFileName = "NoSignal_" + fileId
	# 	data = dict(
	# 		  seqLen = seqLength,
	# 		  numSeq = numberSeq,
	# 		  fastaFile = inpFastaFile,
	# 		  outNoSignalFastaFile = location + "/" + outFileName + ".fa",
	# 	)

	# 	return data;








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










	# def main():
	# 	import sys
	# 	confFile = sys.argv[1]
	# 	confMap = SeqGenUtils.GetConf(confFile)
	# 	targetDir = getTargetDir(confMap)
	# 	signalType, noSignalType = setupVariables(confMap)
	# 	CreateConfFiles(targetDir, noSignalType)

if __name__ == "__main__":
	main();