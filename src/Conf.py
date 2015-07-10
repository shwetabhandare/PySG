import yaml


class Conf():

	confMap = {}
	def __init__(self, configFile):
		self.configFile = configFile;
		f = open(configFile)
		self.confMap = yaml.load(f)
		f.close()

	def GetMinSeqLen(self):
		return self.confMap['sequence'][0]['minLen'][0]

	def GetMaxSeqLen(self):
		return self.confMap['sequence'][1]['maxLen'][0]

	def GetNumSeq(self):
		return self.confMap['sequence'][2]['numSeq'][0]










