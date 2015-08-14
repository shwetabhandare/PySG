import yaml


class Conf():

	confMap = {}
	def __init__(self, configFile):
		self.configFile = configFile;
		f = open(configFile)
		self.confMap = yaml.load(f)
		f.close()

	#{'sequence': {'A': 30, 'maxLen': 350, 'C': 20, 'numSeq': 10, 'T': 30, 'G': 20, 'minLen': 250}}
	def GetMinSeqLen(self):
		return self.confMap['sequence']['minLen']

	def GetMaxSeqLen(self):
		return self.confMap['sequence']['maxLen']

	def GetNumSeq(self):
		return self.confMap['sequence']['numSeq']

	def GetAPercent(self):
		return self.confMap['sequence']['A']

	def GetTPercent(self):
		return self.confMap['sequence']['T']
	def GetGPercent(self):
		return self.confMap['sequence']['G']
	def GetCPercent(self):
		return self.confMap['sequence']['C']

	def GetMotifType(self):
		if 'motif' in self.confMap.keys():
			return self.confMap['motif']['type'];
		else:
			return "";

	def GetMotifLength(self):
		if 'motif' in self.confMap.keys():
			return self.confMap['motif']['length']
		else:
			return 0;

	def GetNumMotifs(self):
		if 'motif' in self.confMap.keys():
			return self.confMap['motif']['numMotifs']
		else:
			return 0;

	def GetDistanceBetweenMotifs(self):
		if 'motif' in self.confMap.keys():
			return self.confMap['motif']['distance']
		else:
			return 0

	def GetMotifLocation(self):
		if 'motif' in self.confMap.keys():
			val =  self.confMap['motif']['location']
			if val == 'end':
				return self.confMap['motif'][4]['location'][0]['end'][0]
			elif val == 'start':
				return self.confMap['motif'][4]['location'][0]['start'][0]
			else:
				return val;
		else:
			return 0;
	def GetMotif(self):
		if 'motif' in self.confMap.keys():
			return self.confMap['motif']['motif'];
		else:
			return ""
