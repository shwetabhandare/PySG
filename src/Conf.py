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

	def GetAPercent(self):
		return self.confMap['sequence'][3]['A'][0]
	def GetTPercent(self):
		return self.confMap['sequence'][4]['T'][0]
	def GetGPercent(self):
		return self.confMap['sequence'][5]['G'][0]
	def GetCPercent(self):
		return self.confMap['sequence'][6]['C'][0]

	def GetMotifType(self):
		return self.confMap['motif'][0]['type'][0];

	def GetMotifLength(self):
		return self.confMap['motif'][1]['length'][0]

	def GetNumMotifs(self):
		return self.confMap['motif'][2]['numMotifs'][0]

	def GetDistanceBetweenMotifs(self):
		return self.confMap['motif'][3]['distance'][0]

	def GetMotifLocation(self):
		val =  self.confMap['motif'][4]['location'][0]
		if val == 'end':
			return self.confMap['motif'][4]['location'][0]['end'][0]
		elif val == 'start':
			return self.confMap['motif'][4]['location'][0]['start'][0]
		else:
			return val;

