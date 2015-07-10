import Conf;
from Conf import *;

conf = Conf("../conf/SeqGen.yaml");
class SeqGen():
	def __init__(self, conf):
		self.numSeqs = conf.GetNumSeq();
		self.minLen = conf.GetMinSeqLen();
		self.maxLen = conf.GetMaxSeqLen();



	def GenerateRandomSequences(self):
		print "Generating ", self.numSeqs, self.numSeqs, self.minLen, self.maxLen;
            #for num in range(0, self.numSeqs):
#randomLength = randint(self.minLen, self.maxLen);

