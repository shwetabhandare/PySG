import Conf;
from Conf import *;
import SeqGenUtils;
from SeqGenUtils import *;
import random;
import numpy;
from random import *;

conf = Conf("../conf/SeqGen.yaml");
class SeqGen():
	def __init__(self, conf):
		self.numSeqs = conf.GetNumSeq();
		self.minLen = conf.GetMinSeqLen();
		self.maxLen = conf.GetMaxSeqLen();



	def GenerateRandomSequences(self):
		print "Generating ", self.numSeqs, self.numSeqs, self.minLen, self.maxLen;
		for num in range(0, self.numSeqs): 
			randomLength = numpy.random.randint(self.minLen, self.maxLen); 
			randomSeq = SeqGenUtils.GetRandomSequence(randomLength); 
			print "Random Sequence: " , num, randomSeq;

