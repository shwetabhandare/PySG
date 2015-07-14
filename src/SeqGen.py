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
		self.aPercent = conf.GetAPercent();
		self.tPercent = conf.GetTPercent();
		self.gPercent = conf.GetGPercent();
		self.cPercent = conf.GetCPercent();


	def GetRandomSequence(self, seqLen):
		seq=""
		for count in range(seqLen):
			seq+=SeqGenUtils.weightedchoice([("C", self.cPercent), ("G", self.gPercent), ("A",
				self.aPercent), ("T", self.tPercent)]);
		return seq;

	def GenerateRandomSequences(self):
		print "Generating ", self.numSeqs, self.numSeqs, self.minLen, self.maxLen;
		for num in range(0, self.numSeqs): 
			randomLength = numpy.random.randint(self.minLen, self.maxLen); 
			randomSeq = self.GetRandomSequence(randomLength); 
			print "Random Sequence " , num, ": ", randomSeq;

