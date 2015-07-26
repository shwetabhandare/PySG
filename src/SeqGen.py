import Conf;
from Conf import *;
import SeqGenUtils;
from SeqGenUtils import *;
import random;
import numpy;
from random import *;

#conf = Conf("../conf/SeqGen.yaml");
#conf = Conf("/tmp/seq.yml")
class SeqGen():
	def __init__(self, conf):
		self.numSeqs = conf.GetNumSeq();
		self.minLen = conf.GetMinSeqLen();
		self.maxLen = conf.GetMaxSeqLen();
		self.aPercent = conf.GetAPercent();
		self.tPercent = conf.GetTPercent();
		self.gPercent = conf.GetGPercent();
		self.cPercent = conf.GetCPercent();

		self.motifType = conf.GetMotifType();
		self.motifLength = conf.GetMotifLength();
		self.numMotifs = conf.GetNumMotifs();
		self.distanceBetweenMotifs = conf.GetDistanceBetweenMotifs();
		self.motifLocation = conf.GetMotifLocation();

		self.negativeSet = []
		self.positiveSet = []

		self.posFileName = "";
		self.negFileName = "";
		self.motif = conf.GetMotif();

		print "Generating ", self.numSeqs, " sequences, Motif Type: ", self.motifType, ", Motif Length: ", self.motifLength, ", Number of Motifs: ", self.numMotifs;

	def SetPosFileName(self, posFileName):
		self.posFileName = posFileName;

	def SetNegFileName(self, negFileName):
		self.negFileName = negFileName;

	def writePositiveFile(self):
		target = open(self.posFileName, 'w')
		for idx, seq in enumerate(self.positiveSet):
			header = ">PySeq_Pos_" + str(idx);
			target.write(header)
			target.write("\n");
			target.write(seq);
			target.write("\n")

	def writeNegativeFile(self):
		target = open(self.negFileName, 'w')
		for idx, seq in enumerate(self.negativeSet):
			header = ">PySeq_Neg_" + str(idx);
			target.write(header)
			target.write("\n");
			target.write(seq);
			target.write("\n")

	def GetRandomSequence(self, seqLen):
		seq=""
		for count in range(seqLen):
			seq+=SeqGenUtils.weightedchoice([("C", self.cPercent), ("G", self.gPercent), ("A",
				self.aPercent), ("T", self.tPercent)]);
		return seq;

	def GetPositiveSet(self):
		return self.positiveSet;

	def GetNegativeSet(self):
		return self.negativeSet;

	def GenerateRandomSequences(self, setType):
		for num in range(0, self.numSeqs): 
			randomLength = numpy.random.randint(self.minLen, self.maxLen); 
			randomSeq = self.GetRandomSequence(randomLength); 
			if (setType == "positive"):
				self.positiveSet.insert(num, randomSeq);
			else:
				self.negativeSet.insert(num, randomSeq);

	def embedMotifInSequence(self):
		motif = self.motif;
		motifLocation = self.motifLocation;
		print "Motif to embed: ", self.motif, ", at location: ", motifLocation;
		for idx, seq in enumerate(self.positiveSet):
			if motifLocation == 'random':
				location = numpy.random.randint(0, len(seq) - len(motif));
				new_seq = seq[:location] + motif + seq[location:]
				self.positiveSet[idx] = new_seq;
