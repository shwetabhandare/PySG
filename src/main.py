import Conf;
import SeqGen;
from Conf import *;
from SeqGen import *;


conf = Conf("../conf/SeqGen.yaml");
seqGen = SeqGen(conf);
seqGen.GenerateRandomSequences("negative");
seqGen.GenerateRandomSequences("positive");

motif = seqGen.GenerateMotif()
print "Motif: ", motif;

seqGen.embedMotifInSequence();

for seq in seqGen.GetPositiveSet():
	print seq;

