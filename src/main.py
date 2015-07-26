import Conf;
import SeqGen;
from Conf import *;
from SeqGen import *;
import generateYaml;
from generateYaml import *;


#conf = Conf("../conf/SeqGen.yaml");

generateYaml("/tmp");
conf = Conf("/tmp/seq.yml");
seqGen = SeqGen(conf);
seqGen.GenerateRandomSequences("negative");
seqGen.GenerateRandomSequences("positive");

motif = seqGen.GenerateMotif()

seqGen.embedMotifInSequence();

print "Positive Set: "
for seq in seqGen.GetPositiveSet():
	print seq;

print "Negative Set: "
for seq in seqGen.GetNegativeSet():
	print seq;
