import Conf;
import SeqGen;
from Conf import *;
from SeqGen import *;


conf = Conf("../conf/SeqGen.yaml");
seqGen = SeqGen(conf);
seqGen.GenerateRandomSequences();

