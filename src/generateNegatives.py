import Conf;
import SeqGen;
from Conf import *;
from SeqGen import *;
import generateYaml;
from generateYaml import *;
import sys;

directory = sys.argv[1]
filename = os.path.splitext(os.path.basename(confFile))[0]
negFastaFile = directory + "/" + filename + "_neg.fa";
seqGen.SetNegFileName(negFastaFile)
seqGen = SeqGen(conf);
seqGen.GenerateRandomSequences("negative", 0);
