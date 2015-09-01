import sys
import yaml
import SplitTrainAndTest;
from SplitTrainAndTest import *

confFile = sys.argv[1]; # yaml file.

f = open(confFile)
confMap = yaml.load(f);
f.close()
print confMap

posDataFile = confMap["posFile"]
negDataFile = confMap["negFile"]

trainPercent = int(confMap["trainPercent"]);
testPercent = 100 - trainPercent;

posTrainFileName = confMap["trainPosFile"]
negTrainFileName = confMap["trainNegFile"]

posTestFileName = confMap["testPosFile"]
negTestFileName = confMap["testNegFile"]

SplitTrainAndTest(posDataFile, negDataFile, trainPercent, posTrainFileName, negTrainFileName, posTestFileName, negTestFileName);

#CreateCombinedFile(posTrainFileName, negTrainFileName)
#CreateCombinedFile(posTestFileName, negTestFileName);
