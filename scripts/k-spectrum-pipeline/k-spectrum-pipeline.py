import sys
import yaml
import SplitTrainAndTest;
from SplitTrainAndTest import *
import seq_counter
from seq_counter import *

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

CreateCombinedFile(posTrainFileName, negTrainFileName)
CreateCombinedFile(posTestFileName, negTestFileName);

## Return the name of the file that was created by CreateCombinedFile()
## Get the number of training examples from the filename : HuR_Train_NonHuR_Train_3277_3277_CompleteSet.txt

#trainCombinedFile, numTrain = CreateCombinedFile(posTrainFileName, negTrainFileName)
#testCombinedFile, numTest = CreateCombinedFile(posTrainFileName, negTrainFileName)


## Use the training file, and number of training examples to create a model.
## Use the model generated, and the test file
