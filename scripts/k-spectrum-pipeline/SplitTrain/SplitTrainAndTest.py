import fasta
import sys
import random
import re
import os


def getAllData(dataFile): 
	allData = []
	for record in fasta.fasta_itr(dataFile):						#Create parallel arrays for negative sequence/headers
		sequence = re.sub('[()\'\']', '', record.sequence.strip());
		header = re.sub('[()\'\']', '', record.header.strip());
		allData.append([header, sequence]);

	return allData;

def getTrainIds(allData, numTrain):
	trainIdx = [];
	for item in allData:
		index = random.randint(0, len(allData) - 1);
		rc = index in trainIdx;
		while (rc != False):
			index = random.randint(0, len(allData) - 1);
			rc = index in trainIdx;
		if (len(trainIdx) < numTrain):
			trainIdx.append(index)
	return trainIdx;


def getTestIds(trainIds, allData):
	testIds = []
	for i in range(0, len(allData)):
		index = random.randint(0, len(allData) - 1);
		rc = index in trainIds;
		while (rc != False):
			index = random.randint(0, len(allData) - 1);
			rc = index in trainIds;
		rc1 = index in testIds;
		if (rc == False and rc1 == False):
			testIds.append(index);

	return testIds;

def getTrainAndTestData(allData, trainPercent, testPercent):
	testData = [];
	trainData = [];
	numTrain = (len(allData) * trainPercent)/100;	
	trainIds = getTrainIds(allData, numTrain);
	testIds = getTestIds(trainIds, allData);
	numAdded = 0;
	for i in range(0, len(trainIds)):
		indexToAdd = trainIds[i];
		trainData.append(allData[indexToAdd]);
	for i in range(0, len(testIds)):
		indexToAdd = testIds[i];
		testData.append(allData[indexToAdd]);
	return trainData, testData;
			

def writeDataToFile(data, outputFile):
	outFile = open(outputFile , "w")
	for item in data:
		outFile.write(">")
		outFile.write(item[0]);
		outFile.write("\n")
		outFile.write(item[1]);
		outFile.write("\n")
	outFile.close();

# main starts here

#Prepare and open files
confFile = sys.argv[1]; // yaml file.
confMap = yaml.load(confFile);

posDataFile = confMap["input"]["posFile"]
negDataFile = confMap["input"]["negFile"]

trainPercent = int(confMap["input"]["trainPercent"]);
testPercent = 100 - trainPercent;

posTrainFileName = confMap["output"]["trainPosFile"]
negTrainFileName = confMap["output"]["trainNegFile"]

posTestFileName = confMap["output"]["testPosFile"]
negTestFileName = confMap["output"]["testNegFile"]



allData = getAllData(dataFile);
trainData, testData = getTrainAndTestData(allData, trainPercent, testPercent);
print len(trainData), len(testData)
writeDataToFile(trainData, trainFileName);
writeDataToFile(testData, testFileName);
