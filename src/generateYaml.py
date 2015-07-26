import yaml


seqMinLen = [25, 50]
seqMaxLen = [55, 100]
aPercent = [10, 20, 30, 40]
tPercent = [10, 20, 30, 40]
gPercent = [10, 20, 30, 40]
cPercent = [10, 20, 30, 40]
MotifType = ['HuR', 'TTP', 'Generated']
HuRMotif = ['ATTTA', 'CTTTTTC']
TtpMotif = ['TTATTTATT']
NumMotifs = [1, 2]
DistanceMotifs = [20, 70, 100]
MotifLocation = ['random', -5, 10]

numberSeq = 10;

aPercent = 30;
tPercent = 30;
gPercent = 20;
cPercent = 20;

def generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent, gPercent, cPercent, motifType, actualMotif):
	data = dict(
		 sequence = dict (
			  minLen = minVal,
			  maxLen = maxVal,
			  numSeq = numberSeq,
			  A = aPercent,
			  T = tPercent,
			  G = gPercent,
			  C = cPercent,
		  ),
		  motif = dict (
			type = motifType,
			length = len(actualMotif),
			numMotifs = 1,
			distance = 0,
			location = 'random',
			motif = actualMotif
		  )
		 )
	yamlFileName = location + "/seq_" + idx + ".yml";
	with open(yamlFileName, 'w') as outfile:
		 outfile.write( yaml.dump(data, default_flow_style=True) )

## Program stats here.
def CreateConfFiles(location):
	for idx, minVal in enumerate(seqMinLen):
		for maxVal in seqMaxLen:
			for motifType in MotifType:
				if motifType == 'HuR':
					for hurMotif in HuRMotif:
						idx = str(minVal) + "_" + str(maxVal) + "_" + hurMotif;
						generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
										 gPercent, cPercent, motifType, hurMotif);
				elif motifType == 'TTP':
					for ttpMotif in TtpMotif:
						idx = str(minVal) + "_" + str(maxVal) + "_" + ttpMotif;
						generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
										 gPercent, cPercent, motifType, ttpMotif);
				else:
					motifType = 'Generated'
					idx = str(minVal) + "_" + str(maxVal);
					generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
										 gPercent, cPercent, motifType, "AAATTTGGGCCC");

#CreateConfFiles("/tmp");