import yaml


seqMinLen = [19]
seqMaxLen = [343]
aPercent = [35]
tPercent = [35]
gPercent = [15]
cPercent = [15]
MotifType = ['']
#MotifType = ['HuR']
#aPercent = [10, 20, 30, 40]
#tPercent = [10, 20, 30, 40]
#gPercent = [10, 20, 30, 40]
#cPercent = [10, 20,s 30, 40]
#aPercent = [10, 20, 30, 40]
#tPercent = [10, 20, 30, 40]
#gPercent = [10, 20, 30, 40]
#cPercent = [10, 20, 30, 40]
#MotifType = ['HuR', 'TTP', 'Generated']
HuRMotif = ['ATTTA', 'CTTTTTC']
TtpMotif = ['TTATTTATT']
NumMotifs = [1]
DistanceMotifs = [20, 70, 100]
MotifLocation = ['random', -5, 10]

numberSeq = 3642;

aPercent = 30;
tPercent = 30;
gPercent = 20;
cPercent = 20;

def generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent, gPercent, cPercent, motifType, actualMotif, numMotifs, distance):
	if numMotifs == 1:
		distance = 0;

	if motifType != "":

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
			numMotifs = numMotifs,
			distance = distance,
			location = 'random',
			motif = actualMotif,
		 )
		)
	else:
		data = dict(
		 sequence = dict (
			  minLen = minVal,
			  maxLen = maxVal,
			  numSeq = numberSeq,
			  A = aPercent,
			  T = tPercent,
			  G = gPercent,
			  C = cPercent,
		  )
		)		
	yamlFileName = location + "/seq_" + str(idx) + ".yml";
	print "Data: ", data
	with open(yamlFileName, 'w') as outfile:
		 outfile.write( yaml.dump(data, default_flow_style=True) )

## Program stats here.
def CreateConfFiles(location):
	for idx, minVal in enumerate(seqMinLen):
		for maxVal in seqMaxLen:
			for numMotifs in NumMotifs:
				for distance in DistanceMotifs:
					for motifType in MotifType:
						print "MotifType: ", motifType
						if motifType == 'HuR':
							for hurMotif in HuRMotif:
								idx = str(minVal) + "_" + str(maxVal) + "_" + hurMotif;
								generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
												 gPercent, cPercent, motifType, hurMotif, numMotifs, distance);
						elif motifType == 'TTP':
							for ttpMotif in TtpMotif:
								idx = str(minVal) + "_" + str(maxVal) + "_" + ttpMotif;
								generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
												 gPercent, cPercent, motifType, ttpMotif, numMotifs, distance);
						elif motifType == 'Generated':
							idx = str(minVal) + "_" + str(maxVal);
							generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
												 gPercent, cPercent, motifType, "AAATTTGGGCCC", numMotifs, distance);

						else: 
							generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
										 gPercent, cPercent, "", "", 0, 0);
#CreateConfFiles("/tmp");
