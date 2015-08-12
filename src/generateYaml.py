import yaml


seqMinLen = [50, 100]
seqMaxLen = [105, 150]
aPercent = [10, 20, 30, 40]
tPercent = [10, 20, 30, 40]
gPercent = [10, 20, 30, 40]
cPercent = [10, 20, 30, 40]
# generate A, T, G, C using s = np.random.dirichlet((0.1,0.1,0.1,0.1),1)
# expect sparse distribution, there is no apriori information.
# >>> tp = np.random.dirichlet((0.1,0.1,0.1,0.1),1)[0]
#>>> tp
#array([  4.61012266e-01,   5.37356121e-01,   1.59085119e-03,
#         4.07617397e-05])
#>>> hur_tp = np.random.dirichlet(tp,1)[0]
#>>> hur_tp
#array([ 0.06761885,  0.93238115,  0.        ,  0.        ])
#>>> hur_tp = np.random.dirichlet(tp*100,1)[0]
#>>> hur_tp
#array([  4.13762328e-01,   5.85305180e-01,   9.32492390e-04,
#         1.97868086e-59])
#>>> ttp_tp = np.random.dirichlet(tp*100,1)[0]
#>>> ttp_tp
#array([  5.14829673e-01,   4.85166046e-01,   4.28061017e-06,
#         1.59080129e-58])
# Feed hur_tp, ttp_tp to the sequence generator for A, T, G, C percents.
# when gamma is large, the ATGC percentages will be similar, and when its smaller they will differ more.
# here gammma = 100.

MotifType = ['HuR', 'TTP', 'Generated']
HuRMotif = ['ATTTA', 'CTTTTTC']
TtpMotif = ['TTATTTATT']
NumMotifs = [1]
DistanceMotifs = [20, 70, 100]
MotifLocation = ['random', -5, 10]

numberSeq = 10;

aPercent = 30;
tPercent = 30;
gPercent = 20;
cPercent = 20;

def generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent, gPercent, cPercent, motifType, actualMotif, numMotifs, distance):
	if numMotifs == 1:
		distance = 0;

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
			for numMotifs in NumMotifs:
				for distance in DistanceMotifs:
					for motifType in MotifType:
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
						else:
							motifType = 'Generated'
							idx = str(minVal) + "_" + str(maxVal);
							generateYaml(idx, location, minVal, maxVal, numberSeq, aPercent, tPercent,
												 gPercent, cPercent, motifType, "AAATTTGGGCCC", numMotifs, distance);

#CreateConfFiles("/tmp");
