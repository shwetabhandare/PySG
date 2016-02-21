import yaml

#numSeq = [1000, 2000, 3000, 4000, 5000]
#seqLen = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

numSeq = [1000, 2000, 3000, 4000, 5000]
seqLen = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
# generate A, T, G, C using s = np.random.dirichlet((0.1,0.1,0.1,0.1),1)
# expect sparse distribution, there is no apriori information.
# tp = np.random.dirichlet((0.1,0.1,0.1,0.1),1)[0]
# tp
#array([  4.61012266e-01,   5.37356121e-01,   1.59085119e-03,
#         4.07617397e-05])
# hur_tp = np.random.dirichlet(tp,1)[0]
# hur_tp
#array([ 0.06761885,  0.93238115,  0.        ,  0.        ])
# hur_tp = np.random.dirichlet(tp*100,1)[0]
# hur_tp
#array([  4.13762328e-01,   5.85305180e-01,   9.32492390e-04,
#         1.97868086e-59])
# ttp_tp = np.random.dirichlet(tp*100,1)[0]
# ttp_tp
#array([  5.14829673e-01,   4.85166046e-01,   4.28061017e-06,
#         1.59080129e-58])
# Feed hur_tp, ttp_tp to the sequence generator for A, T, G, C percents.
# when gamma is large, the ATGC percentages will be similar, and when its smaller they will differ more.
# here gammma = 100.


def generateYaml(location, numberSeq, seqLength, inpFastaFile):
	outFileName = "NoSignal_" + str(numberSeq) + "_" + str(seqLength)
	print outFileName
	data = dict(
	 sequence = dict (
		nosignal = dict(
		  seqLen = seqLength,
		  numSeq = numberSeq,
		  fastaFile = inpFastaFile,
		  outFastaFile = location + "/" + outFileName + ".fa",
		  )
		)
	)
	yamlFileName = location + "/" + outFileName + ".yml";
	print "Data: ", data
	with open(yamlFileName, 'w') as outfile:
		 outfile.write( yaml.dump(data, default_flow_style=True) )

## Program stats here.
def CreateConfFiles(location):
	inpFastaFile = "/projects/bhandare/workspace/scripts/NegFileCreator/3UTR_transcripts_Human.txt"
	for numSeqIdx, i in enumerate(numSeq):
		for numSeqLenIdx, j in enumerate(seqLen):
			print i, j
			generateYaml(location, i, j, inpFastaFile);


if __name__ == "__main__":
	import sys
	location = sys.argv[1]
	CreateConfFiles(location)
