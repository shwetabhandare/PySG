def getPredictedKmerDict(dremeResultFile):
	return getPredictedDremeMotifs(dremeResultFile);

def getReadKmerDict(realKmersCsvFile):
	return GetRealKmerDict(realKmersCsvFile)


def getKmerLengthFromREString(kmerREString):
	kmers = kmerREString.split('|')
	kmerStr = "".join([str(x) for x in kmers] )
	kmerStr = kmerStr[1:-1]
	return len(kmerStr);

def getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer):
	numTP = 0;
	numFP = 0;

	predKmerIndex = 0;
	#print "Start Index: ", str(startIndex), ", End Index: ", str(endIndex)
	for index in range(startIndex, endIndex):
		if seq[index] == predKmer[predKmerIndex]:
			numTP = numTP + 1;
		else:
			numFP = numFP + 1;
		predKmerIndex = predKmerIndex + 1;
	#print "After KMER: TP: ", str(numTP), ", FP: " , str(numFP)

	return numTP, numFP;

def getPredictedInfo(m):
	return int(m.start()), int(m.end()), m.group(1)

def getNumbersForKmerNotAddedToSeq(kmerREString, seq):
	numFP = 0;

	for m in re.finditer(kmerREString, seq):
		predictedStart, predictedEnd, predKmer  = getPredictedInfo(m);
		numFP = numFP + (predictedEnd - predictedStart);
	return numFP;

def getStartIndexAndUpdateNumbers(predictedStart, predictedEnd, realStart, realEnd, predKmer):
	numFP = 0;
	numFN = 0;

	if predictedStart > realStart:
		if predictedStart > realEnd:
			numFP = numFP + len(predKmer)
			numFN = numFN + (realEnd - realStart)
			#There is no k-mer to compare, so we set startIndex to realEnd
			# so that we never compare the k-mer. 
			startIndex = realEnd;
		else:
			startIndex = predictedStart;
			numFN = predictedStart - realStart;
	else:
		#Predicted k-mer starts before the real k-mer.
		if predictedEnd < realStart:
			#Predicted k-mer ends before the real k-mer starts,
			# There is no real k-mer to compare.
			startIndex = realStart;
			numFP = numFP + len(predKmer)
			# We can't update the numFN here as there might be more k-mers
			# that are after the real start.
		else:
			#Predicted k-mer ends after the real-kmer starts.
			startIndex = realStart;
			numFP = numFP + (realStart - predictedStart);

	return startIndex, numFP, numFN;

def getNumbersForSeq(kmerREString, realStart, realEnd, seq):
	predictedStart = predictedEnd = 0;
	notFoundFP = 0;
	notFoundFN = 0;
	numTP = 0;
	numFN = 0;
	numFP = 0;

	for m in re.finditer(kmerREString, seq):

		predictedStart, predictedEnd, predKmer  = getPredictedInfo(m);

		#print "Predicted Kmer: ", predKmer, ", Predicted Start: ", str(predictedStart), ", End: ", str(predictedEnd)
		startIndex, startFP, startFN = getStartIndexAndUpdateNumbers(predictedStart, predictedEnd, realStart, realEnd, predKmer);

		numFN = numFN + startFN;
		numFP = numFP + startFP;

		#print "Based on Start Location: FP: " , str(startFP), ", FN: ", str(startFN)
		endIndex = getEndIndex(predictedEnd, realEnd)	

		# This can happen if the predicted kmer ends before the real kmer starts.
		if endIndex < startIndex:
			endIndex = startIndex;
	
		kmerTP, kmerFP = getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		#print "After KMER TP: ", str(kmerTP), str(kmerFP)
		
		numTP = numTP + kmerTP
		numFP = numFP + kmerFP

		# Update realStart Index after k-mer is parsed.
		if endIndex >= realStart:
			realStart = endIndex;
		
		
	# We did not find the predicted k-mer in the sequence.
	if predictedStart == 0 and predictedEnd == 0:
		notFoundFP = notFoundFP + getKmerLengthFromREString(kmerREString)

	# We did not find any predicted k-mers in the sequence. 
	if realStart < realEnd:
		numFN = numFN + (realEnd - realStart);

	numFP = numFP + notFoundFP
	numFN = numFN + notFoundFN

	return numTP, numFP, numFN;

def getEndIndex(predictedEnd, realEnd):
	if predictedEnd < realEnd:
		endIndex = predictedEnd
	else:
		endIndex = realEnd;

	return endIndex;

def getRealKmerDetails(realKmerDict, seqid):
	if seqid in realKmerDict:
		return realKmerDict[seqid][0], int(realKmerDict[seqid][1]), (int(realKmerDict[seqid][1]) + len(realKmerDict[seqid][0]))
	else:
		return "", 0, 0

def getTotalsForSeq(realKmer, realStart, realEnd, kmerREString, seq):
	numTP = numFP = numFN = 0;
	
	if realKmer == "" and realStart == 0 and realEnd == 0:
		print "Sequence does not contain kmer."
		numFP = getNumbersForKmerNotAddedToSeq(kmerREString, seq);
	else:
		numTP, numFP, numFN = getNumbersForSeq(kmerREString, realStart, realEnd, seq);

	return numTP, numFP, numFN;

def getTotalNumbersForSeqDict(SeqDict, realKmerDict, pssmList, pwm, predictedMotifs):

	numTotalTP = numTotalFP = numTotalFN = 0;

	for seqid, seq in SeqDict.iteritems():
		realKmer, realStart, realEnd = getRealKmerDetails(realKmerDict, seqid);
		
		if predictedMotifs != None:
			kmerREString = getKmerRE(predictedMotifs, seq)
		elif pssmList != None:
			kmerREString = getKmerFromPSSM(pssmList, seq)
		else:
			kmerREString = getKmerFromPWM(pwm, seq)

		print "SEQ ID: ", seqid, "REAL KMER: ", realKmer, ", KMER RE: ", kmerREString;
		numTP, numFP, numFN = getTotalsForSeq(realKmer, realStart, realEnd, kmerREString, seq);
		print seqid, ": TP: ", str(numTP), ", FP: ", str(numFP), ", FN: ", str(numFN)

		numTotalTP = numTotalTP + numTP;
		numTotalFP = numTotalFP + numFP;
		numTotalFN = numTotalFN + numFN;
	
	return numTotalTP, numTotalFP, numTotalFN;

def GetTotalNumbers(realKmerDict, posFile, negFile, pssmList, pwm, predictedMotifs):
	PosSeqDict = SeqGenUtils.fasta_read(posFile);
	NegSeqDict = SeqGenUtils.fasta_read(negFile);

	numPosTP = 	numPosFP = 	numPosFN = 0;	
	numNegTP = 	numNegFP = 	numNegFP = 0;	

	numPosTP, numPosFP, numPosFN = getTotalNumbersForSeqDict(PosSeqDict, realKmerDict, pssmList, pwm, predictedMotifs);
	numNegTP, numNegFP, numNegFN = getTotalNumbersForSeqDict(NegSeqDict, realKmerDict, pssmList, pwm, predictedMotifs);

	print "Pos File: TP: ", str(numPosTP), ", FP: ", numPosFP, ", FN: ", numPosFN
	print "Neg File: TP: ", str(numNegTP), ", FP: ", numNegFP, ", FN: ", numNegFN
	return (numPosTP + numNegTP), (numPosFP + numNegFP) , (numPosFN + numNegFN);


def getTotalFNKmerNotFound(realKmerDict):
	totalFN = 0
	for seqid, kmerDetails in realKmerDict.iteritems():
		kmer = realKmerDict[seqid][0]
		totalFN = totalFN + len(kmer)
	return totalFN;