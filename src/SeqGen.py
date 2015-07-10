class SeqGen:
	def __init__(self, numSeqs, minLen, maxLen):
		self.numSeqs = numSeqs;
		self.minLen = minLen;
		self.maxLen = maxLen;



	def GenerateRandomSequences():
		for num in range(0, self.numSeqs):
			randomLength = randint(self.minLen, self.maxLen);

