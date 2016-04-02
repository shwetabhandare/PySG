import unittest
import compareKmers
import parseDreme
import TAMO_Motif
import splitKmerInDict

class TestCompareKmers(unittest.TestCase):
	def testPredictedStartEndBeforeRealStart(self):
		predictedStart = 0;
		predictedEnd = 10;
		realStart = 15;
		realEnd = 20;
		predKmer = "ATTTA"

		startIndex, numFP, numFN = compareKmers.getStartIndexAndUpdateNumbers(predictedStart, 
			predictedEnd, realStart, realEnd, predKmer)
		self.assertEqual(startIndex, realStart)
		self.assertEqual(numFP, 5)
		self.assertEqual(numFN, 0)

	def testPredictedStartEndAfterRealStart(self):
		predictedStart = 25;
		predictedEnd = 35;
		realStart = 15;
		realEnd = 20;
		predKmer = "ATTTA"

		startIndex, numFP, numFN = compareKmers.getStartIndexAndUpdateNumbers(predictedStart, 
			predictedEnd, realStart, realEnd, predKmer)
		self.assertEqual(startIndex, realEnd)
		self.assertEqual(numFP, 5)
		self.assertEqual(numFN, 5)		

	def testPredictedStartBeforeRealStartButEndBeforeRealEnd(self):
		predictedStart = 0;
		predictedEnd = 20;
		realStart = 10;
		realEnd = 30;
		predKmer = "ATTTA"

		startIndex, numFP, numFN = compareKmers.getStartIndexAndUpdateNumbers(predictedStart, 
			predictedEnd, realStart, realEnd, predKmer)
		self.assertEqual(startIndex, realStart)
		self.assertEqual(numFP, (realStart - predictedStart))
		self.assertEqual(numFN, 0)			

	def testPredictedStartBeforeRealStartButEndAfterRealEnd(self):
		predictedStart = 0;
		predictedEnd = 35;
		realStart = 10;
		realEnd = 30;
		predKmer = "ATTTA"

		startIndex, numFP, numFN = compareKmers.getStartIndexAndUpdateNumbers(predictedStart, 
			predictedEnd, realStart, realEnd, predKmer)
		self.assertEqual(startIndex, realStart)
		self.assertEqual(numFP, (realStart - predictedStart))
		self.assertEqual(numFN, 0)	

	def testPredictedStartAfterRealStartButEndBeforeRealEnd(self):
		predictedStart = 15;
		predictedEnd = 25;
		realStart = 10;
		realEnd = 30;
		predKmer = "ATTTA"

		startIndex, numFP, numFN = compareKmers.getStartIndexAndUpdateNumbers(predictedStart, 
			predictedEnd, realStart, realEnd, predKmer)
		self.assertEqual(startIndex, predictedStart)
		self.assertEqual(numFN, (predictedStart - realStart))
		self.assertEqual(numFP, 0)		

	def testPredictedStartAfterRealStartButEndAfterRealEnd(self):
		predictedStart = 15;
		predictedEnd = 35;
		realStart = 10;
		realEnd = 30;
		predKmer = "ATTTA"

		startIndex, numFP, numFN = compareKmers.getStartIndexAndUpdateNumbers(predictedStart, 
			predictedEnd, realStart, realEnd, predKmer)
		self.assertEqual(startIndex, predictedStart)
		self.assertEqual(numFN, (predictedStart - realStart))
		self.assertEqual(numFP, 0)		

	def testPredictedEndBeforeRealEnd(self):
		predictedEnd = 10;
		realEnd = 20;

		endIndex = compareKmers.getEndIndex(predictedEnd, realEnd);
		self.assertEqual(endIndex, predictedEnd)		
	
	def testPredictedEndEqualRealEnd(self):
		predictedEnd = 30;
		realEnd = 30;

		endIndex = compareKmers.getEndIndex(predictedEnd, realEnd);
		self.assertEqual(endIndex, realEnd)	

	def testGetRealKmerDetails(self):
		realKmerDict = dict()

		realKmerDict["1"] = ["ATTA", 10]
		realKmerDict["2"] = ["ATTAAAA", 20]
		realKmerDict["3"] = ["ATTAAAAATTTT", 30]

		realKmer, realStart, realEnd = compareKmers.getRealKmerDetails(realKmerDict, "1")

		self.assertEqual(realKmer, "ATTA")
		self.assertEqual(realStart, 10)
		self.assertEqual(realEnd, (10 + len("ATTA")))

		realKmer, realStart, realEnd = compareKmers.getRealKmerDetails(realKmerDict, "2")

		self.assertEqual(realKmer, "ATTAAAA")
		self.assertEqual(realStart, 20)
		self.assertEqual(realEnd, (20 + len("ATTAAAA")))

		realKmer, realStart, realEnd = compareKmers.getRealKmerDetails(realKmerDict, "3")


		self.assertEqual(realKmer, "ATTAAAAATTTT")
		self.assertEqual(realStart, 30)
		self.assertEqual(realEnd, (30 + len("ATTAAAAATTTT")))

	def testPredictedEndAfterRealEnd(self):
		predictedEnd = 30;
		realEnd = 20;

		endIndex = compareKmers.getEndIndex(predictedEnd, realEnd);
		self.assertEqual(endIndex, realEnd)	


	def testKmerComparisonExactKmer(self):
		seq = "GGAACCGCGTTCGGGGGGGGGGGGGCCGAACCCTTCCAGCATTGAGCTCCTGCCGCTAGCTTATGCGGCCTCCCATCCAGTCGGCCGAGACGCACGACTT"
		predKmer = "GAACCCTTC"

		startIndex = seq.find(predKmer)
		endIndex = startIndex + len(predKmer)

		numTP, numFP = compareKmers.getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		self.assertEqual(numTP, len(predKmer))
		self.assertEqual(numFP, 0)		

	def testKmerComparisonKmerOffByOneCharacterEnd(self):
		seq = "GGAACCGCGTTCGGGGGGGGGGGGGCCGAACCCTTCCAGCATTGAGCTCCTGCCGCTAGCTTATGCGGCCTCCCATCCAGTCGGCCGAGACGCACGACTT"
		predKmer = "GAACCCTTG"

		startIndex = 27
		endIndex = startIndex + len(predKmer)

		numTP, numFP = compareKmers.getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		self.assertEqual(numTP, len(predKmer) - 1)
		self.assertEqual(numFP, 1)		

	def testKmerComparisonKmerOffByOneCharacterStart(self):
		seq = "GGAACCGCGTTCGGGGGGGGGGGGGCCGAACCCTTCCAGCATTGAGCTCCTGCCGCTAGCTTATGCGGCCTCCCATCCAGTCGGCCGAGACGCACGACTT"
		predKmer = "CAACCCTTC"

		startIndex = 27
		endIndex = startIndex + len(predKmer)

		numTP, numFP = compareKmers.getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		self.assertEqual(numTP, len(predKmer) - 1)
		self.assertEqual(numFP, 1)		

	def testKmerComparisonKmerOffMiddle(self):
		seq = "GGAACCGCGTTCGGGGGGGGGGGGGCCGAACCCTTCCAGCATTGAGCTCCTGCCGCTAGCTTATGCGGCCTCCCATCCAGTCGGCCGAGACGCACGACTT"
		predKmer = "GAATTTTTC"

		startIndex = 27
		endIndex = startIndex + len(predKmer)

		numTP, numFP = compareKmers.getNumbersAfterKmerComparison(startIndex, endIndex, seq, predKmer)
		self.assertEqual(numTP, len(predKmer) - 3)
		self.assertEqual(numFP, 3)		

	def test_kmerREString(self):
		seq = "GGAACCGCGTTCGGGGGGGGGGGGGCCGAACCCTTCCAGCATTGAGCTCCTGCCGCTAGCTTATGCGGCCTCCCATCCAGTCGGCCGAGACGCACGACTT"
		predictedMotifs = ['MCCCGA', 'AAAAAA']
		kmerReString = compareKmers.getKmerRE(predictedMotifs, seq)
		self.assertEqual(kmerReString, '(TCGGGG|CTTATG)')

	def test_getNumbersForSeqOneKmer(self):
		seq = "GGAACCGCXXXXXXGGGGGCGACXXXXXXGGGCATTGAGCTC"
		predKmer = "XXXXXX"

		realStart = 4;
		realEnd = 44
		kmerREString = "(XXXXXX)"

		numTP, numFP, numFN = compareKmers.getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		self.assertEqual(numTP, 12)
		self.assertEqual(numFP, 0)
		self.assertEqual(numFN, 28)


	def test_getNumbersForSeqTwoKmers(self):
		seq = "GGAACCGCXXXXXXGGGGGCGACYYYYYYGGGCATTGAGCTC"
		realStart = 4;
		realEnd = 44

		kmerREString = "(XXXXXX|YYYYYY)"

		numTP, numFP, numFN = compareKmers.getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		self.assertEqual(numTP, 12)
		self.assertEqual(numFP, 0)
		self.assertEqual(numFN, 28)

	def test_getLengthOfKmersFromKmerREString(self):
		kmerREString = "(XXXXXX|YYYYYY)"

		kmerStringLen = compareKmers.getKmerLengthFromREString(kmerREString);
		self.assertEqual(kmerStringLen, len(kmerREString) - 3)

	def test_getNumbersForSeqKmerNotFound(self):
		seq = "XXXXXXGCCCCCXXXXXXGGGGGCGACYYYYYYGGGCATTGAGCTC"

		realStart = 4;
		realEnd = 44

		kmerREString = "(ABCDEF|XYZFFF)"

		numTP, numFP, numFN = compareKmers.getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		self.assertEqual(numTP, 0)
		self.assertEqual(numFP, 12)
		self.assertEqual(numFN, 40)

	def test_getNumbersForSeqKmerInNegFile(self):
		seq = "GTCTTACAAGAGCCCCGACGCCCCGGCCACTGCGCGCGACTAGCCCTATGTCAGGAAAAAACGCGCACAATGTCCTCCTGCAGGACAGTCGGCTGGCGCTACCGATACGGAT"

		realStart = 0;
		realEnd = 0;

		kmerREString = "(CCCCGA|AAAAAA)"

		numTP, numFP, numFN = compareKmers.getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		
		#self.assertEqual(numTP, 0)
		#self.assertEqual(numFP, 12)
		#self.assertEqual(numFN, 40)

	def test_getNumbersForSeqOneKmerInOneKmerOut(self):
		seq = "ATCCCTAACTCCGGCAAAAAAAAAACCGGAAACTACATCGCTCTCCACCGGTGCAGACGTCGCCTCGCGCCCCGAAACCGGTGCTGGCAGGGTACGTAAT"
		kmerREString = "(CCCCGA|AAAAAA)"
		realKmer = "CCGGCAAAAAAAAAACCGG"
		realStart = 10
		realEnd = realStart + len(realKmer)

		numTP, numFP, numFN = compareKmers.getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		self.assertEqual(numTP, 6)
		self.assertEqual(numFP, 6)
		self.assertEqual(numFN, 13)
	
	def test_getNumbersForSeqTwoKmersInSeq(self):
		seq = "CTGTCCCTTTTCGGGTTTTTTTTTTCCGAGCGGCCTCGGTGGGTGAAATGAACGACACTCATGCGAGCGACACTAGGGCGCCGTTCGTTCTGTGCACCCA"
		kmerREString = "(TCGGGT|TTTTTT)"
		realKmer = "TCGGGTTTTTTTTTTCCGA"
		realStart = 10;
		realEnd = realStart + len(realKmer)
		numTP, numFP, numFN = compareKmers.getNumbersForSeq(kmerREString, realStart, realEnd, seq);
		self.assertEqual(numTP, 12)
		self.assertEqual(numFP, 0)
		self.assertEqual(numFN, 7)
	
	def test_getKmerREFromPSSM(self):
		predictedDremeFile = "/projects/bhandare/workspace/PySG/src/resources/dreme.txt"
		seq = "CTGTCCCTTTTCGGGTTTTTTTTTTCCGAGCGGCCTCGGTGGGTGAAATGAACGACACTCATGCGAGCGACACTAGGGCGCCGTTCGTTCTGTGCACCCA"
		pssmList = parseDreme.getPSSMListFromDremeFile(predictedDremeFile);
		kmerREString = compareKmers.getKmerFromPSSM(pssmList, seq)

		self.assertEqual(kmerREString, '(TCGGGT|TTTTTT)')

	def test_compareKmersTextMotif(self):
		predictedDremeFile = "/projects/bhandare/workspace/PySG/src/resources/dreme.txt"
		realCsvFile = "/projects/bhandare/workspace/PySG/src/resources/Signal.kmers"
		posFile = "/projects/bhandare/workspace/PySG/src/resources/Signal.fa"
		negFile = "/projects/bhandare/workspace/PySG/src/resources/NoSignal.fa"

		numTP, numFP, numFN = compareKmers.CompareDremeKmers(realCsvFile, predictedDremeFile, posFile, negFile, 1);
		self.assertEqual(numTP, 170)
		self.assertEqual(numFP, 304)
		self.assertEqual(numFN, 207)

	def test_compareKmersTextMotifNotAllSeqContainKmer(self):
		predictedDremeFile = "/projects/bhandare/workspace/PySG/src/resources/dreme90.txt"
		realCsvFile = "/projects/bhandare/workspace/PySG/src/resources/Signal90.kmers"
		posFile = "/projects/bhandare/workspace/PySG/src/resources/Signal90.fa"
		negFile = "/projects/bhandare/workspace/PySG/src/resources/NoSignal90.fa"

		numTP, numFP, numFN = compareKmers.CompareDremeKmers(realCsvFile, predictedDremeFile, posFile, negFile, 1);
		self.assertEqual(numTP, 0)
		self.assertEqual(numFP, 490)
		self.assertEqual(numFN, 0)

	def test_compareKmersPSSM(self):
		predictedDremeFile = "/projects/bhandare/workspace/PySG/src/resources/dreme.txt"
		realCsvFile = "/projects/bhandare/workspace/PySG/src/resources/Signal.kmers"
		posFile = "/projects/bhandare/workspace/PySG/src/resources/Signal.fa"
		negFile = "/projects/bhandare/workspace/PySG/src/resources/NoSignal.fa"

		numTP, numFP, numFN = compareKmers.CompareDremeKmers(realCsvFile, predictedDremeFile, posFile, negFile, 0);
		self.assertEqual(numTP, 174)
		self.assertEqual(numFP, 300)
		self.assertEqual(numFN, 206)

	def test_compareKmersPSSMNotAllSeqContainKmer(self):
		predictedDremeFile = "/projects/bhandare/workspace/PySG/src/resources/dreme90.txt"
		realCsvFile = "/projects/bhandare/workspace/PySG/src/resources/Signal90.kmers"
		posFile = "/projects/bhandare/workspace/PySG/src/resources/Signal90.fa"
		negFile = "/projects/bhandare/workspace/PySG/src/resources/NoSignal90.fa"

		numTP, numFP, numFN = compareKmers.CompareDremeKmers(realCsvFile, predictedDremeFile, posFile, negFile, 0);
		self.assertEqual(numTP, 0)
		self.assertEqual(numFP, 490)
		self.assertEqual(numFN, 0)


	def test_getRealKmerDetailsKeyDoesNotExists(self):
		realKmerDict  = dict();

		realKmerDict['0'] = ['ATTTA', 10]
		realKmerDict['1'] = ['ATTTTA', 10]
		realKmerDict['2'] = ['ATTTTTA', 10]
		realKmerDict['3'] = ['ATTTAAAA', 10]

		realKmer, realStart, realEnd = compareKmers.getRealKmerDetails(realKmerDict, '0')
		self.assertEqual(realKmer, 'ATTTA');
		self.assertEqual(realStart, 10);
		self.assertEqual(realEnd, 15)

		# does not exist.
		realKmer, realStart, realEnd = compareKmers.getRealKmerDetails(realKmerDict, '4')
		self.assertEqual(realKmer, '');
		self.assertEqual(realStart, 0);
		self.assertEqual(realEnd, 0)

	def test_getTotalFPWhenPredictedIsEmpty(self):
		realKmerDict  = dict();

		realKmerDict['0'] = ['ATTTA', 10]
		realKmerDict['1'] = ['ATTTTA', 10]
		realKmerDict['2'] = ['ATTTTTA', 10]
		realKmerDict['3'] = ['ATTTAAAA', 10]

		totalFN = compareKmers.getTotalFNKmerNotFound(realKmerDict);
		self.assertEqual(totalFN, 26);
		
	def test_getKmerFromPWM(self):
		predictedKspectrumKmers = "/projects/bhandare/workspace/PySG/src/resources/Signal90_Features.dat"
		pwm = splitKmerInDict.GetKspectrumPWM(predictedKspectrumKmers);
		seq = "GATCTCCCCGTATTTATTTCCTCGACTACCCCCTCTCGCTAAGTTGCAACACAACAACCCGACCCGTTATAACTATGAGAGAAACAAATCGCTCGGACCC"
		kmerRE = compareKmers.getKmerFromPWM(pwm, seq);
		self.assertEqual(kmerRE, "(TATTTATTT)");

	def test_compareKspectrumKmers(self):
		predictedKspectrumKmers = "/projects/bhandare/workspace/PySG/src/resources/Signal90_Features.dat"
		realCsvFile = "/projects/bhandare/workspace/PySG/src/resources/Signal90.kmers"
		posFile = "/projects/bhandare/workspace/PySG/src/resources/Signal.fa"
		negFile = "/projects/bhandare/workspace/PySG/src/resources/NoSignal.fa"

		numTP, numFP, numFN = compareKmers.CompareKspectrumKmers(realCsvFile, predictedKspectrumKmers,  posFile, negFile);

		self.assertEqual(numTP, 50)
		self.assertEqual(numFP, (63 + 180))
		self.assertEqual(numFN, 112)
