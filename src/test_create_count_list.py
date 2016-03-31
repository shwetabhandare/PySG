import unittest
import parseKspectrum
import splitKmerInDict

class TestCreateCountList(unittest.TestCase):
	def testGetInitializedLists(self):
		aList, tList, gList, cList = splitKmerInDict.getInitializedLists(10);
		self.assertEqual(len(aList), 10);
		self.assertEqual(len(tList), 10);
		self.assertEqual(len(gList), 10);
		self.assertEqual(len(cList), 10);
		self.assertTrue(all(el==0 for el in aList)); 
		self.assertTrue(all(el==0 for el in tList)); 
		self.assertTrue(all(el==0 for el in gList)); 
		self.assertTrue(all(el==0 for el in cList)); 
		
	def testCreateKmerDictList(self):
		kmerList = ['TTAATATTT', 'TTTTTATTT']
		maxKmerLen, kmerDictList = splitKmerInDict.createKmerDictList(kmerList)
		self.assertEqual(maxKmerLen, 9);
		self.assertEqual(len(kmerDictList), 2);
		kmerDict1 = {0: 'T', 1: 'T', 2: 'A', 3: 'A', 4: 'T', 5: 'A', 6: 'T', 7: 'T', 8: 'T'} 
		kmerDict2 = {0: 'T', 1: 'T', 2: 'T', 3: 'T', 4: 'T', 5: 'A', 6: 'T', 7: 'T', 8: 'T'}


		shared_items = set(kmerDict1.items()) & set(kmerDictList[0].items())
		self.assertEqual(len(shared_items), len(kmerDict1))

		shared_items = set(kmerDict2.items()) & set(kmerDictList[1].items())
		self.assertEqual(len(shared_items), len(kmerDict2))

	def testGetUpdateCounts(self):
		kmerList = ['TTAATATTT', 'TTTTTATTT']
		maxKmerLen, kmerDictList = splitKmerInDict.createKmerDictList(kmerList)
		aList, tList, gList, cList = splitKmerInDict.getUpdatedCounts(maxKmerLen, kmerDictList);
		expectedAList = [ 0, 0, 1, 1, 0, 2, 0, 0, 0]
		expectedTList = [2, 2, 1, 1, 2, 0, 2, 2, 2]

		self.assertTrue(len(set(expectedAList).intersection(aList)) > 0)
		self.assertTrue(len(set(expectedTList).intersection(tList)) > 0)


	def testGetUpdateCountsLargeKmerList(self):
		kmerList = ['TTTTATTAT', 'TTCTATTTA', 'ATTTATTTT', 'TTAATATTT', 'TTTATTTTT', 'TAATTATTT', 'TTATTTTTA', 'TATTATTTT', 'TATTTATTA', 'CTTATTTAT', 'ATTTTATTT', 'TTTTTATTT', 'TATATTTAT', 'TATTTTATT']
		maxKmerLen, kmerDictList = splitKmerInDict.createKmerDictList(kmerList)
		aList, tList, gList, cList = splitKmerInDict.getUpdatedCounts(maxKmerLen, kmerDictList);
		expectedAList = [2, 5, 3, 4, 4, 5, 1, 3, 3] 
		expectedTList = [11, 9, 10, 10, 10, 9, 13, 11, 11] 
		expectedGList = [0, 0, 0, 0, 0, 0, 0, 0, 0] 
		expectedCList = [1, 0, 1, 0, 0, 0, 0, 0, 0]

		self.assertTrue(len(set(expectedAList).intersection(aList)) > 0)
		self.assertTrue(len(set(expectedTList).intersection(tList)) > 0)
		self.assertTrue(len(set(expectedGList).intersection(gList)) > 0)
		self.assertTrue(len(set(expectedCList).intersection(cList)) > 0)

	def testCreateNormalizedList(self):
		kmerList = ['TTTTATTAT', 'TTCTATTTA', 'ATTTATTTT', 'TTAATATTT', 'TTTATTTTT', 'TAATTATTT', 'TTATTTTTA', 'TATTATTTT', 'TATTTATTA', 'CTTATTTAT', 'ATTTTATTT', 'TTTTTATTT', 'TATATTTAT', 'TATTTTATT']
		aList = [2, 5, 3, 4, 4, 5, 1, 3, 3] 
		tList = [11, 9, 10, 10, 10, 9, 13, 11, 11] 
		gList = [0, 0, 0, 0, 0, 0, 0, 0, 0] 
		cList = [1, 0, 1, 0, 0, 0, 0, 0, 0]

		expectedNormalizedAList = [0.143, 0.357, 0.214, 0.286, 0.286, 0.357, 0.071, 0.214, 0.214] 
		expectedNormalizedTList = [0.786, 0.643, 0.714, 0.714, 0.714, 0.643, 0.929, 0.786, 0.786] 
		expectedNormalizedGList = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] 
		expectedNormalizedCList = [0.071, 0.0, 0.071, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

		aList, tList, gList, cList = splitKmerInDict.createNormalizedLists(aList, tList, gList, cList, len(kmerList));

		self.assertTrue(len(set(expectedNormalizedAList).intersection(aList)) > 0)
		self.assertTrue(len(set(expectedNormalizedTList).intersection(tList)) > 0)
		self.assertTrue(len(set(expectedNormalizedGList).intersection(gList)) > 0)
		self.assertTrue(len(set(expectedNormalizedCList).intersection(cList)) > 0)
