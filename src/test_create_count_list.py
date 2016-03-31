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

