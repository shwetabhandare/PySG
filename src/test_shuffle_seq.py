import unittest
import NoSignal
import Kmer
import os, fnmatch
import yaml
import SeqGenUtils
import parseResults
import uuid;
import Kmer;
import shuffle_utils;

class TestShuffleSeq(unittest.TestCase):
	conFile = "";
	posSeqFile = ""
	negSeqFile = ""

	def setUp(self):
		self.confFile = "/projects/bhandare/workspace/PySG/src/resources/test_shuffle.yml"
		self.posSeqFile = "/projects/bhandare/workspace/PySG/src/resources/tmp/test_shuffle_signal.fa"
		self.negSeqFile = "/projects/bhandare/workspace/PySG/src/resources/tmp/test_shuffle_nosignal.fa"

	def test_create_negative_seq_dict(self):
		print self.confFile;
		PosSeqDict = Kmer.CreateFastaWithSignal(self.confFile);
		NegSeqDict = NoSignal.ShuffleToCreateNoSignalSequences(PosSeqDict, self.confFile)
		self.assertEqual(len(PosSeqDict), len(NegSeqDict))
		self.assertTrue(os.path.isfile(self.posSeqFile))
		self.assertTrue(os.path.isfile(self.negSeqFile))

	def test_compute_dinucleotide_distribution(self):
		posSeqs, gc_list, fg_lengths = shuffle_utils.get_seqs(self.posSeqFile)
		dinuc_distrib = shuffle_utils.compute_dinuc_distrib(posSeqs)
		shuffle_utils.print_dinuc_distrib(dinuc_distrib, "dinuc.txt");


