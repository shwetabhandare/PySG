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
import generateGraphs;

class TestShuffleSeq(unittest.TestCase):
	conFile = "";
	posSeqFile = ""
	negSeqFile = ""

	def setUp(self):
		self.confFile = "/projects/bhandare/workspace/PySG/src/resources/test_shuffle.yml"
		self.posSeqFile = "/projects/bhandare/workspace/PySG/src/resources/tmp/test_shuffle_signal.fa"
		self.negSeqFile = "/projects/bhandare/workspace/PySG/src/resources/tmp/test_shuffle_nosignal.fa"

	def test_create_negative_seq_dict(self):
		#print self.confFile;
		PosSeqDict = Kmer.CreateFastaWithSignal(self.confFile);
		NegSeqDict = NoSignal.ShuffleToCreateNoSignalSequences(PosSeqDict, self.confFile)
		self.assertEqual(len(PosSeqDict), len(NegSeqDict))
		self.assertTrue(os.path.isfile(self.posSeqFile))
		self.assertTrue(os.path.isfile(self.negSeqFile))

	def test_compute_dinucleotide_distribution(self):
		PosSeqDict = Kmer.CreateFastaWithSignal(self.confFile);
		NegSeqDict = NoSignal.ShuffleToCreateNoSignalSequences(PosSeqDict, self.confFile)

		posSeqs, gc_list, fg_lengths = shuffle_utils.get_seqs(self.posSeqFile)
		pos_dinuc_distrib = shuffle_utils.compute_dinuc_distrib(posSeqs, True)
		print "Positive Set: \n", pos_dinuc_distrib;

		negSeqs, gc_list, fg_lengths = shuffle_utils.get_seqs(self.negSeqFile)
		neg_dinuc_distrib = shuffle_utils.compute_dinuc_distrib(negSeqs, True)
		print "Shuffled: Negative Set: \n", neg_dinuc_distrib;

	def test_compute_dinucleotide_distribution_random_files(self):

		posSeqFile = "/projects/bhandare/workspace/PySG/src/resources/tmp/af36a99f-d018-4e2c-9a4c-db59c53b5ff3/100_100_100_100/Signal_100_100_100_100_RF000037.fa.fa"
		negSeqFile = "/projects/bhandare/workspace/PySG/src/resources/tmp/af36a99f-d018-4e2c-9a4c-db59c53b5ff3/100_100_100_100/NoSignal_100_100_100_100_RF000037.fa.fa"

		#shuffle_utils.print_dinuc_distrib(dinuc_distrib, "dinuc.txt");
		posSeqs, gc_list, fg_lengths = shuffle_utils.get_seqs(posSeqFile)
		pos_dinuc_distrib = shuffle_utils.compute_dinuc_distrib(posSeqs, True)
		print "Dirichlet: Positive Set: \n", pos_dinuc_distrib;

		negSeqs, gc_list, fg_lengths = shuffle_utils.get_seqs(negSeqFile)
		neg_dinuc_distrib = shuffle_utils.compute_dinuc_distrib(negSeqs, True)
		print "Dirichlet: Negative Set: \n", neg_dinuc_distrib;


	def test_compute_dinucleotide_distribution_total(self):
		resultDir = "/projects/bhandare/workspace/PySG/src/resources/tmp"
		generateGraphs.computeDiNucleotideDistribution(resultDir)

	def test_compute_dinucleotide_distribution_3utrs(self):
		three_utr_file_dir = "/projects/bhandare/workspace/scripts/NegFileCreator/";
		three_utr_file = "3UTR_transcripts_Human.txt"

		threeUtrDist = generateGraphs.GetDistributionForFiles(three_utr_file_dir, three_utr_file);
		#graphFileName = "/projects/bhandare/workspace/PySG/src/tmp_results/3utr.png"
		#generateGraphs.PlotGraphs(threeUtrDist, threeUtrDist, graphFileName);
