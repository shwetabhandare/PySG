import unittest
import generateYaml
import glob
import os, fnmatch
import yaml
import SeqGenUtils
import parseResults
import uuid;
import Kmer;

class TestGenerateYaml(unittest.TestCase):
	def test_ConstructObject_experiment(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments.yaml"
		uuid_to_append = str(uuid.uuid4())

		generator = generateYaml.YamlFastaGenerator(confFile, uuid_to_append);
		targetDir = generator.GetTargetDir();

		self.assertTrue(os.path.isdir(targetDir))
		self.assertEqual(generator.GetNumSeqList(), [100, 200])
		self.assertEqual(generator.GetNumSeqLenList(), [50, 100])
		self.assertEqual(generator.GetPwmFilesList(), ["/projects/bhandare/workspace/PySG/data/pwm/Nhp6b.pwm"])
		self.assertEqual(generator.GetNoSignalType(), "dirichlet");
		self.assertEqual(generator.GetAlpha(), [100, 1000]);
		self.assertEqual(generator.GetSignalType(), "pwmFiles");
		self.assertEqual(generator.GetSignalPercentList(), [75, 90]);

	def test_ConstructObject_structure(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/test_structure.yml"
		uuid_to_append = str(uuid.uuid4())

		generator = generateYaml.YamlFastaGenerator(confFile, uuid_to_append);
		targetDir = generator.GetTargetDir();

		self.assertTrue(os.path.isdir(targetDir))
		generator.CreateConfFiles();

	def test_DictFromFasta(self):
		structureAlignmentFile = "/projects/bhandare/workspace/PySG/src/resources/RF000037.fa"
		structureDict = SeqGenUtils.fasta_read(structureAlignmentFile);
		self.assertEqual(len(structureDict), 62);

		for key, value in structureDict.iteritems():
			condition = 'T' in value;
			self.assertFalse(condition);

		structureDict = SeqGenUtils.ChangeUsToTs(structureDict);
		self.assertEqual(len(structureDict), 62);
		for key, value in structureDict.iteritems():
			condition = 'T' in value;
			self.assertTrue(condition);

	def test_KmersFromStructureFile(self):
		structureAlignmentFile = "/projects/bhandare/workspace/PySG/src/resources/RF000037.fa"
		kmers = Kmer.GetKmersFromStructureFile(structureAlignmentFile, 80);
		self.assertEqual(len(kmers), 80);





	# def test_ConstructObject_experiment1(self):
	# 	confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments1.yaml"
	# 	generator = generateYaml.YamlFastaGenerator(confFile);
	# 	targetDir = generator.GetTargetDir();

	# 	self.assertEqual(targetDir, "/projects/bhandare/workspace/PySG/src/resources/tmp")
	# 	self.assertTrue(os.path.isdir(targetDir))
	# 	self.assertEqual(generator.GetNumSeqList(), [100])
	# 	self.assertEqual(generator.GetNumSeqLenList(), [50])
	# 	self.assertEqual(generator.GetPwmDir(), '/projects/bhandare/workspace/PySG/data/pwm')
	# 	self.assertEqual(generator.GetNoSignalType(), "shuffle");
	# 	self.assertEqual(generator.GetAlpha(), []);
	# 	self.assertEqual(generator.GetSignalType(), "pwmDir");
	# 	self.assertEqual(generator.GetSignalPercentList(), [75, 90]);
	# 	self.assertEqual(generator.GetNoSignalFastaFile(), "/projects/bhandare/workspace/scripts/NegFileCreator/3UTR_transcripts_Human.txt");

	# def test_CreateConfFilesExpt(self):
	# 	confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments.yaml"
	# 	generator = generateYaml.YamlFastaGenerator(confFile);
	# 	targetDir = generator.GetTargetDir();
	# 	generator.CreateConfFiles();
	# 	self.assertTrue(os.path.isdir(targetDir))
	# 	numFiles = 0;
	# 	for root, dirs, files in os.walk(targetDir):
	# 		for file in fnmatch.filter(files, "*.yml"):
	# 			numFiles = numFiles + 1;
	# 	self.assertEqual(numFiles, 4);

	# def test_CreateConfFilesExpt1(self):
	# 	confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments1.yaml"
	# 	generator = generateYaml.YamlFastaGenerator(confFile);
	# 	targetDir = generator.GetTargetDir();
	# 	generator.CreateConfFiles();
	# 	self.assertTrue(os.path.isdir(targetDir))
	# 	numFiles = 0;
	# 	for root, dirs, files in os.walk(targetDir):
	# 		for file in fnmatch.filter(files, "*.yml"):
	# 			numFiles = numFiles + 1;
	# 	print "Num Files: ", str(numFiles)
	# 	self.assertEqual(numFiles, 14);

	# def test_CreateConfFilesExpt2(self):
	# 	confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments2.yaml"
	# 	generator = generateYaml.YamlFastaGenerator(confFile);
	# 	targetDir = generator.GetTargetDir();
	# 	generator.CreateConfFiles();
	# 	self.assertTrue(os.path.isdir(targetDir))
	# 	numFiles = 0;
	# 	for root, dirs, files in os.walk(targetDir):
	# 		for file in fnmatch.filter(files, "*.yml"):
	# 			numFiles = numFiles + 1;
	# 	print "Num Files: ", str(numFiles)
	# 	self.assertEqual(numFiles, 8);


	# def test_CreateConfFilesExpt3(self):
	# 	confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments3.yaml"
	# 	generator = generateYaml.YamlFastaGenerator(confFile);
	# 	targetDir = generator.GetTargetDir();
	# 	generator.CreateConfFiles();
	# 	self.assertTrue(os.path.isdir(targetDir))
	# 	numFiles = 0;
	# 	for root, dirs, files in os.walk(targetDir):
	# 		for file in fnmatch.filter(files, "*.yml"):
	# 			numFiles = numFiles + 1;
	# 	print "Num Files: ", str(numFiles)

		#self.assertEqual(numFiles, 8);