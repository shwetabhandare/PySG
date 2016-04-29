import unittest
import generateYaml
import glob
import os, fnmatch
import yaml
import SeqGenUtils

class TestGenerateYaml(unittest.TestCase):
	def test_ConstructObject_experiment(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments.yaml"
		generator = generateYaml.YamlFastaGenerator(confFile);

		targetDir = generator.GetTargetDir();

		self.assertEqual(targetDir, "/projects/bhandare/workspace/PySG/src/resources/tmp")
		self.assertTrue(os.path.isdir(targetDir))
		self.assertEqual(generator.GetNumSeqList(), [10])
		self.assertEqual(generator.GetNumSeqLenList(), [50])
		self.assertEqual(generator.GetPwmFilesList(), ["/projects/bhandare/workspace/PySG/data/pwm/Nhp6b.pwm"])
		self.assertEqual(generator.GetNoSignalType(), "dirichlet");
		self.assertEqual(generator.GetAlpha(), [100, 1000]);
		self.assertEqual(generator.GetSignalType(), "pwmFiles");
		self.assertEqual(generator.GetSignalPercentList(), [75, 90]);

	def test_ConstructObject_experiment1(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments1.yaml"
		generator = generateYaml.YamlFastaGenerator(confFile);
		targetDir = generator.GetTargetDir();

		self.assertEqual(targetDir, "/projects/bhandare/workspace/PySG/src/resources/tmp")
		self.assertTrue(os.path.isdir(targetDir))
		self.assertEqual(generator.GetNumSeqList(), [10])
		self.assertEqual(generator.GetNumSeqLenList(), [50])
		self.assertEqual(generator.GetPwmDir(), '/projects/bhandare/workspace/PySG/data/pwm')
		self.assertEqual(generator.GetNoSignalType(), "shuffle");
		self.assertEqual(generator.GetAlpha(), []);
		self.assertEqual(generator.GetSignalType(), "pwmDir");
		self.assertEqual(generator.GetSignalPercentList(), [75, 90]);
		self.assertEqual(generator.GetNoSignalFastaFile(), "/my/tmp/fasta.file");

	def test_CreateConfFilesExpt(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments.yaml"
		generator = generateYaml.YamlFastaGenerator(confFile);
		targetDir = generator.GetTargetDir();
		generator.CreateConfFiles();
		self.assertTrue(os.path.isdir(targetDir))
		numFiles = 0;
		for root, dirs, files in os.walk(targetDir):
			for file in fnmatch.filter(files, "*.yml"):
				numFiles = numFiles + 1;
		self.assertEqual(numFiles, 4);

	# def test_setupVariables_PwmFileDir(self):
	# 	confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments1.yaml"
	# 	confMap = SeqGenUtils.GetConf(confFile);
	# 	targetDir = generateYaml.getTargetDir(confMap);	
	# 	signalType, noSignalType = generateYaml.setupVariables(confMap)
	# 	self.assertEqual(signalType, "pwmDir");
	# 	self.assertEqual(noSignalType, "shuffle");

	# def test_createConfFiles_Dirichlet_PwmFile(self):
	# 	confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments.yaml"

	# 	confMap = SeqGenUtils.GetConf(confFile);
	# 	targetDir = generateYaml.getTargetDir(confMap);	
	# 	signalType, noSignalType = generateYaml.setupVariables(confMap)
	# 	self.assertEqual(signalType, "pwmFiles");
	# 	self.assertEqual(noSignalType, "dirichlet");

	# 	generateYaml.CreateConfFiles(targetDir, noSignalType);