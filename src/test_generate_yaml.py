import unittest
import generateYaml
import glob
import os
import yaml
import SeqGenUtils

class TestGenerateYaml(unittest.TestCase):
	def test_getTargetDir(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments.yaml"
		confMap = SeqGenUtils.GetConf(confFile);
		targetDir = generateYaml.getTargetDir(confMap);
		self.assertEqual(targetDir, "/projects/bhandare/workspace/PySG/src/resources/tmp")
		self.assertTrue(os.path.isdir(targetDir))

	def test_setupVariables(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments.yaml"
		confMap = SeqGenUtils.GetConf(confFile);
		targetDir = generateYaml.getTargetDir(confMap);	
		signalType, noSignalType = generateYaml.setupVariables(confMap)
		self.assertEqual(signalType, "pwmFiles");
		self.assertEqual(noSignalType, "dirichlet");

	def test_setupVariables_PwmFileDir(self):
		confFile = "/projects/bhandare/workspace/PySG/src/resources/experiments1.yaml"
		confMap = SeqGenUtils.GetConf(confFile);
		targetDir = generateYaml.getTargetDir(confMap);	
		signalType, noSignalType = generateYaml.setupVariables(confMap)
		self.assertEqual(signalType, "pwmDir");
		self.assertEqual(noSignalType, "shuffle");