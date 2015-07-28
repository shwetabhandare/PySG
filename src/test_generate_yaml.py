import generateYaml
import glob
import os

class TestGenerateYaml:
	def test_conf(self):
		generateYaml.CreateConfFiles("/tmp");
		files = glob.glob("/tmp/*.yml")
		assert len(files) > 0

	def test_generate_yaml(self):
		generateYaml.generateYaml("100", "/tmp", 10, 100, 10, 20, 30, 10, 40, "HuR", "ATTTA", 1, 0);
		expectedContent = "{motif: {distance: 0, length: 5, location: random, motif: ATTTA, numMotifs: 1, type: HuR}, sequence: {A: 20, C: 40, G: 10, T: 30, maxLen: 100, minLen: 10, numSeq: 10}}"
		assert os.path.isfile('/tmp/seq_100.yml');
		content = ""
		with open('/tmp/seq_100.yml', 'r') as content_file:
			content = " ".join(line.strip() for line in content_file);
    		assert expectedContent.strip() == content.strip();