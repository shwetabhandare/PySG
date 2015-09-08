import generateYaml
import glob
import os

class TestGenerateYaml:
	def test_conf(self):
		generateYaml.CreateConfFiles("/tmp");
		files = glob.glob("/tmp/*.yml")
		assert len(files) > 0

	def test_generate_yaml_hur(self):
		generateYaml.generateYaml("100", "/tmp", 10, 100, 10, 20, 30, 10, 40, "HuR", "ATTTA", 1, 0);
		expectedContent = "{motif: {distance: 0, length: 5, location: random, motif: ATTTA, numMotifs: 1, type: HuR}, sequence: {A: 20, C: 40, G: 10, T: 30, maxLen: 100, minLen: 10, numSeq: 10}}"
		assert os.path.isfile('/tmp/seq_100.yml');
		content = ""
		with open('/tmp/seq_100.yml', 'r') as content_file:
			content = " ".join(line.strip() for line in content_file);
    		assert expectedContent.strip() == content.strip();

	def test_generate_yaml_ttp(self):
		generateYaml.generateYaml("200", "/tmp", 10, 100, 10, 20, 30, 10, 40, "TTP", "TTATTTATT", 1, 0);
		expectedContent = "{motif: {distance: 0, length: 9, location: random, motif: TTATTTATT, numMotifs: 1, type: TTP}, sequence: {A: 20, C: 40, G: 10, T: 30, maxLen: 100, minLen: 10, numSeq: 10}}"
		assert os.path.isfile('/tmp/seq_200.yml');
		content = ""
		with open('/tmp/seq_200.yml', 'r') as content_file:
			content = " ".join(line.strip() for line in content_file);
    		assert expectedContent.strip() == content.strip();   	

 	def test_generate_yaml_generated(self):
		generateYaml.generateYaml("300", "/tmp", 10, 100, 10, 20, 30, 10, 40, "Generated", "AAATTTGGGCCC", 1, 0);
		expectedContent = "{motif: {distance: 0, length: 12, location: random, motif: AAATTTGGGCCC, numMotifs: 1, type: Generated}, sequence: {A: 20, C: 40, G: 10, T: 30, maxLen: 100, minLen: 10, numSeq: 10}}"
		assert os.path.isfile('/tmp/seq_300.yml');
		content = ""
		with open('/tmp/seq_300.yml', 'r') as content_file:
			content = " ".join(line.strip() for line in content_file);
    		assert expectedContent.strip() == content.strip();   	s