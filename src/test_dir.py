import numpy as np
import sys

from SeqGenUtils import *
np.set_printoptions(precision=2)

def generate_pwm(scale_factor, uniform_dist=[0.25, 0.25, 0.25, 0.25]):
	alpha = scale_factor * np.array(uniform_dist);	
	N = 10;
	samples = np.random.dirichlet(alpha, N);
	return np.transpose(samples)

def stats(scale_factor, utr_dist=[.27, .21, .23, .29], N=25):
	alpha = scale_factor * np.array(utr_dist);
	samples = np.random.dirichlet(alpha, N);
	print "                          alpha:", scale_factor
	print "              per element mean:", samples.mean(axis=0)
	print "per element standard deviation:", samples.std(axis=0)
	print 
	#print samples
	return samples;

def writePwmFile(pwm_array, pwmFile, scale_factor):
	f = open(pwmFile, 'w');

	titleLine = "Gene: Scale Factor: " + str(scale_factor) + "\n"
	f.write(titleLine);
	lineCount = 0;
	for pwmLine in pwm_array:
		if lineCount == 0:
			lineToAdd = "A:"
		elif lineCount == 1:
			lineToAdd = "C:"
		elif lineCount == 2:
			lineToAdd = "G:"
		elif lineCount == 3:
			lineToAdd = "T:"
		else:
			print "Error : Invalid line count", str(lineCount)
		for value in pwmLine:
			lineToAdd = lineToAdd + "\t" + '%.2f'%(value);
		lineToAdd = lineToAdd + "\n";

		f.write(lineToAdd);
		lineCount = lineCount + 1;
	print type(pwm_array);




# for scale in [0.1, 1, 10, 100, 1000]:
# 	samples = stats(scale)
# 	seqLen = 25;
# 	keys = ['A', 'C', 'T', 'G']

# 	for sample in samples:
# 		b = list()
# 		seq = ""
# 		for count in range(seqLen):
# 			b = [(keys[0], sample[0]*100), (keys[1], sample[1]*100),  
# 			     (keys[2], sample[2]*100),  (keys[3], sample[3]*100)]
# 			seq += GetRandomNucleotide(b)
# 		print "Sequence generated: ", seq;



scale_factor = int(sys.argv[1])
pwmFile = sys.argv[2]
pwm_array = generate_pwm(scale_factor);
writePwmFile(pwm_array, pwmFile, scale_factor)
