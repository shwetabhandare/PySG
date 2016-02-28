import numpy as np
from SeqGenUtils import *
np.set_printoptions(precision=2)

def generate_pwm(uniform_dist=[0.25, 0.25, 0.25, 0.25]):
	scale_factor = 1;
	alpha = scale_factor * np.array(uniform_dist);	
	N = 10;
	samples = np.random.dirichlet(alpha, N);

def stats(scale_factor, utr_dist=[.27, .21, .23, .29], N=25):
	alpha = scale_factor * np.array(utr_dist);
	samples = np.random.dirichlet(alpha, N);
	print "                          alpha:", scale_factor
	print "              per element mean:", samples.mean(axis=0)
	print "per element standard deviation:", samples.std(axis=0)
	print 
	#print samples
	return samples;

for scale in [0.1, 1, 10, 100, 1000]:
	samples = stats(scale)
	seqLen = 25;
	keys = ['A', 'C', 'T', 'G']

	for sample in samples:
		b = list()
		seq = ""
		for count in range(seqLen):
			b = [(keys[0], sample[0]*100), (keys[1], sample[1]*100),  
			     (keys[2], sample[2]*100),  (keys[3], sample[3]*100)]
			seq += GetRandomNucleotide(b)
		print "Sequence generated: ", seq;


#generate_pwm();
