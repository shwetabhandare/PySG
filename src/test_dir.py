import numpy as np
import sys

from SeqGenUtils import *
np.set_printoptions(precision=2)

def generatePWM(scale_factor, N=10, uniform_dist=[0.25, 0.25, 0.25, 0.25]):
	alpha = scale_factor * np.array(uniform_dist);	
	samples = np.random.dirichlet(alpha, N);
	return np.transpose(samples)

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

	f.close();

if __name__ == "__main__":
	import sys
	scale_factor = int(sys.argv[1])
	pwmLength = int(sys.argv[2])
	pwmFile = sys.argv[3]
	pwm_array = generatePWM(scale_factor, pwmLength);
	writePwmFile(pwm_array, pwmFile, scale_factor)
