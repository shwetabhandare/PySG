import sys, os, math, errno, getopt, re
import TAMO
from   TAMO		import MotifTools
from   TAMO.seq import Fasta

def Read_PWM(filename):
	pwm = []
	f = open(filename)
	lines = f.readlines()
	fields = lines[0].split()
	#print "fields:", fields
	name = fields[4]
	print "Motif Name: ", name;
	
	vals_A = lines[1].split('\t')
	vals_C = lines[2].split('\t')
	vals_G = lines[3].split('\t')
	vals_T = lines[4].split('\t')

	for i in range(1,len(vals_A)):
		vals = [float(vals_A[i]), float(vals_C[i]), float(vals_G[i]), float(vals_T[i])]
		pwm.append(vals)
	
	f.close()
	
	#print "PWM:"
	#for v in pwm:
		#print "\t",v
		
	return name, pwm
#end Read_PWM
	
	
##########################################################################################

def Make_PWM_Motif(filename):

	print "# Reading PWM from: [%s]"%filename
	name, pwm = Read_PWM(filename)
	
	print "Building motif:", name
	m = MotifTools.toDict(pwm)
	motif = MotifTools.Motif_from_counts(m)
	motif.source = name
	
	print "Motif:", motif.source
	print "Max Motif Score:", motif.maxscore
	print "Motif Summary:", motif.summary()
	motif.printlogo(2.3,10)
	
	return motif
#end Make_PWM_Motif

##########################################################################################

if __name__ == "__main__":
	motiffile = sys.argv[1]
	motif = Make_PWM_Motif(motiffile)
	print motif.random_kmer()
