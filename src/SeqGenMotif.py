import sys, os, math, errno, getopt, re
import TAMO
from   TAMO		import MotifTools
from   TAMO.seq import Fasta

def ReadPWMFile(filename):
	pwm = []
	f = open(filename)
	lines = f.readlines()
	f.close()
	
	line1 = lines[1].rstrip().split('\t')
	line2 = lines[2].rstrip().split('\t')
	line3 = lines[3].rstrip().split('\t')
	line4 = lines[4].rstrip().split('\t')

	print line1, line2, line3, line4
	pwm.append(line1)
	pwm.append(line2)
	pwm.append(line3)
	pwm.append(line4)

	
	print "PWM:"
	for v in pwm:
		print "\t",v
		
	name = "dreme";
	return name, pwm

def MakePWMMotif(filename):

	print "# Reading PWM from: [%s]"%filename
	name, pwm = ReadPWMFile(filename)
	
	print "Building motif:", name
	m = MotifTools.toDict(pwm)
	print m
	#motif = MotifTools.Motif_from_ll(m)
	#A 0.307 C 0.141 G 0.060 T 0.493
	#motif = MotifTools.Motif_from_counts(m, beta=0.01, bg={'A': 0.307, 'C': 0.141, 'G': 0.060, 'T': 0.493})
	motif = MotifTools.Motif_from_counts(m)
	motif.source = name
	
	print "Motif:", motif.source
	print "Max Motif Score:", motif.maxscore
	print "Motif Summary:", motif.summary()
	motif.printlogo(2.3,10)
	
	return motif

if __name__ == "__main__":
	import sys
	motifFile = sys.argv[1]
	print "Motif file: ", motifFile
	name, pwm = MakePWMMotif(motifFile)
	#name, pwn = ReadPWMFile(motifFile)
