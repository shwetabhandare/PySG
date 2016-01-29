import sys, os, math, errno, getopt, re
import TAMO
from   TAMO		import MotifTools
from   TAMO.seq import Fasta

def ReadPWMFile(filename):
	pwm = []
	f = open(filename)
	lines = f.readlines()
	fields = lines[0].split()
	print "fields:", fields[4]
	#print "fields:", fields
	name = fields[4]
	print "Name: ", name
	
	vals_A = lines[1].rstrip().split('\t')
	vals_C = lines[2].rstrip().split('\t')
	vals_T = lines[3].rstrip().split('\t')
	vals_G = lines[4].rstrip().split('\t')

	for i in range(0,len(vals_A)):
		vals = [float(vals_A[i]), float(vals_C[i]), float(vals_G[i]), float(vals_T[i])]
		pwm.append(vals)
	
	f.close()
	
	#print "PWM:"
	#for v in pwm:
	#	print "\t",v
		
	return name, pwm

def MakePWMMotif(filename):

	print "# Reading PWM from: [%s]"%filename
	name, pwm = ReadPWMFile(filename)
	
	print "Building motif:", name
	#m = MotifTools.toDict(pwm)
	#print m
	#motif = MotifTools.Motif_from_counts(m)
	#print vals_A, vals_C, vals_T, vals_G
	motif = MotifTools.Motif_from_text("TGASTCA")
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