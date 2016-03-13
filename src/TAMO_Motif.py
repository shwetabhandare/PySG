import sys, os, math, errno, getopt, re
import TAMO
from   TAMO		import MotifTools
from   TAMO.seq import Fasta

def GetKmerForSeq(motifText, seq):
	motif = MotifTools.Motif_from_text(motifText);
	return motif.bestscanseq(seq);

def Read_Dreme_PSSM(filename):
	pwm = []
	f = open(filename)
	lines = f.readlines()
	fields = lines[0].split()
	
	name = fields[1]

	vals = []
	for item in lines[1].split():
		vals.append(float(item))
	pwm.append(vals)
	vals = []
	for item in lines[2].split():
		vals.append(float(item))
	pwm.append(vals)
	vals = []
	for item in lines[3].split():
		vals.append(float(item))
	pwm.append(vals)
	vals = []
	for item in lines[4].split():
		vals.append(float(item))
	pwm.append(vals)
	vals = []
	for item in lines[5].split():
		vals.append(float(item))
	pwm.append(vals)
	vals = []
	for item in lines[6].split():
		vals.append(float(item))
	pwm.append(vals)

	#print pwm

	m = MotifTools.toDict(pwm)
	motif = MotifTools.Motif_from_text('MCCCGA')
	#motif = MotifTools.Motif_from_counts(m)
	#print m;

	seq1 = "GGGGGCAGCGTCGGGTTTTTTTTTTCCGAGGAACGCGGCATATCAAGAGGTATAAGGTCTATCTAACCTTACCCGCTATAAGCAATAGCCAAAAAACGAC"

	seq2 = "GGAACCGCGTTCGGGGGGGGGGGGGCCGAACCCTTCCAGCATTGAGCTCCTGCCGCTAGCTTATGCGGCCTCCCATCCAGTCGGCCGAGACGCACGACTT"

	seq3 = "CTAGCCACGATCGGGTTTTTTTTTTCCGACTTTCACTCCGCATAGTTCGCACTGACCCAGGTGGTCCTAAACACTCGCATCGGTATCCCGTCCTAGTCTA"
	print "Seq 1: ", motif.bestscanseq(seq1)
	print "Seq 2: ", motif.bestscanseq(seq2)
	print "Seq 3: ", motif.bestscanseq(seq3)

	print "Best k-mers: "
	print motif.bestseqs()
	return pwm;

def Read_PWM(filename):
	pwm = []
	f = open(filename)
	lines = f.readlines()
	fields = lines[0].split()
	#print "fields:", fields
	name = fields[1]
	
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
	#	print "\t",v
		
	return name, pwm
#end Read_PWM
	
	
##########################################################################################

def Make_PWM_Motif(filename, motifBackGround=""):

	#print "# Reading PWM from: [%s]"%filename
	name, pwm = Read_PWM(filename)
	
	#print "Building motif:", name
	m = MotifTools.toDict(pwm)
	print m
	#print "Motif BackGround: ", motifBackGround
	if motifBackGround != "":
		motif = MotifTools.Motif_from_counts(m, bg=motifBackGround)
	else:
		motif = MotifTools.Motif_from_counts(m)
	motif.source = name
	
	#print "Motif:", motif.source
	#print "Max Motif Score:", motif.maxscore
	#print "Motif Summary:", motif.summary()
	#motif.printlogo(2.3,10)
	
	return motif
#end Make_PWM_Motif

##########################################################################################

def Make_Text_Motif(textMotif):
	return MotifTools.Motif_from_text(textMotif)

if __name__ == "__main__":
	motiffile = sys.argv[1]
	dreme = int(sys.argv[2])
	if (dreme == 1):
		Read_Dreme_PSSM(motiffile)
	else:
		motif = Make_PWM_Motif(motiffile)
		print motif;
		print motif.random_kmer()
	
