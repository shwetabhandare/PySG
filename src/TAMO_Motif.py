import sys, os, math, errno, getopt, re
import TAMO
from   TAMO		import MotifTools
from   TAMO.seq import Fasta
import parseDreme

def GetKmerFromTextMotifForSeq(motifText, seq):
	motif = MotifTools.Motif_from_text(motifText);
	return motif.bestscanseq(seq);

def Read_Dreme_PSSM(lines):
	pwm = []
	name = "Dreme Motif";

	vals = []
	for line in lines.split('\n'):
		for item in line.split():
			vals.append(float(item))
		pwm.append(vals)
		vals = [];
	#print pwm

	m = MotifTools.toDict(pwm)
	motif = MotifTools.Motif_from_counts(m)
	return motif;


def GetKmerFromMotifFromPWM(pwm, seq):
	m = MotifTools.toDict(pwm)
	#print m
	motif = MotifTools.Motif_from_counts(m)
	return motif.bestscanseq(seq);

def GetKmerFromMotifFromPSSM(lines, seq):
	motif = Read_Dreme_PSSM(lines);
	return motif.bestscanseq(seq);

def Read_PWM(filename):
	pwm = []
	f = open(filename)
	lines = f.readlines()
	#print lines;
	fields = lines[0].split()
	#print "fields:", fields
	name = fields[1]
	
	#print lines[1]

	vals_A = lines[1].split('\t')
	vals_C = lines[2].split('\t')
	vals_G = lines[3].split('\t')
	vals_T = lines[4].split('\t')

	#print vals_A;

	for i in range(1,len(vals_A)):
		vals = [float(vals_A[i]), float(vals_C[i]), float(vals_G[i]), float(vals_T[i])]
		pwm.append(vals)
	
	f.close()
	
	# print "PWM:"
	# for v in pwm:
	# 	print "\t",v
		
	return name, pwm
#end Read_PWM
	
	
##########################################################################################

def Make_PWM_Motif(filename, motifBackGround=""):

	#print "# Reading PWM from: [%s]"%filename
	name, pwm = Read_PWM(filename)
	
	m = MotifTools.toDict(pwm)
	#print m
	motif = MotifTools.Motif_from_ll(m);
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
	dremeFile = sys.argv[1]
	pssmList = parseDreme.getPSSMListFromDremeFile(dremeFile);
	for pssmLines in pssmList:
		Read_Dreme_PSSM(pssmLines);

if __name__ == "__main__":
	import sys
	pwmFile = sys.argv[1]
	motif = Make_PWM_Motif(pwmFile)
	print motif;


