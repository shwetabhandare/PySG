import SeqGenUtils;
from   TAMO		import MotifTools
import TAMO_Motif;


def getMotifList(motifFile):
	lines = [line.rstrip() for line in open(motifFile)]
	return lines;


def findMotifInSequences(motifList, seqDict):
	matches = 0;
	motifMatchDict = dict();
	motifCountDict = dict();
	totalMisMatches = 0;

	for header, sequence in seqDict.iteritems():
		motifFound = False;

		for motif in motifList:
			m = TAMO_Motif.Make_Text_Motif(motif)
			match = m.bestscan(sequence)
			if match > m.maxscore * 0.7:
				motifFound = True;
				if header in motifMatchDict.keys():
					motifMatchDict[header].append(motif)
				else:
					motifMatchDict[header] = [motif]
				#print "header: ", header, ", match: ", match;

				if motif in motifCountDict.keys():
					motifCountDict[motif] = motifCountDict[motif] + 1
				else:
					motifCountDict[motif] = 1;

		if motifFound == False:
			totalMisMatches = totalMisMatches + 1;

	print "Total Mismatches: ", totalMisMatches;
	return motifMatchDict, motifCountDict;

if __name__ == "__main__":	
	import sys
	motifFile = sys.argv[1]
	faFile = sys.argv[2]

	motifList = getMotifList(motifFile)
	print motifList
	seqDict  = SeqGenUtils.fasta_read(faFile)
	motifMatchDict, motifCountDict = findMotifInSequences(motifList, seqDict)

	singleMotif = 0
	multipleMotif = 0;
	for header, value in motifMatchDict.iteritems():
		if len(value) > 1:
			multipleMotif = multipleMotif + 1;
		if len(value) == 1:
			singleMotif = singleMotif + 1;

	print "Sequences that matches single motif: ", singleMotif;
	print "Sequences that matches multiple motifs: ", multipleMotif;
	print "Total Sequences: ", len(seqDict)

	# for item in sorted(motifMatchDict):
	# 	header = item[0]
	# 	motif = item[1]
	# 	print header, motif, motifMatchDict[item]