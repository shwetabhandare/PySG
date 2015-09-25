import fasta
import sys
import random

PositiveFileName = sys.argv[1]
NegativeFileName = sys.argv[2]
OutputFileName = sys.argv[3];
#OutputFileName = PositiveFileName.lstrip("TTP_PARCLIP_ConversionSpecificity")
#OutputFileName = OutputFileName.rstrip(".txt")
#OutputFileName = "Non_TTPPARCLIP" + OutputFileName + ".txt"
PositiveFile = open(PositiveFileName, "r")
NegativeFile = open(NegativeFileName, "r")
OutputFile = open(OutputFileName , "w")

PosLengths = []
for record in fasta.fasta_itr(PositiveFileName):						#Remember lengths of positive sequences
	seqLength = len(record.sequence)
	PosLengths.append(seqLength)

NegSequences = []
NegHeaders = []
for record in fasta.fasta_itr(NegativeFileName):
	sequence = record.sequence
	header = record.header
	
	NegSequences.append(sequence)
	NegHeaders.append(header)


for number in PosLengths:
	gotLength = False
	seqLength = 0
	index = 0
	allNucs = ""
	while(gotLength != True):
		index = random.randint(0, len(NegSequences))
		#print "Index is: ",index," Len NegSequences",len(NegSequences);
		if index == len(NegSequences):
			index = index - 1;
		allNucs = NegSequences[index]
		seqLength = len(allNucs)
		if seqLength > number:
			gotLength = True
	startPoint = random.randint(0, seqLength - (number+1))
	counter = 0
	selectedNucs = ""
	while (counter < number):
		selectedNucs = selectedNucs + allNucs[startPoint + counter]
		counter = counter + 1

	OutputFile.write(">")
	OutputFile.write(NegHeaders[index])
	OutputFile.write("\n")
	OutputFile.write(selectedNucs)
	OutputFile.write("\n")

PositiveFile.close()
NegativeFile.close()
OutputFile.close()
